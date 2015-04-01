
#include <zlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "vcf.h"
#include "memutil.h"


const char *vcf_fix_headers[] =
  {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
   "FILTER", "INFO", "FORMAT"};



void vcf_read_header(gzFile vcf_fh, VCFInfo *vcf_info) {
  char *line, *cur, *token;
  int tok_num;
  int n_fix_header;
  const char delim[] = " \t";

  n_fix_header = sizeof(vcf_fix_headers) / sizeof(const char *);

  fprintf(stderr, "there are %d fixed headers\n", n_fix_header);

  vcf_info->n_header_lines = 0;
  
  while(1) {
    line = util_gzgets_line(vcf_fh);

    if(line == NULL) {
      my_err("%s:%d: could not read header information from file",
	     __FILE__, __LINE__);

    }

    if(util_str_starts_with(line, "##")) {
      /* header line */
      vcf_info->n_header_lines += 1;
    }
    else if(util_str_starts_with(line, "#CHROM")) {
      /* this should be last header line that contains list of fixed fields */
      vcf_info->n_header_lines += 1;
	
      cur = line;
      tok_num = 0;
      while((token = strsep(&cur, delim)) != NULL) {
	if(tok_num < n_fix_header) {
	  if(strcmp(token, vcf_fix_headers[tok_num]) != 0) {
	    my_warn("expected token %d to be %s but got '%s'",
		    tok_num, vcf_fix_headers[tok_num], token);
	  }
	}
	tok_num += 1;
      }
      vcf_info->n_samples = tok_num - n_fix_header;
      fprintf(stderr, "there are %d samples\n", vcf_info->n_samples);
      break;
    } else {
      my_err("expected last line in header to start with #CHROM");
    }
    my_free(line);
  }
}




void parse_geno_probs(VCFInfo *vcf_info, double *geno_probs,
		      char *cur) {
  char delim[] = " \t";
  char inner_delim[] = ":";
  char *tok, *inner_tok, *inner_cur;
  char fmt[VCF_MAX_FORMAT];
  char gtype[VCF_MAX_FORMAT];
  long gl_idx, i, n, n_geno_probs, expect_geno_probs;
  double like_homo_ref, like_het, like_homo_alt;
  double prob_homo_ref, prob_het, prob_homo_alt, prob_sum;

  expect_geno_probs = vcf_info->n_samples * 3;
  
  /* copy format string before tokenizing (which modifies it) */
  strcpy(fmt, vcf_info->format);

  /* find index of GL in format specifier */
  gl_idx = -1;
  i = 0;
  inner_cur = fmt;
  while((tok = strsep(&inner_cur, inner_delim)) != NULL) {
    if(strcmp(tok, "GL") == 0) {
      /* found genotype likelihood format string */
      gl_idx = i;
      break;
    }
    i++;
  }

  if(gl_idx == -1) {
    my_err("%s:%d: VCF format string does not specify GL token, cannot "
	   "obtain genotype probabilities", __FILE__, __LINE__);
  }

  n_geno_probs = 0;
  
  while((tok = strsep(&cur, delim)) != NULL) {
    /* each genotype string is delimited by ':'
     * each GL portion is delimited by ','
     */
    util_strncpy(gtype, tok, sizeof(gtype));

    i = 0;
    inner_cur = gtype;
    while((inner_tok = strsep(&inner_cur, inner_delim)) != NULL) {
      if(i == gl_idx) {
	n = sscanf(inner_tok, "%lg,%lg,%lg", &like_homo_ref, &like_het,
		   &like_homo_alt);

	if(n != 3) {
	  if(strcmp(inner_tok, ".") == 0) {
	    /* '.' indicates missing data
	     * set all likelihoods to log(0.333) = -0.477
	     */
	    like_homo_ref = like_het = like_homo_alt = -0.477;
	  } else {
	    my_err("%s:%d: failed to parse genotype likelihoods from "
		   "string '%s'", __FILE__, __LINE__, inner_tok);
	  }
	}

	/* convert log10(prob) to prob */
	prob_homo_ref = pow(10.0, like_homo_ref);
	prob_het = pow(10.0, like_het);
	prob_homo_alt = pow(10.0, like_homo_alt);

	if((n_geno_probs + 3) > expect_geno_probs) {
	  my_err("%s:%d: more genotype likelihoods per line than expected",
		 __FILE__, __LINE__);
	}
	
	/* most of time probs sum to 1.0, but sometimes they do not
	 * possibly reflects different likelihoods used for indel 
	 * calling but not sure. Normalize probs so they sum to 1.0
	 * This is like getting posterior assuming uniform prior.
	 */
	prob_sum = prob_homo_ref + prob_het + prob_homo_alt;
	prob_homo_ref = prob_homo_ref / prob_sum;
	prob_het = prob_het / prob_sum;
	prob_homo_alt = prob_homo_alt / prob_sum;
     	
	geno_probs[n_geno_probs] = prob_homo_ref;
	geno_probs[n_geno_probs + 1] = prob_het;
	geno_probs[n_geno_probs + 2] = prob_homo_alt;

	n_geno_probs += 3;
      }
    }
  }

  if(n_geno_probs != expect_geno_probs) {
    my_err("%s:%d: expected %ld genotype likelihoods per line, but got "
	   "%ld", __FILE__, __LINE__, expect_geno_probs, n_geno_probs);
  }
  
}



/**
 * Gets next line of VCF file and parses it into VCFInfo datastructure.
 *
 * If geno_probs array is non-null genotype likelihoods are parsed and
 * stored in the provided array. The array must be of length
 * n_samples*3.
 *
 * Returns 0 on success, -1 if at EOF.
 */
int vcf_read_line(gzFile vcf_fh, VCFInfo *vcf_info,
		  double *geno_probs) {
  char *line, *cur, *token;
  int n_fix_header, ref_len, alt_len;
  size_t tok_num;
  const char delim[] = " \t";

  n_fix_header = sizeof(vcf_fix_headers) / sizeof(const char *);

  /* read a line */
  line = util_gzgets_line(vcf_fh);

  if(!line) {
    return -1;
  }
  
  cur = line;
  tok_num = 0;

  /* chrom */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->chrom, token, sizeof(vcf_info->chrom));
  
  
  /* pos */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  vcf_info->pos = util_parse_long(token);
  
  /* ID */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->id, token, sizeof(vcf_info->id));
  
  /* ref */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  vcf_info->ref_len = strlen(token);
  ref_len = util_strncpy(vcf_info->ref_allele, token,
			 sizeof(vcf_info->ref_allele));

  if(ref_len != vcf_info->ref_len) {
    my_warn("truncating long allele (%ld bp) to %ld bp\n",
	    vcf_info->ref_len, ref_len);
  }
  
  /* alt */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  vcf_info->alt_len = strlen(token);
  alt_len = util_strncpy(vcf_info->alt_allele, token,
			 sizeof(vcf_info->alt_allele));

  if(alt_len != vcf_info->alt_len) {
    my_warn("truncating long allele (%ld bp) to %ld bp\n",
	    vcf_info->alt_len, alt_len);
  }

  /* qual */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->qual, token, sizeof(vcf_info->qual));

  /* filter */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->filter, token, sizeof(vcf_info->filter));


  /* info */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->info, token, sizeof(vcf_info->info));

  
  /* format */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", n_fix_header);
  }
  util_strncpy(vcf_info->format, token, sizeof(vcf_info->format));


  /* now parse genotypes and or genotype likelihoods */
  if(geno_probs) {
    parse_geno_probs(vcf_info, geno_probs, cur);
  }
}
