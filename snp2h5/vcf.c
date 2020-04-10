
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "util.h"
#include "vcf.h"
#include "memutil.h"

#define VCF_GTYPE_MISSING -1



const char *vcf_fix_headers[] =
  {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
   "FILTER", "INFO", "FORMAT"};



/**
 * Allocates memory for VCFInfo structure
 * which is used for parsing VCF files.
 */
VCFInfo *vcf_info_new() {
  VCFInfo *vcf_info;

  vcf_info = my_malloc(sizeof(VCFInfo));

  /* dynamic buffer used for reading lines (can grow)  */
  vcf_info->buf_size = 1024;
  vcf_info->buf = my_malloc(vcf_info->buf_size);

  vcf_info->n_sample = 0;
  vcf_info->cur_line = 0;
  vcf_info->has_format = 0;
  vcf_info->sample_names = NULL;

  return vcf_info;
}


/** 
 * free memory allocated for reading lines
 */
void vcf_info_free(VCFInfo *vcf_info) {
  int i;
  
  if(vcf_info->sample_names) {
    for(i = 0; i < vcf_info->n_sample; i++) {
      my_free(vcf_info->sample_names[i]);
    }
    if(vcf_info->sample_names) {
      my_free(vcf_info->sample_names);
    }
  }

  
  my_free(vcf_info->buf);
  my_free(vcf_info);
}



/**
 * returns the name of the chromosome that is present on the
 * first non-header line of the VCF. If the file has no 
 * non-header lines, NULL is returned.
 * The returned string should be freed once it is no longer needed.
 */
char *vcf_get_chrom_name(const char *filename) {
  gzFile gzf;
  VCFInfo *vcf_info;
  SNP snp;
  int ret;
  char *chrom_name;


  gzf = util_must_gzopen(filename, "rb");

  vcf_info = vcf_info_new();

  /* read header */
  vcf_read_header(gzf, vcf_info);
  
  if(vcf_info->n_sample == 0) {
    return NULL;
  }
  
  /* read first line */
  ret = vcf_read_line(gzf, vcf_info, &snp, NULL, NULL, NULL);

  if(ret == -1) {
    /* at end of file */
    my_warn("VCF file %s contained no data lines\n", filename);
    return NULL;
  }

  chrom_name = util_str_dup(snp.chrom);
  
  vcf_info_free(vcf_info);
  
  gzclose(gzf);
  
  return chrom_name;
}

void vcf_read_header(gzFile vcf_fh, VCFInfo *vcf_info) {
  char *line, *cur, *token;
  int tok_num;
  int n_fix_header, i;
  
  /* const char delim[] = " \t"; */
  const char delim[] = "\t";

  n_fix_header = sizeof(vcf_fix_headers) / sizeof(const char *);

  vcf_info->n_header_line = 0;
  
  while(util_gzgetline(vcf_fh, &vcf_info->buf, &vcf_info->buf_size) != -1) {
    line = vcf_info->buf;
    vcf_info->cur_line += 1;

    if(line[0] == '\n' || line[0] == '\0') {
      /* empty line, assume part of header */
      vcf_info->n_header_line += 1;
    }
    else if(util_str_starts_with(line, "##")) {
      /* header line */
      vcf_info->n_header_line += 1;
    }
    else if(util_str_starts_with(line, "#CHROM")) {
      /* this should be last header line that contains list of fixed fields */
      vcf_info->n_header_line += 1;
	
      cur = vcf_info->buf;
      line = util_str_dup(vcf_info->buf);
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

      if(tok_num < n_fix_header) {
	my_warn("VCF contains fewer header columns in #CHROM line "
		"than expected\n");
	
	fprintf(stderr, "missing headers: ");
	for(i = tok_num; i < n_fix_header; i++) {
	  fprintf(stderr, "%s", vcf_fix_headers[i]);
	  (i < n_fix_header-1) ? fprintf(stderr, ",") : fprintf(stderr, "\n");
	}
	vcf_info->n_sample = 0;
      } else {
	vcf_info->has_format = TRUE;
	vcf_info->n_sample = tok_num - n_fix_header;
      }

      if(vcf_info->n_sample == 0) {
	my_warn("VCF header contains no sample identifiers in #CHROM line\n");
      } else {      
	/*
	 * read sample names from remaining part of header
	 */
	vcf_info->sample_names = my_malloc(sizeof(char *) * vcf_info->n_sample);
	cur = line;
	tok_num = 0;
	i = 0;
	while((token = strsep(&cur, delim)) != NULL) {
	  if(tok_num >= n_fix_header) {
	    vcf_info->sample_names[i] = util_str_dup(token);
	    i += 1;
	  }
	  tok_num += 1;
	}
      }
      my_free(line);
      break;
    } else {
      /* VCF did not appear to contain any header lines starting with '#' */
      my_warn("VCF file did not contain any header lines. "
	      "Expected final header line to start with '#CHROM' "
	      "and to contain sample names.\n");
      vcf_info->n_sample = 0;
      vcf_info->has_format = FALSE;
      break;
    }
  }
}



/**
 * Parses ':'-delimited format string and 
 * returns index of token that matches. 
 * Returns -1 on failure.
 * 
 */
int get_format_index(const char *format_str, const char *label) {
  char fmt[VCF_MAX_FORMAT];
  char delim[] = ":";
  char *cur, *tok;
  int idx, i;

  /* copy format string before tokenizing (which modifies it) */
  strcpy(fmt, format_str);

  /* find index of label format specifier */
  idx = -1;
  i = 0;
  cur = fmt;
  while((tok = strsep(&cur, delim)) != NULL) {
    if(strcmp(tok, label) == 0) {
      /* found label in format string */
      idx = i;
      break;
    }
    i++;
  }

  return idx;
}


void vcf_parse_haplotypes(VCFInfo *vcf_info, char *haplotypes,
			  char *cur, char *haplotypes_phase) {
  int gt_idx, hap1, hap2, i, n, phase;
  static int warn_phase = TRUE;
  static int warn_parse = TRUE;
  long expect_haps, n_haps, expect_phase, n_phase;
  char gt_str[VCF_MAX_FORMAT];
  
  /* char delim[] = " \t"; */
  char delim[] = "\t";
  char inner_delim[] = ":";
  char *inner_cur, *tok, *inner_tok;

  /* get index of GT token in format string*/
  gt_idx = get_format_index(vcf_info->format, "GT");
  if(gt_idx == -1) {
    my_err("%s:%d: VCF format string does not specify GT token "
	   "so cannot obtain haplotypes. Format string: '%s'.\n"
	   "To use this file, you must run snp2h5 without "
	   "the --haplotype option.",
	   __FILE__, __LINE__, vcf_info->format);
  }
  
  expect_haps = vcf_info->n_sample * 2;
  expect_phase = vcf_info->n_sample;
  
  n_haps = 0;
  n_phase = 0;
  
  while((tok = strsep(&cur, delim)) != NULL) {
    /* Each genotype string is delimited by ':'
     * The GT portions of the string are delimited by '/' or '|'
     * '|' indicates phased, '/' indicates unphased.
     */
    util_strncpy(gt_str, tok, sizeof(gt_str));
    
    i = 0;
    inner_cur = gt_str;
    while((i <= gt_idx) && (inner_tok = strsep(&inner_cur, inner_delim)) != NULL) {
      if(i == gt_idx) {
	n = sscanf(inner_tok, "%d|%d", &hap1, &hap2);
  phase = 1;
	if(n != 2) {
	  /* try with '/' separator instead */
	  n = sscanf(inner_tok, "%d/%d", &hap1, &hap2);

	  if(n == 2) {
      phase = 0;
	    if(warn_phase) {
	      my_warn("%s:%d: some genotypes are unphased (delimited "
		      "with '/' instead of '|')\n", __FILE__, __LINE__,
		      inner_tok);
	      warn_phase = FALSE;
	    }
	  } else {
	    if(warn_parse) {
	      my_warn("%s:%d: could not parse some genotype "
		      "strings that look like: '%s'\n", __FILE__, __LINE__,
		      inner_tok);
	      warn_parse = FALSE;
	    }
	    hap1 = VCF_GTYPE_MISSING;
	    hap2 = VCF_GTYPE_MISSING;
	  }
	}

	if((hap1 != VCF_GTYPE_MISSING && hap1 != 0 && hap1 != 1)  ||
	   (hap2 != VCF_GTYPE_MISSING && hap2 != 0 && hap2 != 1)) {

	  /* Copy number polymorphisms and multi-allelic SNPs
	   * can have values other than 0 and 1 (e.g. 3, 4, ...).
	   * Combined haplotype test does not currently deal with 
	   * these. Set the genotypes to MISSING (-1)
	   */
	  hap1 = VCF_GTYPE_MISSING;
	  hap2 = VCF_GTYPE_MISSING;
	}

	if((n_haps + 2) > expect_haps) {
	  my_err("%s:%d: more genotypes per line than expected (line: %ld)",
		 __FILE__, __LINE__, vcf_info->cur_line);
	}
	haplotypes[n_haps] = hap1;
	haplotypes[n_haps+1] = hap2;
  haplotypes_phase[n_phase] = phase;

	n_haps += 2;
  n_phase += 1;
      }
      
      i++;
    }
  }

  if(n_haps != expect_haps) {
    my_err("%s:%d: expected %ld genotype values per line, but got "
	   "%ld (line: %ld)", __FILE__, __LINE__,
	   expect_haps, n_haps, vcf_info->cur_line);
  }
  if(n_phase != expect_phase) {
    my_err("%s:%d: expected %ld phase values per line, but got "
     "%ld (line: %ld)", __FILE__, __LINE__,
     expect_phase, n_phase, vcf_info->cur_line);
  }
}



/**
 * get genotype probabilities by parsing and converting genotype likelihoods
 * (GL) from VCF line
 */
void vcf_parse_gl(VCFInfo *vcf_info, float *geno_probs, char *cur, long gl_idx) {
  char delim[] = "\t";
  char inner_delim[] = ":";
  char *tok, *inner_tok, *inner_cur;
  char gtype[VCF_MAX_FORMAT];
  long  i, n, n_geno_probs, expect_geno_probs;
  float like_homo_ref, like_het, like_homo_alt;
  float prob_homo_ref, prob_het, prob_homo_alt, prob_sum;

  expect_geno_probs = vcf_info->n_sample * 3;
  
  n_geno_probs = 0;
  
  while((tok = strsep(&cur, delim)) != NULL) {
    /* each genotype string is delimited by ':'
     * each GL portion is delimited by ','
     */
    util_strncpy(gtype, tok, sizeof(gtype));

    i = 0;
    inner_cur = gtype;
    while((i <= gl_idx) && (inner_tok = strsep(&inner_cur, inner_delim)) != NULL) {
      if(i == gl_idx) {
	n = sscanf(inner_tok, "%g,%g,%g", &like_homo_ref, &like_het,
		   &like_homo_alt);

	if(n != 3) {
	  if(strcmp(inner_tok, ".") == 0) {
	    /* '.' indicates missing data
	     * set all likelihoods to log(0.333) = -0.477
	     */
	    like_homo_ref = like_het = like_homo_alt = -0.477;
	  } else {
	    my_err("%s:%d: failed to parse genotype likelihoods from "
		   "string '%s' (line: %ld)", __FILE__, __LINE__,
		   inner_tok, vcf_info->cur_line);
	  }
	}

	/* convert log10(prob) to prob */
	prob_homo_ref = pow(10.0, like_homo_ref);
	prob_het = pow(10.0, like_het);
	prob_homo_alt = pow(10.0, like_homo_alt);

	if((n_geno_probs + 3) > expect_geno_probs) {
	  my_err("%s:%d: more genotype likelihoods per line "
		 "than expected (line: %ld)",
		 __FILE__, __LINE__, vcf_info->cur_line);
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

      i++;
    }
  }

  if(n_geno_probs != expect_geno_probs) {
    my_err("%s:%d: expected %ld genotype likelihoods per line, but got "
	   "%ld (line: %ld)", __FILE__, __LINE__, expect_geno_probs,
	   n_geno_probs, vcf_info->cur_line);
  }
}  


/**
 * get genotype probabilities by parsing GP token from VCF line
 */
void vcf_parse_gp(VCFInfo *vcf_info, float *geno_probs, char *cur, long gp_idx) {
  char delim[] = "\t";
  char inner_delim[] = ":";
  char *tok, *inner_tok, *inner_cur;
  char gtype[VCF_MAX_FORMAT];
  long  i, n, n_geno_probs, expect_geno_probs;
  float prob_homo_ref, prob_het, prob_homo_alt, prob_sum;

  expect_geno_probs = vcf_info->n_sample * 3;
  
  n_geno_probs = 0;
  
  while((tok = strsep(&cur, delim)) != NULL) {
    /* each genotype string is delimited by ':'
     * each GP portion is delimited by ','
     */
    util_strncpy(gtype, tok, sizeof(gtype));

    i = 0;
    inner_cur = gtype;
    while((i <= gp_idx) && (inner_tok = strsep(&inner_cur, inner_delim)) != NULL) {
      if(i == gp_idx) {
	n = sscanf(inner_tok, "%g,%g,%g", &prob_homo_ref, &prob_het,
		   &prob_homo_alt);

	if(n != 3) {
	  if(strcmp(inner_tok, ".") == 0) {
	    /* '.' indicates missing data
	     * set all probabilities to 0.333
	     */
	    prob_homo_ref = prob_het = prob_homo_alt = 0.333;
	  } else {
	    my_err("%s:%d: failed to parse genotype probabilities from "
		   "string '%s' (line: %ld)",
		   __FILE__, __LINE__, inner_tok, vcf_info->cur_line);
	  }
	}
	
	/* check that probs sum to 1.0, normalize if they don't */
	prob_sum = prob_homo_ref + prob_het + prob_homo_alt;
	if((prob_sum > 1.001) || (prob_sum < 0.999)) {
	  prob_homo_ref = prob_homo_ref / prob_sum;
	  prob_het = prob_het / prob_sum;
	  prob_homo_alt = prob_homo_alt / prob_sum;
     	}
	geno_probs[n_geno_probs] = prob_homo_ref;
	geno_probs[n_geno_probs + 1] = prob_het;
	geno_probs[n_geno_probs + 2] = prob_homo_alt;

	n_geno_probs += 3;
      }

      i++;
    }
  }

  if(n_geno_probs != expect_geno_probs) {
    my_err("%s:%d: expected %ld genotype probabilities per line, but got "
	   "%ld (line: %ld)", __FILE__, __LINE__, expect_geno_probs,
	   n_geno_probs, vcf_info->cur_line);
  }
}  


void vcf_parse_geno_probs(VCFInfo *vcf_info, float *geno_probs, char *cur) {
  long gl_idx, gp_idx;

  /* get index of GP and GL tokens in format string */
  gp_idx = get_format_index(vcf_info->format, "GP");
  gl_idx = get_format_index(vcf_info->format, "GL");

  if(gl_idx == -1) {
    /* PL is same as GL, but rounded to nearest integer */
    gl_idx = get_format_index(vcf_info->format, "PL");
  }

  if((gl_idx == -1) && (gp_idx == -1)) {
    my_err("%s:%d: VCF format string does not specify GL, GP, or PL token "
	   "so cannot obtain genotype probabilities. Format string: '%s'.\n"
	   "To use this file, you must run snp2h5 without "
	   "the --geno_prob option.", __FILE__, __LINE__,
	   vcf_info->format);
  }

  if(gp_idx > -1) {
    vcf_parse_gp(vcf_info, geno_probs, cur, gp_idx);
    return;
  }

  vcf_parse_gl(vcf_info, geno_probs, cur, gl_idx);  
}



/**
 * Gets next line of VCF file and parses it into VCFInfo datastructure.
 *
 * If geno_probs array is non-null genotype likelihoods are parsed and
 * stored in the provided array. The array must be of length
 * n_sample*3.
 *
 * If haplotypes array is non-null phased genotypes are parsed and
 * stored in the provided array. The array must be of length
 * n_sample*2.
 *
 * Returns 0 on success, -1 if at EOF.
 */
int vcf_read_line(gzFile vcf_fh, VCFInfo *vcf_info, SNP *snp,
		  float *geno_probs, char *haplotypes, char *haplotypes_phase) {
  char *cur, *token;
  int n_fix_header, ref_len, alt_len;
  size_t tok_num;

  /* Used to allow space or tab delimiters here but now only allow
   * tab.  This is because VCF specification indicates that fields
   * should be tab-delimited, and occasionally some fields contain
   * spaces.
   */
  /* const char delim[] = " \t";*/
  const char delim[] = "\t";

  n_fix_header = sizeof(vcf_fix_headers) / sizeof(const char *);

  /* read a line */
  if(util_gzgetline(vcf_fh, &vcf_info->buf, &vcf_info->buf_size) == -1) {
    return -1;
  }
  vcf_info->cur_line += 1;
    
  cur = vcf_info->buf;
  tok_num = 0;

  /* chrom */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }

  util_strncpy(snp->chrom, token, sizeof(snp->chrom));
  
  
  /* pos */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }
  snp->pos = util_parse_long(token);
  
  /* ID */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }
  util_strncpy(snp->name, token, sizeof(snp->name));
  
  /* ref */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n", n_fix_header,
	   vcf_info->cur_line);
  }
  ref_len = util_strncpy(snp->allele1, token, sizeof(snp->allele1));

  /* used to warn about truncations, but makes program too
   * chatty if there are a lot of them
   */
  vcf_info->ref_len = 0;
  /* vcf_info->ref_len = strlen(token); */
  /* if(ref_len != vcf_info->ref_len) { */
  /*   my_warn("truncating long allele (%ld bp) to %ld bp\n", */
  /* 	    vcf_info->ref_len, ref_len); */
  /* } */
  
  /* alt */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }
  alt_len = util_strncpy(snp->allele2, token, sizeof(snp->allele2));
  
  vcf_info->alt_len = 0;
  /* vcf_info->alt_len = strlen(token); */
  /* if(alt_len != vcf_info->alt_len) { */
  /*   my_warn("truncating long allele (%ld bp) to %ld bp\n", */
  /* 	    vcf_info->alt_len, alt_len); */
  /* } */

  /* qual */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }
  util_strncpy(vcf_info->qual, token, sizeof(vcf_info->qual));

  /* filter */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n", n_fix_header,
	   vcf_info->cur_line);
  }
  util_strncpy(vcf_info->filter, token, sizeof(vcf_info->filter));


  /* info */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line (line: %ld)\n",
	   n_fix_header, vcf_info->cur_line);
  }
  util_strncpy(vcf_info->info, token, sizeof(vcf_info->info));

  if(vcf_info->has_format) {
    /* format */
    token = strsep(&cur, delim);
    if(token == NULL) {
      my_err("expected at least %d tokens per line (line: %ld)\n",
	     n_fix_header, vcf_info->cur_line);
    }
    util_strncpy(vcf_info->format, token, sizeof(vcf_info->format));

    /* now parse haplotypes and/or genotype likelihoods */
    if(geno_probs && haplotypes) {
      char *cur_copy;    
      /* Both genotype probs and haplotypes requested.
       * Need to copy string because it is modified
       * by the tokenizing in the parsing functions.
       *
       * This could be made more efficient by doing the parsing
       * of both types of data at same time
       */
      cur_copy = my_malloc(strlen(cur)+1);
      strcpy(cur_copy, cur);
    
      vcf_parse_geno_probs(vcf_info, geno_probs, cur_copy);
      my_free(cur_copy);

      vcf_parse_haplotypes(vcf_info, haplotypes, cur, haplotypes_phase);
    } else if(geno_probs) {
      vcf_parse_geno_probs(vcf_info, geno_probs, cur);
    } else if(haplotypes) {
      vcf_parse_haplotypes(vcf_info, haplotypes, cur, haplotypes_phase);
    }
  } else {
    /* header specifies no FORMAT string */
    if(geno_probs || haplotypes) {
      my_err("cannot parse genotypes or haplotypes without FORMAT field");
    }
  }

  /* my_free(line); */

  return 0;
}
