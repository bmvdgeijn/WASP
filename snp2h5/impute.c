

#include <zlib.h>
#include <string.h>

#include "memutil.h"
#include "util.h"
#include "impute.h"
#include "snp.h"
#include "err.h"



ImputeInfo *impute_info_new() {
  ImputeInfo *impute_info;

  impute_info = my_malloc(sizeof(ImputeInfo));
  impute_info->buf_size = 1024;
  impute_info->buf = my_malloc(impute_info->buf_size);

  return impute_info;
}

void impute_info_free(ImputeInfo *impute_info) {
  my_free(impute_info->buf);
  my_free(impute_info);
}


/* 
 * Peaks at next line in file and counts number of fields.
 * Can be used to determine number of samples in file.
 * File is rewound to beginning after.
 * Returns -1 on failure.
*/
long impute_count_fields(gzFile fh) {
  char *line, *cur, *tok;
  char delim[] = " \t";
  long n;
  
  /* read a line */
  line = util_gzgets_line(fh);

  if(line == NULL) {
    my_err("%s:%d: no lines left in file, cannot count fields\n",
	   __FILE__, __LINE__);
  }

  cur = line;
  n = 0;
  while((tok = strsep(&cur, delim)) != NULL) {
    n++;
  }

  my_free(line);

  /* rewind to beginning of file */
  if(gzseek(fh, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }
  
  return n;
}


void impute_parse_haplotypes(char *haplotypes, char *cur, long n_sample) {
  long expect_n, i;
  char *tok;
  char delim[] = " \t";
  
  expect_n = n_sample * 2;

  i = 0;
  while((tok = strsep(&cur, delim)) != NULL) {
    if(i >= expect_n) {
      my_err("%s:%d: more haplotype values per row than expected (%ld)."
	     " Current token: '%s'", __FILE__, __LINE__, expect_n, tok);
    }

    /* expect values to be either 0 or 1 */
    if(strcmp(tok, "0") == 0) {
      haplotypes[i] = 0;
    }
    else if(strcmp(tok, "1") == 0) {
      haplotypes[i] = 1;
    }
    else {
      my_err("%s:%d: expected haplotype values to be either 0 or 1 but got "
	     "%s", __FILE__, __LINE__, tok);
    }
    
    i++;
  }

  if(i != expect_n) {
    my_err("%s:%d: fewer haplotype values (%ld) per row than expected (%ld)",
	   __FILE__, __LINE__, expect_n);
  }
}



void impute_parse_geno_probs(float *geno_probs, char *cur, long n_sample) {
  long expect_n, i;
  char delim[] = " \t";
  char *tok;
  
  expect_n = n_sample * 3;

  i = 0;
  while((tok = strsep(&cur, delim)) != NULL) {
    if(i >= expect_n) {
      my_err("%s:%d: more genotype_probs per row than expected (%ld)",
	     __FILE__, __LINE__, expect_n);
    }
    geno_probs[i] = util_parse_double(tok);
    i++;
  }

  if(i != expect_n) {
    my_err("%s:%d: fewer genotype_probs (%ld) per row than expected (%ld)",
	   __FILE__, __LINE__, expect_n);
  }
}



/**
 * Gets next line of IMPUTE file and parses it into ImputeInfo datastructure.
 *
 * IMPUTE files are described here:
 * http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
 * 
 * example line: 
 *     --- rs149201999 16050408 T C 0.966 0.034 0 0.395 0.467 ....
 * 
 * If geno_probs array is non-null genotype probabilities are parsed and
 * stored in the provided array. The array must be of length
 * n_sample*3.
 *
 * If haplotypes array is non-null phased genotypes are parsed and
 * stored in the provided array. The array must be of length
 * n_sample*2.
 *
 * IMPUTE files contain EITHER haplotypes OR genotypes so only
 * one of geno_probs or haplotypes should be non-null (at most).
 *
 * Returns 0 on success, -1 if at EOF.
 */
int impute_read_line(gzFile fh, ImputeInfo *impute_info, SNP *snp,
		     float *geno_probs, char *haplotypes) {
  char *cur, *token;
  int alt_len;
  size_t tok_num;
  const char delim[] = " \t";

  /* read a line */
  if(util_gzgetline(fh, &impute_info->buf, &impute_info->buf_size) == -1) {
    return -1;
  }
  
  cur = impute_info->buf;
  tok_num = 0;

  /* SNP name, often just set to "---" */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", IMPUTE_FIX_HEADER);
  }

  /* SNP identifier (rs_id) */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", IMPUTE_FIX_HEADER);
  }
  util_strncpy(snp->name, token, sizeof(snp->name));
    
  /* pos */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", IMPUTE_FIX_HEADER);
  }
  snp->pos = util_parse_long(token);
  
  /* ref allele */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", IMPUTE_FIX_HEADER);
  }
  util_strncpy(snp->allele1, token, sizeof(snp->allele1));
    
  /* alt allele */
  token = strsep(&cur, delim);
  if(token == NULL) {
    my_err("expected at least %d tokens per line\n", IMPUTE_FIX_HEADER);
  }
  alt_len = util_strncpy(snp->allele2, token, sizeof(snp->allele2));
  
  /* now parse haplotypes and/or genotype likelihoods */
  if(geno_probs && haplotypes) {
    my_err("impute2 files contain EITHER genotypes or haplotypes, but "
	   "both requested\n");
  }
  else if(geno_probs) {
    impute_parse_geno_probs(geno_probs, cur, impute_info->n_sample);
  }
  else if(haplotypes) {
    impute_parse_haplotypes(haplotypes, cur, impute_info->n_sample);
  }

  return 0;
}
