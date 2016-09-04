#ifndef __IMPUTE_H__
#define __IMPUTE_H__

#include <zlib.h>

#include "snp.h"


#define IMPUTE_FIX_HEADER 5



typedef struct {
  int n_sample;

  /* used for reading lines */
  size_t buf_size;
  char *buf;
  
  /* could store lots of header info here */
} ImputeInfo;



ImputeInfo *impute_info_new();
void impute_info_free(ImputeInfo *impute_info);

int impute_read_line(gzFile fh, ImputeInfo *impute_info, SNP *snp,
		     float *geno_probs, char *haplotypes);

long impute_count_fields(gzFile fh);

#endif
