#ifndef __VCF_H__
#define __VCF_H__


#include <zlib.h>


#define VCF_MAX_ID 1024
#define VCF_MAX_CHROM 1024
#define VCF_MAX_ALLELE 1024
#define VCF_MAX_QUAL 1024
#define VCF_MAX_FILTER 1024
#define VCF_MAX_FORMAT 1024

typedef struct {
  int n_samples;
  long n_header_lines;

  char chrom[VCF_MAX_CHROM];
  char id[VCF_MAX_ID]; /* rsid */
  long pos;

  /* limit allele buffer to fixed size MAX_ALLELE
   * for efficiency
   */
  char ref_allele[VCF_MAX_ALLELE];
  char alt_allele[VCF_MAX_ALLELE];

  /* records true length of ref / alt alleles, which can be
   * truncated by limited buffer size
   */
  size_t ref_len;
  size_t alt_len;

  char qual[VCF_MAX_QUAL];
  char filter[VCF_MAX_FILTER];
  char info[VCF_MAX_FILTER];
  char format[VCF_MAX_FORMAT];

  /* could store lots of header info here */
} VCFInfo;


void vcf_read_header(gzFile vcf_fh, VCFInfo *vcf_info);

int vcf_read_line(gzFile vcf_fh, VCFInfo *vcf_info, double *geno_probs);


#endif
