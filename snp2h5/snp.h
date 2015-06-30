#ifndef __SNP_H__
#define __SNP_H__


#define SNP_MAX_ALLELE 32
#define SNP_MAX_CHROM 32
#define SNP_MAX_NAME 16


typedef struct {
  char name[SNP_MAX_NAME];
  char chrom[SNP_MAX_CHROM];
  long pos;
  char allele1[SNP_MAX_ALLELE];
  char allele2[SNP_MAX_ALLELE];
} SNP;



#endif
