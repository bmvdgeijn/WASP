#ifndef __CHR_H__
#define __CHR_H__

typedef struct {
  int id;
  char *name;
  long len;
  char *assembly;
} Chromosome;


Chromosome *chrom_guess_from_file(const char *filename,
				  Chromosome *chroms,
				  int n_chrom);

void chrom_array_free(Chromosome *chrom, int n_chrom);
Chromosome *chrom_copy(const Chromosome *chrom);
void chrom_free(Chromosome *chrom);
Chromosome *chrom_read_file(const char *filename, int *n_chrom);

#endif
