#ifndef __SNPTAB_H__
#define __SNPTAB_H__

#include <hdf5.h>
#include "vcf.h"


#define SNPTAB_NFIELDS 5
#define SNPTAB_MAX_ALLELE 32
#define SNPTAB_MAX_NAME 16
#define SNPTAB_N_FIELDS 4
#define SNPTAB_CHUNK_SIZE 1000


/* SNPDesc holds data for a single record in table */
typedef struct {
  char name[SNPTAB_MAX_NAME];
  long pos;
  char allele1[SNPTAB_MAX_ALLELE];
  char allele2[SNPTAB_MAX_ALLELE];
} SNPDesc;


/* SNPTab holds information about SNP table 
 * and datatypes of each field.
 */

typedef struct {
  hid_t h5file;
  size_t record_size;
  char *name;
  char *title;
  hid_t field_type[SNPTAB_N_FIELDS];
  size_t field_size[SNPTAB_N_FIELDS];
  size_t field_offset[SNPTAB_N_FIELDS];
  
  hid_t name_type;
  hid_t allele_type;
  
  int compress;
  size_t chunk_size;
} SNPTab;



SNPTab *snp_tab_new(hid_t h5file, const char *chrom_name);
void snp_tab_free(SNPTab *tab);
void snp_tab_append_row(SNPTab *tab, SNPDesc *data);

void snp_tab_append_vcf_row(SNPTab *tab, VCFInfo *vcf_info);

#endif
