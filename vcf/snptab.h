#ifndef __SNPTAB_H__
#define __SNPTAB_H__

#include <hdf5.h>

#include "snp.h"


#define SNPTAB_NFIELDS 5
#define SNPTAB_N_FIELDS 4
#define SNPTAB_CHUNK_SIZE 1000



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

void snp_tab_append_row(SNPTab *tab, SNP *data);

#endif
