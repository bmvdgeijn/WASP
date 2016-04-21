#ifndef __SNPTAB_H__
#define __SNPTAB_H__

#include <hdf5.h>

#include "snp.h"


#define SNPTAB_N_FIELDS 4

/* chunk size affects performance a lot
 * small chunks = much faster writing of tables, but
 * worse compression
 */

#define SNPTAB_CHUNK_SIZE 12



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
  hid_t chrom_type;
  hid_t allele_type;
  
  int compress;
  size_t chunk_size;


  /* all records in table are created at once (might be faster?) 
   * then overwritten. max_record is the number in the HDF5 
   * table and n_record is the number that have actually been written.
   */
  size_t n_record;
  size_t max_record;
  
  
} SNPTab;



SNPTab *snp_tab_new(hid_t h5file, const char *chrom_name,
		    size_t max_record);

void snp_tab_free(SNPTab *tab);

void snp_tab_append_row(SNPTab *tab, SNP *data);

#endif
