#ifndef __SAMPLETAB_H__
#define __SAMPLETAB_H__

#include <hdf5.h>

#include "sample.h"


#define SAMPLETAB_N_FIELDS 1

/* chunk size affects performance a lot
 * small chunks = much faster writing of tables, but
 * worse compression
 */

#define SAMPLETAB_CHUNK_SIZE 100


/* SampleTab holds information about Samples
 * and datatypes of each field. Currently
 * there is only a single field which is the 
 * name of the sample, however you could imagine
 * adding other information such as sex, 
 * age, population, etc.
 */

typedef struct {
  hid_t h5file;
  size_t record_size;
  char *name;
  char *title;
  hid_t field_type[SAMPLETAB_N_FIELDS];
  size_t field_size[SAMPLETAB_N_FIELDS];
  size_t field_offset[SAMPLETAB_N_FIELDS];
  
  hid_t name_type;
  
  int compress;
  size_t chunk_size;

  size_t n_record;
} SampleTab;


SampleTab *sample_tab_new(hid_t h5file, const char *chrom_name,
			  size_t n_record);

SampleTab *sample_tab_create(hid_t h5file, const char *chrom_name,
			     Sample *samples, size_t n_sample);

SampleTab *sample_tab_from_names(hid_t h5file, const char *chrom_name,
				 char **sample_names,
				 size_t n_sample);

void sample_tab_free(SampleTab *tab);

void sample_tab_append_row(SampleTab *tab, Sample *data);

#endif
