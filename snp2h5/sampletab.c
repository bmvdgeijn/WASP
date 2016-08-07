#include <string.h>
#include <stdio.h>

#include <hdf5.h>
#include <hdf5_hl.h>


#include "memutil.h"
#include "sampletab.h"
#include "err.h"
#include "util.h"


/**
 * Allocates memory for and initializes data structure
 * that describes the HDF5 table used to store Samples data.
 *
 */
SampleTab *sample_tab_new(hid_t h5file, size_t n_record) {
  herr_t status;
  SampleTab *tab;
  Sample sample_desc;

  const char *field_names[] = {"name"};

  tab = my_malloc(sizeof(SampleTab));
  
  tab->h5file = h5file;  
  
  /* set datatypes for each field */
  tab->name_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(tab->name_type, SAMPLE_MAX_NAME);
  tab->field_type[0] = tab->name_type; /* name */

  /* sizes of record and each field */
  tab->record_size = sizeof(Sample);
  tab->field_size[0] = sizeof(sample_desc.name);
  
  /* offsets of each field */
  tab->field_offset[0] = HOFFSET(Sample, name);
    
  /* title and name of table */
  tab->title = util_str_dup("samples");
  tab->name = util_str_dup("samples");

  /* set chunk size and compression */
  tab->chunk_size = SAMPLETAB_CHUNK_SIZE;
  tab->compress = 1;

  tab->n_record = 0;
  
  status = H5TBmake_table(tab->title, tab->h5file, tab->name,
			  SAMPLETAB_N_FIELDS, n_record,
			  tab->record_size, field_names,
			  tab->field_offset, tab->field_type,
			  tab->chunk_size, NULL, tab->compress, NULL);
  
  if(status < 0) {
    my_err("%s:%d: could not create samples table\n", __FILE__, __LINE__);
  }

  return tab;
}


void sample_tab_free(SampleTab *tab) {
  int i;
  
  H5Tclose(tab->name_type);
  my_free(tab->title);
  my_free(tab->name);
  
  my_free(tab);
}


/**
 * Creates a new samples table HDF5 file pointed to by h5file handle
 * and populates it using provided array of samples
 */
SampleTab *sample_tab_create(hid_t h5file, Sample *samples, size_t n_sample) {
  SampleTab *tab;
  int i;

  tab = sample_tab_new(h5file, n_sample);
  
  for(i = 0; i < n_sample; i++) {
    sample_tab_append_row(tab, &samples[i]);
  }
  
  return tab;
}



/**
 * Appends row to table described by provided SampleTab datastructure.
 */
void sample_tab_append_row(SampleTab *tab, Sample *data) {
  herr_t status;
  hsize_t n_to_write;

  n_to_write = 1;

  /* fprintf(stderr, "writing sample record:\n" */
  /* 	  "data: %s\n"  */
  /* 	  "n_record: %d\n" */
  /* 	  "record_size: %d\n" */
  /* 	  "field_offset[0]: %d\n" */
  /* 	  "field_size[0]: %d\n", data->name, tab->n_record, tab->record_size, */
  /* 	  tab->field_offset[0], tab->field_size[0]); */
  
  status = H5TBwrite_records(tab->h5file, tab->name,
			     tab->n_record, n_to_write,
			     tab->record_size, tab->field_offset,
			     tab->field_size, data);

  tab->n_record += 1;
  
  if(status < 0) {
    my_err("%s:%d: failed to write record to Sample table\n",
	   __FILE__, __LINE__);
  }
}

