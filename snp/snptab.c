#include <string.h>
#include <stdio.h>

#include <hdf5.h>
#include <hdf5_hl.h>


#include "memutil.h"
#include "snptab.h"
#include "vcf.h"
#include "err.h"
#include "util.h"



/**
 * Allocates memory for and initializes data structure
 * that describes the HDF5 table used to store SNP data.
 * A single SNPTab data structure should be created for 
 * each chromosome, and be freed after use.
 *
 */
SNPTab *snp_tab_new(hid_t h5file, const char *chrom_name) {
  herr_t status;
  SNPTab *tab;
  SNPDesc snp_desc;

  const char *field_names[] =
    {"name", "pos", "allele1", "allele2"};

  tab = my_malloc(sizeof(SNPTab));
  
  tab->h5file = h5file;  
  
  /* set datatypes for each field */
  tab->name_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(tab->name_type, SNPTAB_MAX_NAME);
  tab->allele_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(tab->allele_type, SNPTAB_MAX_ALLELE);    
  tab->field_type[0] = tab->name_type; /* name */
  tab->field_type[1] = H5T_NATIVE_LONG; /* pos */
  tab->field_type[2] = tab->allele_type; /* allele1 */
  tab->field_type[3] = tab->allele_type; /* allele2 */

  /* sizes of record and each field */
  tab->record_size = sizeof(SNPTab);
  tab->field_size[0] = sizeof(snp_desc.name);
  tab->field_size[1] = sizeof(snp_desc.pos);
  tab->field_size[2] = sizeof(snp_desc.allele1);
  tab->field_size[3] = sizeof(snp_desc.allele2);

  /* offsets of each field */
  tab->field_offset[0] = HOFFSET(SNPDesc, name);
  tab->field_offset[1] = HOFFSET(SNPDesc, pos);
  tab->field_offset[2] = HOFFSET(SNPDesc, allele1);
  tab->field_offset[3] = HOFFSET(SNPDesc, allele2);
    
  /* title and name of table */
  tab->title = util_str_concat(chrom_name, " SNPs", NULL);  
  tab->name = util_str_dup(chrom_name);

  /* set chunk size and compression */
  tab->chunk_size = SNPTAB_CHUNK_SIZE;
  tab->compress = 1;
  
  status = H5TBmake_table(tab->title, tab->h5file, tab->name,
			  SNPTAB_N_FIELDS, 0,
			  tab->record_size, field_names,
			  tab->field_offset, tab->field_type,
			  tab->chunk_size, NULL, tab->compress, NULL);
  
  if(status < 0) {
    my_err("%s:%d: could not create SNP table\n", __FILE__, __LINE__);
  }

  
  return tab;
}


void snp_tab_free(SNPTab *tab) {
  int i;
  
  H5Tclose(tab->allele_type);
  H5Tclose(tab->name_type);
  my_free(tab->title);
  my_free(tab->name);
  
  my_free(tab);
}


/**
 * Appends row to table described by provided SNPTab datastructure.
 */
void snp_tab_append_row(SNPTab *tab, SNPDesc *data) {
  herr_t status;
  hsize_t n_records;

  n_records = 1;

  status = H5TBappend_records(tab->h5file, tab->name,
			      n_records, tab->record_size,
			      tab->field_offset, tab->field_size,
			      data);

  if(status < 0) {
    my_err("%s:%d: failed to write record to SNP table\n",
	   __FILE__, __LINE__);
  }
}


/**
 * Copies data from VCF data structure and appends row to SNP table
 */
void snp_tab_append_vcf_row(SNPTab *tab, VCFInfo *vcf) {
  SNPDesc data;
  
  strncpy(data.name, vcf->id, SNPTAB_MAX_NAME);
  data.pos = vcf->pos;

  /* TODO: could warn when long alleles are truncated here... */
  strncpy(data.allele1, vcf->ref_allele, SNPTAB_MAX_ALLELE);
  strncpy(data.allele2, vcf->alt_allele, SNPTAB_MAX_ALLELE);
  
  snp_tab_append_row(tab, &data);
}
