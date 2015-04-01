
#include <zlib.h>
#include <stdlib.h>

#include <hdf5.h>

#include "vcf.h"
#include "util.h"
#include "memutil.h"


/* 
 * These values control size of chunks that are used for 
 * data transfer and compression. Changing them
 * affects effciency of reading/writing time and compression.
 * Not sure what optimal values are...
 */
#define ROW_CHUNK 10
#define COL_CHUNK 1000


int main(int argc, char **argv) {
  VCFInfo vcf;
  hid_t h5file, dprop, file_dataspace, mem_dataspace, dataset;
  hsize_t stride[2] = {1,1};
  hsize_t n_col, n_row, row;
  hsize_t file_dims[2], mem_dims[2], max_dims[2], chunk[2], size[2], count[2], offset[2];
  herr_t status;
  int rank = 2;
  gzFile gzf;
  double *geno_probs;
  size_t n_lines;
  int i;

  if(argc < 2) {
    fprintf(stderr, "usage: %s <vcf_file1> [<vcf_file2> ...]\n", argv[0]);
    exit(2);
  }

  /* create new HDF5 file */
  h5file = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(h5file < 0) {
    my_err("%s:%d: could not create file 'test.h5'", __FILE__, __LINE__);
    return -1;
  }
  
  for(i = 1; i < argc; i++) {    
    gzf = util_must_gzopen(argv[i], "rb");

    fprintf(stderr, "counting lines in file\n");
    n_lines = util_gzcount_lines(gzf);
    fprintf(stderr, "there are %d lines\n", n_lines);
    
    vcf_read_header(gzf, &vcf);
    fprintf(stderr, "num header lines: %ld\n", vcf.n_header_lines);
    
    n_row = n_lines - vcf.n_header_lines;
    
    n_col = vcf.n_samples * 3;
    geno_probs = my_malloc(n_col * sizeof(double));

    /* define a matrix dataspace  which can hold genotype probs */
    file_dims[0] = n_row;
    file_dims[1] = n_col;
    max_dims[0] = n_row;
    max_dims[1] = n_col;
    file_dataspace = H5Screate_simple(rank, file_dims, max_dims);
    if(file_dataspace < 0) {
      my_err("%s:%d: failed to create file dataspace", __FILE__,
	     __LINE__);
    }

    mem_dims[0] = 1;
    mem_dims[1] = n_col;
    mem_dataspace = H5Screate_simple(rank, mem_dims, NULL);
    if(mem_dataspace < 0) {
      my_err("%s:%d: failed to create mem dataspace", __FILE__,
	     __LINE__);
    }
    
    /* create property list for dataset */
    dprop = H5Pcreate(H5P_DATASET_CREATE);
    if(dprop < 0) {
      my_err("%s:%d: failed to create dataset property list",
	     __FILE__, __LINE__);
    }
    /* set chunking properties of dataset */
    chunk[0] = ROW_CHUNK;
    chunk[1] = COL_CHUNK;
    status = H5Pset_chunk(dprop, rank, chunk);
    if(status < 0) {
      my_err("%s:%d: failed to set chunksize", __FILE__, __LINE__);
    }
    /* use zlib compression level 6 */
    status = H5Pset_deflate(dprop, 6);
    if(status < 0) {
      my_err("%s:%d: failed to set compression filter", __FILE__, __LINE__);
    }
    
    /*** TODO: need to use gz compression filter here? ****/
    
    /* create new dataset */
    dataset = H5Dcreate(h5file, "test", H5T_NATIVE_FLOAT, file_dataspace, dprop);
        
    row = 0;
    while(vcf_read_line(gzf, &vcf, geno_probs) != -1) {
      
       /* fprintf(stderr, "chrom: %s, id: %s, pos: %ld, ref: %s, alt: %s " */
       /* 	      "qual: %s, filter: %s, info: %s format: %s\n", */
       /* 	      vcf.chrom, vcf.id, vcf.pos, vcf.ref_allele, vcf.alt_allele, */
       /* 	      vcf.qual, vcf.filter, vcf.info, vcf.format); */

      
      /* select a hyperslab that corresponds to length of 
       * row we want to write to
       */
      offset[0] = row;
      offset[1] = 0;
      count[0] = 1;
      count[1] = n_col;
      status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset,
				   NULL, count, NULL);
      
      if(status < 0) {
	my_err("%s:%d: failed to select hyperslab for row %ld\n",
	       __FILE__, __LINE__, row);
      }

      /* write the genotype probs */
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_dataspace,
			file_dataspace, H5P_DEFAULT, geno_probs);
      
      if(status < 0) {
	my_err("%s:%d: failed to write data", __FILE__, __LINE__, row);
      }
      
      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }

    fprintf(stderr, "\n");

    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Dclose(dataset);
	
    my_free(geno_probs);
    gzclose(gzf);
  }

  
  H5Fclose(h5file);

  
  return 0;
}
