
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <hdf5.h>


#include "vcf.h"
#include "util.h"
#include "memutil.h"
#include "chrom.h"


/* 
 * These values control size of chunks that are used for 
 * data transfer and compression. Changing them
 * affects effciency of reading/writing time and compression.
 * Not sure what optimal values are...
 */
#define ROW_CHUNK 10
#define COL_CHUNK 1000


typedef struct {
  /* chromInfo.txt.gz file that chromosomes are read from */
  char *chrom_file;

  /* HDF5 file that genotype probabilities are written to */
  char *geno_prob_file;
  char **vcf_files;

  int n_vcf_files;
} Arguments;



typedef struct {
  hid_t h5file;  /* HDF5 file handle */
  hid_t dataset_prop; /* dataset properties list */
  hid_t file_dataspace; /* dataspace describing file layout */
  hid_t mem_dataspace;   /* dataspace describing memory layout */
  hid_t dataset;        /* dataset */

  hsize_t n_col; /* number of columns */
  hsize_t n_row; /* number of rows */
  hsize_t file_dims[2]; /* file dataset dimensions */
  hsize_t mem_dims[2]; /* dimensions stored in mem */
  hsize_t max_dims[2]; /* maximum dimensions */
  hsize_t chunk[2]; /* chunk size */
  hsize_t count[2];
  hsize_t offset[2];
} H5MatrixInfo;


void usage(char **argv) {
  fprintf(stderr, "usage: %s --chrom <chromInfo.txt.gz> "
	  "--geno_prob <geno_probs_output_file.h5> "
	  "<chr1.vcf> [<chr2.vcf> ...]\n", argv[0]);
  exit(-1);
}


void parse_args(Arguments *args, int argc, char **argv) {
  int c, i;
  
   static struct option loptions[] = {
     {"chrom", required_argument, 0, 'c'},
     {"geno_prob", required_argument, 0, 'p'},
     {0,0,0,0}
   };

   args->chrom_file = NULL;
   args->geno_prob_file = NULL;

   
   while(1) {
     c = getopt_long(argc, argv, "c:p:", loptions, NULL);
     
     if(c == -1) {
       break;
     }
     
     switch (c) {
       case 'c':
	 args->chrom_file = util_str_dup(optarg);
	 break;
       case 'p':
	 args->geno_prob_file = util_str_dup(optarg);
	 break;
       default:
	 usage(argv);
	 break;
     }
    
   }
      
   if(optind < argc) {
     args->n_vcf_files = argc - optind;
     args->vcf_files = &argv[optind];
   } else {
     fprintf(stderr, "must specify at least one VCF file\n");
     usage(argv);
   }

   if(!args->chrom_file) {
     fprintf(stderr, "--chrom argument is missing\n");
     usage(argv);
   }

   if(!args->geno_prob_file) {
     fprintf(stderr, "--geno_prob argument is missing\n");
     usage(argv);
   }   
}



/**
 * Returns appropriate chromosome for filename by looking for match
 * with chromosome name in filename. Returns NULL if no match found.
 */
Chromosome *guess_chrom(char *filename, Chromosome *chroms, int n_chrom) {
  int i, j, longest_match, n1, n2;
  Chromosome *match_chrom;

  longest_match = 0;
  match_chrom = NULL;
  
  for(i = 0; i < n_chrom; i++) {
    n1 = strlen(chroms[i].name);

    if(n1 > longest_match) {
      /* is this longest-matching chromosome name?
       * Use longest as best match because otherwise "chr1" will match
       * filenames that contain "chr10", etc.
       */
      n2 = strlen(filename);
      for(j = 0; j < n2 - n1 + 1; j++) {
	if(strncmp(chroms[i].name, &filename[j], n1) == 0) {
	  /* chromosome name is present in this filename */
	  match_chrom = &chroms[i];
	  longest_match = n1;
	  break;
	}
      }
    }
  }

  return match_chrom;
}



hid_t create_h5file(const char *filename) {
  hid_t h5file;
  
  h5file = H5Fcreate(filename, H5F_ACC_TRUNC,
		     H5P_DEFAULT, H5P_DEFAULT);
  if(h5file < 0) {
    my_err("%s:%d: could not create file '%s'", __FILE__, __LINE__,
	   filename);
    return -1;
  }

  return h5file;
}



void init_h5matrix(H5MatrixInfo *info,
		   hsize_t n_row, hsize_t n_col,
		   const char *name) {
  int rank = 2;
  herr_t status;

  info->n_row = n_row;
  info->n_col = n_col;
  
  /* define a matrix dataspace which can hold genotype probs, etc... */
  info->file_dims[0] = n_row;
  info->file_dims[1] = n_col;
  info->max_dims[0] = n_row;
  info->max_dims[1] = n_col;
  info->file_dataspace = H5Screate_simple(rank, info->file_dims,
					  info->max_dims);
  if(info->file_dataspace < 0) {
    my_err("%s:%d: failed to create file dataspace", __FILE__,
	   __LINE__);
  }

  /* in memory we will work with 1 entire row at a time */
  info->mem_dims[0] = 1;
  info->mem_dims[1] = n_col;
  info->mem_dataspace = H5Screate_simple(rank, info->mem_dims, NULL);
  if(info->mem_dataspace < 0) {
    my_err("%s:%d: failed to create mem dataspace", __FILE__,
	   __LINE__);
  }
    
  /* create property list for dataset */
  info->dataset_prop = H5Pcreate(H5P_DATASET_CREATE);
  if(info->dataset_prop < 0) {
    my_err("%s:%d: failed to create dataset property list",
	   __FILE__, __LINE__);
  }
  
  /* set chunking properties of dataset */
  info->chunk[0] = ROW_CHUNK;
  info->chunk[1] = COL_CHUNK;
  status = H5Pset_chunk(info->dataset_prop, rank, info->chunk);
  if(status < 0) {
    my_err("%s:%d: failed to set chunksize", __FILE__, __LINE__);
  }
  /* use zlib compression level 6 */
  status = H5Pset_deflate(info->dataset_prop, 6);
  if(status < 0) {
    my_err("%s:%d: failed to set compression filter", __FILE__, __LINE__);
  }
  
  /* create new dataset */
  info->dataset = H5Dcreate(info->h5file, name, H5T_NATIVE_FLOAT,
			    info->file_dataspace,
			    info->dataset_prop);
  if(info->dataset < 0) {
    my_err("%s:%d failed to create dataset\n", __FILE__, __LINE__);
  }  
}



void write_h5matrix_row(H5MatrixInfo *info, hsize_t row_num, double *row_data) {
  hsize_t offset[2];
  hsize_t count[1];
  herr_t status;
  
  /* select a hyperslab that corresponds to a single row (all columns) */
  offset[0] = row_num;
  offset[1] = 0;
  count[0] = 1;
  count[1] = info->n_col;
  status = H5Sselect_hyperslab(info->file_dataspace, H5S_SELECT_SET, offset,
			       NULL, count, NULL);
      
  if(status < 0) {
    my_err("%s:%d: failed to select hyperslab for row %ld\n",
	   __FILE__, __LINE__, row_num);
  }

  /* write a row of data */
  status = H5Dwrite(info->dataset, H5T_NATIVE_FLOAT, info->mem_dataspace,
		    info->file_dataspace, H5P_DEFAULT, row_data);
  
  if(status < 0) {
    my_err("%s:%d: failed to write data for row %ld", __FILE__, __LINE__,
	   row_num);
  }
}




void close_h5matrix(H5MatrixInfo *info) {
    H5Sclose(info->file_dataspace);
    H5Sclose(info->mem_dataspace);
    H5Dclose(info->dataset);
}



int main(int argc, char **argv) {
  VCFInfo vcf;
  Arguments args;
  Chromosome *all_chroms, *chrom;
  int n_chrom;
  H5MatrixInfo gprob_info;
  hsize_t n_col, n_row, row;
  gzFile gzf;
  double *geno_probs;
  size_t n_lines;
  int i;

  parse_args(&args, argc, argv);

  /* read chromosomes */
  all_chroms = chrom_read_file(args.chrom_file, &n_chrom);
  
  /* create new HDF5 file */
  gprob_info.h5file = create_h5file(args.geno_prob_file);
    
  for(i = 0; i < args.n_vcf_files; i++) {
    chrom = guess_chrom(args.vcf_files[i], all_chroms, n_chrom);
    if(chrom == NULL) {
      my_err("%s:%d: could not guess chromosome from filename '%s'\n",
	     __FILE__, __LINE__, args.vcf_files[i]);
    }
    fprintf(stderr, "chromosome: %s, length: %ldbp\n", chrom->name, chrom->len);
    
    fprintf(stderr, "reading from VCF file %s\n", args.vcf_files[i]);
    gzf = util_must_gzopen(args.vcf_files[i], "rb");
    
    fprintf(stderr, "counting lines in file\n");
    n_lines = util_gzcount_lines(gzf);
    fprintf(stderr, "there are %d lines\n", n_lines);
    
    vcf_read_header(gzf, &vcf);
    fprintf(stderr, "num header lines: %ld\n", vcf.n_header_lines);
    
    n_row = n_lines - vcf.n_header_lines;
    
    n_col = vcf.n_samples * 3;
    geno_probs = my_malloc(n_col * sizeof(double));

    init_h5matrix(&gprob_info, n_row, n_col, chrom->name);
    
    row = 0;
    while(vcf_read_line(gzf, &vcf, geno_probs) != -1) {
      
       /* fprintf(stderr, "chrom: %s, id: %s, pos: %ld, ref: %s, alt: %s " */
       /* 	      "qual: %s, filter: %s, info: %s format: %s\n", */
       /* 	      vcf.chrom, vcf.id, vcf.pos, vcf.ref_allele, vcf.alt_allele, */
       /* 	      vcf.qual, vcf.filter, vcf.info, vcf.format); */
      write_h5matrix_row(&gprob_info, row, geno_probs);
      
      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }
    fprintf(stderr, "\n");

    if(row != n_row) {
      my_warn("%s:%d: expected %ld data rows, but only read %ld\n", __FILE__,
	      __LINE__, n_row, row);
    }
    
    my_free(geno_probs);
    close_h5matrix(&gprob_info);
    gzclose(gzf);
  }
  
  chrom_array_free(all_chroms, n_chrom);
  H5Fclose(gprob_info.h5file);

  
  return 0;
}
