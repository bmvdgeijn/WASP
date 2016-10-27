
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <hdf5.h>
#include <stdint.h>

#include "impute.h"
#include "vcf.h"
#include "sample.h"
#include "sampletab.h"
#include "snptab.h"
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
#define VEC_CHUNK 100000

#define GENO_PROB_DATATYPE H5T_NATIVE_FLOAT
#define HAPLOTYPE_DATATYPE H5T_NATIVE_CHAR
#define SNP_INDEX_DATATYPE H5T_NATIVE_LONG

#define SNP_INDEX_NONE -1
#define GENO_PROB_DEFAULT_VAL -1.0
#define HAPLOTYPE_DEFAULT_VAL -1

#define FORMAT_VCF 1
#define FORMAT_IMPUTE 2

typedef struct {
  /* chromInfo.txt.gz file that chromosomes are read from */
  char *chrom_file;

  /* flag indicating format of input files (FORMAT_VCF or FORMAT_IMPUTE) */
  int format;

  /* HDF5 files that SNP info written to */
  char *geno_prob_file; /* genotype probabilities */ 
  char *haplotype_file; /* haplotypes & phase */
  char *snp_index_file; /* base position => SNP table row lookup */
  char *snp_tab_file; /* SNP table with id, alleles, etc */

  char *sample_file;
  char **input_files;
  
  int n_input_files;
} Arguments;



typedef struct {
  hid_t h5file;  /* HDF5 file handle */
  hid_t dataset_prop; /* dataset properties list */
  hid_t file_dataspace; /* dataspace describing file layout */
  hid_t mem_dataspace;   /* dataspace describing memory layout */
  hid_t dataset;        /* dataset */

  hid_t datatype; /* float, int8, etc */

  hsize_t n_col; /* number of columns */
  hsize_t n_row; /* number of rows */
  hsize_t file_dims[2]; /* file dataset dimensions */
  hsize_t mem_dims[2]; /* dimensions stored in mem */
  hsize_t max_dims[2]; /* maximum dimensions */
  hsize_t chunk[2]; /* chunk size */
  hsize_t count[2];
  hsize_t offset[2];
} H5MatrixInfo;


typedef struct {
  hid_t h5file;  /* HDF5 file handle */
  hid_t dataset_prop; /* dataset properties list */
  hid_t file_dataspace; /* dataspace describing file layout */
  hid_t dataset;        /* dataset */
  hid_t datatype; /* float, int8, etc */

  hsize_t len; /* number of columns */
  hsize_t chunk; /* chunk size */
} H5VectorInfo;



typedef struct {
  long n_line; /* number of lines in file (including header lines) */
  long n_row;   /* number of data rows in file */
  long n_sample; /* number of samples in file */

  long n_geno_prob_col; /* number of genotype prob columns */
  long n_haplotype_col; /* number of haplotype columns */
} FileInfo;


void usage(char **argv) {
  fprintf(stderr, "\nusage: %s OPTIONS INPUT_FILE1 [INPUT_FILE2 ...]\n"
	  "\n"
	  "Description:\n"
	  "  This program converts VCF or IMPUTE files containing genetic\n"
	  "  polymorphism data into HDF5 files. HDF5 files are\n"
	  "  platform-independent binary files that provide very\n"
	  "  efficient programmatic access to data in C, Python (via\n"
	  "  PyTables) and several other programming languages\n"
	  "\n"
	  "Input Files:\n"
	  "     A separate VCF or IMPUTE input file must be provided\n"
	  "     for each chromosome. The filename must contain the name\n"
	  "     of the chromosome. Chromosome names should match those in\n"
	  "     the CHROM_FILE, which is provided by the --chrom option.\n"
	  "     IMPUTE input files are expected to end with .impute2.gz\n"
	  "     (for genotype probability files) or .impute2_haps.gz (for\n"
	  "     haplotype files)"
	  "\n"
	  "Input Options:\n"
	  "  --chrom CHROM_FILE [required]\n"
	  "     Path to chromInfo.txt file (may be gzipped) with list of\n"
	  "     chromosomes for the relevant genome assembly. Each line\n"
	  "     in file should contain tab-separated chromosome name and\n"
	  "     chromosome length (in basepairs). chromInfo.txt files can\n"
	  "     be downloaded from the UCSC genome browser. For example,\n"
	  "     a chromInfo.txt.gz file for hg19 can be downloaded from\n"
	  "     http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/\n"
	  "\n"
	  "  --format vcf|impute [required]\n"
	  "     Specifies the format of the input files. Currently supported\n"
	  "     formats are 'impute' or 'vcf'\n"
	  "\n"
	  "  --samples SAMPLES_FILE\n"
	  "     Input file containing ordered names of samples, one sample\n"
	  "     per line. This is only required for impute-formatted input\n"
	  "     files. Sample names for VCF input files are read from \n"
	  "     header lines in the VCF files.\n"
	  "\n"
	  "Output Options:\n"
	  "  --geno_prob GENO_PROB_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write genotype probabilities to.\n"
	  "     This option can only be used for impute2 files or VCF files\n"
	  "     that provide genotype likelihoods or posterior probabilities\n"
	  "     (GL or GP in the VCF FORMAT specifier).\n"
	  "\n"
	  "  --haplotype HAPLOTYPE_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write haplotypes to. This option can only\n"
	  "     be used if impute2_haps files or VCF files with genotypes (GT\n"
	  "     in the VCF FORMAT specifier) are provided. The genotypes\n"
	  "     should be phased (| haplotype separator). If the genotypes\n"
	  "     are unphased (/ haplotype separator) a warning is printed.\n"
	  "\n"
	  "  --snp_index SNP_INDEX_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write SNP index to. The SNP index can\n"
	  "     be used to convert the genomic position of a SNP to its\n"
	  "     corresponding row in the geno_prob, haplotype and snp_tab\n"
	  "     HDF5 files.\n"
	  "\n"
	  "  --snp_tab SNP_TABLE_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write SNP table to. Each row of SNP\n"
	  "     table contains SNP name (rs_id), position, allele1, allele2.\n"
	  "\n"
	  "Examples:\n"
	  "  # read 1000 genomes VCF files and write haplotype, snp_index \n"
	  "  # and snp_tab to HDF5 files\n"
	  "  snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \\\n"
	  "         --format vcf \\\n"
	  "         --haplotype haplotypes.h5 \\\n"
	  "         --snp_index snp_index.h5 \\\n"
	  "         --snp_tab   snp_tab.h5 \\\n"
	  "         data/1000G/ALL.chr*.vcf.gz\n"
	  "\n"
	  "  # read 1000 genomes VCF files that contain genotype likelihoods\n"
	  "  # and write genotype probabilties to HDF5 file\n"
	  "  snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \\\n"
	  "         --format vcf \\\n"
	  "         --geno_prob geno_probs.h5 \\\n"
	  "         1000G/supporting/genotype_likelihoods/shapeit2/ALL.chr*.gl.vcf.gz\n"
	  "\n"
	  "  # read genotype probabilities and haplotypes from IMPUTE2\n"
	  "  # output files and write to HDF5 files\n"
	  " snp2h5 --chrom data/ucsc/hg19/chromInfo.txt\n"
	  "    --format impute\n"
	  "    --geno_prob geno_probs.h5 \n"
	  "    --snp_index snp_index.h5 \n"
	  "    --snp_tab snp_tab.h5 \n"
	  "    --haplotype haps.h5 \n"
	  "    --samples samples_names.txt\n"
	  "    genotypes/chr*.hg19.impute2.gz \n"
	  "    genotypes/chr*.hg19.impute2_haps.gz\n"
	  "\n"
	  "\n", argv[0]);
}


void parse_args(Arguments *args, int argc, char **argv) {
  int c;
  char *format_str = NULL;
  
   static struct option loptions[] = {
     {"chrom",     required_argument, 0, 'c'},
     {"format",    required_argument, 0, 'f'},
     {"geno_prob", required_argument, 0, 'p'},
     {"haplotype", required_argument, 0, 'h'},
     {"snp_index", required_argument, 0, 'i'},
     {"snp_tab", required_argument, 0, 't'},
     {"samples", required_argument, 0, 's'},
     {0,0,0,0}
   };
   args->chrom_file = NULL;
   args->format = 0;
   args->geno_prob_file = NULL;
   args->haplotype_file = NULL;
   args->snp_index_file = NULL;
   args->snp_tab_file = NULL;
   args->sample_file = NULL;

   while(1) {
     c = getopt_long(argc, argv, "c:f:p:h:i:t:s:", loptions, NULL);
     
     if(c == -1) {
       break;
     }     
     switch (c) {
       /** TODO: fix this very small mem leak: filenames never freed */
       case 'c': args->chrom_file = util_str_dup(optarg); break;
       case 'f': format_str = util_str_dup(optarg); break;
       case 'p': args->geno_prob_file = util_str_dup(optarg); break;
       case 'h': args->haplotype_file = util_str_dup(optarg); break;
       case 'i': args->snp_index_file = util_str_dup(optarg); break;
       case 't': args->snp_tab_file = util_str_dup(optarg); break;
       case 's': args->sample_file = util_str_dup(optarg); break;
       default: usage(argv); break;
     }
   }
      
   if(optind < argc) {
     args->n_input_files = argc - optind;
     args->input_files = &argv[optind];
   } else {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "must specify at least one input file\n");
     exit(-1);
   }
   if(!args->chrom_file) {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "--chrom argument is missing\n");
     exit(-1);
   }
   if(!args->geno_prob_file && !args->haplotype_file &&
      !args->snp_index_file && !args->snp_tab_file) {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "at least one of --geno_prob, --haplotype, "
	     "--snp_index should be specified\n");
     exit(-1);
   }

   if(format_str == NULL) {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "--format argument is missing\n");
     exit(-1);
   } else {
     util_str_lc(format_str);
     if(strcmp(format_str, "vcf") == 0) {
       args->format = FORMAT_VCF;
     }
     else if(strcmp(format_str, "impute") == 0) {
       args->format = FORMAT_IMPUTE;
     }
     else {
       usage(argv);
       fprintf(stderr, "Error: ");
       fprintf(stderr, "--format should be one of 'vcf' or 'impute'");
       exit(-1);
     }
   }

   if(args->sample_file && (args->format == FORMAT_VCF)) {
     my_warn("ignoring sample names from --samples input file "
	     "because using sample information from VCF headers instead\n");
   }
   
}



Sample *read_sample_info(Arguments *args, int *n_sample) {
  gzFile gzf;
  Sample *samples;
  char *line;
  int i;

  *n_sample = 0;
  samples = NULL;
  
  if(args->sample_file) {
    /* read sample information from samples file */
    *n_sample = util_count_lines(args->sample_file);

    samples = my_malloc(sizeof(Sample) * *n_sample);

    gzf = util_must_gzopen(args->sample_file, "rb");
    i = 0;
    while((line = util_gzgets_line(gzf)) != NULL) {
      if(i >= *n_sample) {
	my_err("%s:%d: more sample lines than expected in file %s\n",
	       __FILE__, __LINE__, args->sample_file);
      }
      util_str_strip(line);
      util_strncpy(samples[i].name, line, sizeof(samples[i].name));
      my_free(line);
      i += 1;
    }
    if(i != *n_sample) {
      my_err("%s:%d: expected %d lines in file, but got %d\n",
	     __FILE__, __LINE__, *n_sample, i);
    }
  }

  if(args->format == FORMAT_VCF) {    
    /* Use sample info from VCF headers, rather than this input
     * file. Do this because the number of samples can differ across
     * VCF files (e.g. chrY VCF from 1000 genomes only has male
     * samples)
     */    
    my_warn("ignoring sample names from --samples input file "
	    "using sample information from VCF headers instead\n");
  }
  
  return samples;
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


void init_h5vector(H5VectorInfo *info, hsize_t len, hid_t datatype,
		   const char *name) {
  int rank = 1;
  herr_t status;
  hsize_t dim[1];
  hsize_t chunk[1];

  info->len = len;
  info->datatype = datatype;

  dim[0] = len;
  chunk[0] = (VEC_CHUNK < len) ? VEC_CHUNK : len;
    
  /* define a matrix dataspace which can hold genotype probs, etc... */

  info->file_dataspace = H5Screate_simple(rank, dim, NULL);
  if(info->file_dataspace < 0) {
    my_err("%s:%d: failed to create file dataspace", __FILE__,
	   __LINE__);
  }
    
  /* create property list for dataset */
  info->dataset_prop = H5Pcreate(H5P_DATASET_CREATE);
  if(info->dataset_prop < 0) {
    my_err("%s:%d: failed to create dataset property list",
	   __FILE__, __LINE__);
  }
  
  /* set chunking properties of dataset */
  status = H5Pset_chunk(info->dataset_prop, rank, chunk);
  if(status < 0) {
    my_err("%s:%d: failed to set chunksize", __FILE__, __LINE__);
  }
  /* use zlib compression level 6 */
  status = H5Pset_deflate(info->dataset_prop, 6);
  if(status < 0) {
    my_err("%s:%d: failed to set compression filter",
	   __FILE__, __LINE__);
  }
  
  /* create new dataset */
  info->dataset = H5Dcreate(info->h5file, name, info->datatype,
			    info->file_dataspace,
			    info->dataset_prop);
  if(info->dataset < 0) {
    my_err("%s:%d failed to create dataset\n", __FILE__, __LINE__);
  }
  
}


void init_h5matrix(H5MatrixInfo *info,
		   hsize_t n_row, hsize_t n_col,
		   hid_t datatype,
		   const char *name) {
  int rank = 2;
  herr_t status;

  fprintf(stderr, "initializing HDF5 matrix with dimension: "
	  "(%zu, %zu)\n", (size_t)n_row, (size_t)n_col);
  
  info->n_row = n_row;
  info->n_col = n_col;
  info->datatype = datatype;
  
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
  info->chunk[0] = (ROW_CHUNK < n_row) ? ROW_CHUNK : n_row;
  info->chunk[1] = (COL_CHUNK < n_col) ? COL_CHUNK : n_col;
  status = H5Pset_chunk(info->dataset_prop, rank, info->chunk);
  if(status < 0) {
    my_err("%s:%d: failed to set chunksize", __FILE__, __LINE__);
  }
  /* use zlib compression level 6 */
  status = H5Pset_deflate(info->dataset_prop, 6);
  if(status < 0) {
    my_err("%s:%d: failed to set compression filter",
	   __FILE__, __LINE__);
  }
  
  /* create new dataset */
  info->dataset = H5Dcreate(info->h5file, name, info->datatype,
			    info->file_dataspace,
			    info->dataset_prop);
  if(info->dataset < 0) {
    my_err("%s:%d failed to create dataset\n", __FILE__, __LINE__);
  }  
}



void write_h5matrix_row(H5MatrixInfo *info, hsize_t row_num,
			void *row_data) {
  hsize_t offset[2];
  hsize_t count[1];
  herr_t status;
  
  /* select a hyperslab that corresponds to a 
   * single row (all columns) 
   */
  offset[0] = row_num;
  offset[1] = 0;
  count[0] = 1;
  count[1] = info->n_col;
  status = H5Sselect_hyperslab(info->file_dataspace,
			       H5S_SELECT_SET, offset,
			       NULL, count, NULL);
      
  if(status < 0) {
    my_err("%s:%d: failed to select hyperslab for row %ld\n",
	   __FILE__, __LINE__, row_num);
  }

  /* write a row of data */
  status = H5Dwrite(info->dataset, info->datatype,
		    info->mem_dataspace,
		    info->file_dataspace, H5P_DEFAULT, row_data);
  
  if(status < 0) {
    my_err("%s:%d: failed to write data for row %ld",
	   __FILE__, __LINE__, row_num);
  }
}


void write_snp_index(H5VectorInfo *info, long *data) {
  herr_t status;
  
  status = H5Dwrite(info->dataset, info->datatype, H5S_ALL,
		    info->file_dataspace, H5P_DEFAULT, data);

  if(status < 0) {
    my_err("%s:%rd: failed to write data", __FILE__, __LINE__);
  }
}




void close_h5matrix(H5MatrixInfo *info) {
    H5Sclose(info->file_dataspace);
    H5Sclose(info->mem_dataspace);
    H5Dclose(info->dataset);
}

void close_h5vector(H5VectorInfo *info) {
  H5Sclose(info->file_dataspace);
  H5Dclose(info->dataset);
}



/**
 * Counts number of lines in file, reads header information,
 * sets attributes such as number of samples in FileInfo structure,
 * and advances file handle to end of header
 */
void set_file_info(gzFile gzf, char *filename, Arguments *args, FileInfo *file_info,
		   VCFInfo *vcf, ImputeInfo *impute_info) {
  long n_col;
  
  /* count total number lines in file, this tells us number of records */
  fprintf(stderr, "counting lines in file\n");
  file_info->n_line = util_count_lines(filename);
  fprintf(stderr, "  total lines: %ld\n", file_info->n_line);

  if(args->format == FORMAT_VCF) {
    /* parse VCF headers */
    fprintf(stderr, "reading VCF header\n");
    vcf_read_header(gzf, vcf);
    fprintf(stderr, "  VCF header lines: %ld\n", vcf->n_header_line);
    
    file_info->n_row = file_info->n_line - vcf->n_header_line;
    file_info->n_sample = vcf->n_sample;
  }
  else if(args->format == FORMAT_IMPUTE) {
    /* get number of samples from first row of IMPUTE file */
    file_info->n_row = file_info->n_line;
    n_col = impute_count_fields(gzf) - IMPUTE_FIX_HEADER;
    
    /* here we assume that the file is an IMPUTE file
     * with genotypes, and that info from this file is assumed
     * correct. Due to impute bug sometimes lines are missing 
     * from haplotype files.
     */

    if((n_col % 3) != 0) {
      my_err("%s:%d: expected number of columns (%ld) in .impute2 file to be "
	     "multiple of 3", __FILE__, __LINE__, n_col);
    }
    
    file_info->n_sample = n_col / 3;
    file_info->n_geno_prob_col = n_col;

    /* if we used haplotype file, n_sample would be n_col/2 */
    /* file_info->n_sample = n_col / 2; */
    
    impute_info->n_sample = file_info->n_sample;
    
  } else {
    my_err("%s:%d: unknown file format\n", __FILE__, __LINE__);
  }
  
  fprintf(stderr, "  number of samples: %ld\n", file_info->n_sample);    
  file_info->n_geno_prob_col = file_info->n_sample * 3;
  file_info->n_haplotype_col = file_info->n_sample * 2;  
}




/** groups impute files by whether their extension. Filenames with
 * extension .impute2.gz are put into imp_files array, while filenames
 * with extension .impute2_haps.gz are put into hap_files array.
 * files in the hap array are then ordered so that their chromosomal
 * ordering matches the imp array.
 * 
 */
void group_impute_files(char **all_files, int n_files, Chromosome *all_chrom,
			int n_chrom, char **imp_files, char **hap_files,
			int *n_imp_files, int *n_hap_files) {
  Chromosome *hap_chrom, *imp_chrom;
  char **tmp_hap_files;
  int i, j;

  *n_imp_files = 0;
  *n_hap_files = 0;

  tmp_hap_files = my_malloc(sizeof(char *) * n_files);
  
  for(i = 0; i < n_files; i++) {
    if(util_str_ends_with(all_files[i], ".impute2.gz")) {
      imp_files[*n_imp_files] = util_str_dup(all_files[i]);
      *n_imp_files += 1;
    }
    else if(util_str_ends_with(all_files[i], ".impute2_haps.gz")) {
      tmp_hap_files[*n_hap_files] = util_str_dup(all_files[i]);
      *n_hap_files += 1;
    }
    else {
      my_warn("%s:%d: ignoring file '%s'. Expected extension "
	     "'.impute2.gz' or 'impute2_haps.gz'", __FILE__, __LINE__,
	      all_files[i]);
    }
  }

  if(*n_hap_files > 0) {
    if(*n_hap_files != *n_imp_files) {
      my_err("%s:%d: expected same number of 'haps' and 'impute' files"
	     ", but got %d and %d", __FILE__, __LINE__, *n_imp_files,
	     *n_hap_files);
    }

    /* re-order hap files so chromosome order matches imp files
     */
    for(i = 0; i < *n_hap_files; i++) {
      imp_chrom = chrom_guess_from_file(imp_files[i],
					all_chrom, n_chrom);

      if(imp_chrom == NULL) {
	my_err("%s:%d: could not guess chromosome name from "
	       "filename %s\n", __FILE__, __LINE__, imp_files[i]);
      }
      
      hap_files[i] = NULL;
      
      for(j = 0; j < *n_imp_files; j++) {
	hap_chrom = chrom_guess_from_file(tmp_hap_files[j],
					  all_chrom, n_chrom);

	if(hap_chrom == NULL) {
	  my_err("%s:%d: could not guess chromosome name from "
		 "filename %s\n", __FILE__, __LINE__,
		 tmp_hap_files[i]);
	}
	
	if(strcmp(hap_chrom->name, imp_chrom->name) == 0) {
	  /* found match */
	  hap_files[i] = tmp_hap_files[j];
	  fprintf(stderr, "     %s pairs with %s\n", imp_files[i], tmp_hap_files[j]);
	  break;
	}
      }
      if(hap_files[i] == NULL) {
	fprintf(stderr, "     '%s' pairs with ???\n", imp_files[i]);
	my_err("%s:%d: could not find hap file for chrom %s\n",
	       __FILE__, __LINE__, imp_chrom->name);
      }
    }
  }

  my_free(tmp_hap_files);
}



void parse_impute(Arguments *args, Chromosome *all_chroms, int n_chrom,
		  H5MatrixInfo *gprob_info, H5MatrixInfo *haplotype_info,
		  H5VectorInfo *snp_index_info,
		  hid_t snp_tab_h5file) {
  ImputeInfo *impute_info;
  FileInfo file_info;
  long i, j;
  float *geno_probs;
  char *haplotypes;
  long *snp_index;
  SNP snp;
  SNPTab *snp_tab;
  hsize_t row;
  gzFile gzf;
  Chromosome *chrom;		  
  char **imp_files;
  char **hap_files;
  int n_imp_files;
  int n_hap_files;
  int n_sample;
  long missing_geno_probs;
  long n_haplotype_row;
  SampleTab *samp_tab;
  Sample *samples;

  if(args->sample_file) {
    samples = read_sample_info(args, &n_sample);
  } else {
    samples = NULL;
    n_sample = 0;
  }

  impute_info = impute_info_new();
  
  imp_files = my_malloc(sizeof(char *) * args->n_input_files);
  hap_files = my_malloc(sizeof(char *) * args->n_input_files);

  group_impute_files(args->input_files, args->n_input_files,
		     all_chroms, n_chrom, imp_files, hap_files,
		     &n_imp_files, &n_hap_files);

  /* loop over input files, there should be one for each chromosome */
  for(i = 0; i < n_imp_files; i++) {
    /* guess chromosome file filename */
    chrom = chrom_guess_from_file(imp_files[i], all_chroms,
				  n_chrom);
    if(chrom == NULL) {
      my_err("%s:%d: could not guess chromosome from filename '%s'\n",
	     __FILE__, __LINE__, args->input_files[i]);
    }

    fprintf(stderr, "chromosome: %s, length: %ldbp\n",
	    chrom->name, chrom->len);

    /* open file, count number of lines */
    fprintf(stderr, "reading from file %s\n", imp_files[i]);
    gzf = util_must_gzopen(args->input_files[i], "rb");
    set_file_info(gzf, args->input_files[i], args, &file_info, NULL,
		  impute_info);
    
    /* initialize output files and memory to hold genotypes,
     * haplotypes, etc
     */
    if(args->geno_prob_file) {
      geno_probs = my_malloc(file_info.n_geno_prob_col*sizeof(float));
      init_h5matrix(gprob_info, file_info.n_row,
		    file_info.n_geno_prob_col,
		    GENO_PROB_DATATYPE, chrom->name);

      if(samples) {
	if((n_sample*3) != file_info.n_geno_prob_col) {
	  my_warn("%s:%d number of samples*3 (%d*3=%d) does not match "
		  "number of genotype columns for chromosome %s "
		  "(%d)\n",
		  __FILE__, __LINE__, n_sample,  n_sample*3, chrom->name,
		  file_info.n_geno_prob_col);
	}
	
	/* write genotype prob sample names table for this chromosome */
	samp_tab = sample_tab_create(gprob_info->h5file, chrom->name,
				     samples, n_sample);
	sample_tab_free(samp_tab);
      }

      /* fill H5Matrix with default values for genotype probabilities */
      for(j = 0; j < file_info.n_geno_prob_col; j++) {
	geno_probs[j] = GENO_PROB_DEFAULT_VAL;
      }
      for(j = 0; j < file_info.n_row; j++) {
	write_h5matrix_row(gprob_info, j, geno_probs);
      }

    } else {
      geno_probs = NULL;
    }

    /* SNP index vector is same length as chromosome,
     * initialize to SNP_INDEX_NONE (-1)
     */
    snp_index = my_malloc(chrom->len * sizeof(long));
    for(j = 0; j < chrom->len; j++) {
      snp_index[j] = SNP_INDEX_NONE;
    }

    if(args->snp_index_file) {
      init_h5vector(snp_index_info, chrom->len,
		    SNP_INDEX_DATATYPE, chrom->name);      
    }

    if(args->snp_tab_file) {
      snp_tab = snp_tab_new(snp_tab_h5file, chrom->name,
			    file_info.n_row);      
    } else {
      snp_tab = NULL;
    }
    
    row = 0;
    
    fprintf(stderr, "parsing file and writing to HDF5 files\n");
    
    while(impute_read_line(gzf, impute_info,
			   &snp, geno_probs, NULL) != -1) {

      if(geno_probs) {
	write_h5matrix_row(gprob_info, row, geno_probs);
      }

      /*  set snp_index element at this chromosome position
       * to point to row in matrices / SNP table 
       */
      if(snp.pos > chrom->len || snp.pos < 1) {
	my_err("%s:%d: SNP position (%ld) is outside of "
	       "chromomosome %s range:1-%ld", __FILE__, __LINE__,
	       snp.pos, chrom->len);
      }

      /* set value in snp_index array to point to current row */
      snp_index[snp.pos-1] = row;

      /* append row to SNP table */
      if(snp_tab) {
	snp_tab_append_row(snp_tab, &snp);
      }
      
      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }
    if(row != file_info.n_row) {
      my_warn("%s:%d: expected %ld data rows, but only read %ld\n",
	      __FILE__, __LINE__, file_info.n_row, row);
    }
    gzclose(gzf);

    /* clean up geno_prob memory for this chromosome */
    if(geno_probs) {
      my_free(geno_probs);
      close_h5matrix(gprob_info);
    }

    fprintf(stderr, "\n");
    
    if(args->haplotype_file) {
      /* NOW do pass of haplotype file. */
      /* Due to bug? in IMPUTE some SNPs are missing from the haplotype
       * files, so we need to lookup correct index in hap table
       * using the snp_index array we created from the main impute
       * files.
       */
	haplotypes = my_malloc(file_info.n_haplotype_col *
			       sizeof(char));
	init_h5matrix(haplotype_info, file_info.n_row,
		      file_info.n_haplotype_col,
		      HAPLOTYPE_DATATYPE, chrom->name);

	if(samples) {
	  if((n_sample*2) != file_info.n_haplotype_col) {
	    my_warn("%s:%d number of samples*2 (%d*2=%d) does not match "
		    "number of haplotype columns for chromosome %s "
		    "(%d)\n", __FILE__, __LINE__, n_sample,
		    n_sample*2, chrom->name, file_info.n_haplotype_col);
	  }
	
	  /* write haplotype sample names table for this chromosome */
	  samp_tab = sample_tab_create(haplotype_info->h5file, chrom->name,
				       samples, n_sample);
	  sample_tab_free(samp_tab);
	}

	
	/* fill H5Matrix with default values for haplotypes */
	for(j = 0; j < file_info.n_haplotype_col; j++) {
	  haplotypes[j] = HAPLOTYPE_DEFAULT_VAL;
	}
	for(j = 0; j < file_info.n_row; j++) {
	  write_h5matrix_row(haplotype_info, j, haplotypes);
	}

	fprintf(stderr, "reading from file %s\n", hap_files[i]);
	gzf = util_must_gzopen(hap_files[i], "rb");
      
	/* set gzip buffer size to 128k (131072 bytes) */
	/* if(gzbuffer(gzf, 131072) == -1) { */
	/*   my_err("%s:%d: failed to set gzbuffer size\n", */
	/* 	 __FILE__, __LINE__); */
	/* } */
	  
	fprintf(stderr, "parsing file and writing to HDF5 files\n");

	int line_num = 0;
	n_haplotype_row = 0;
	missing_geno_probs = 0;
	while(impute_read_line(gzf, impute_info, &snp,
			       NULL, haplotypes) != -1) {

	  if(snp.pos > chrom->len || snp.pos < 1) {
	    my_err("%s:%d: SNP position (%ld) is outside of "
		   "chromomosome %s range:1-%ld", __FILE__, __LINE__,
		   snp.pos, chrom->len);
	  }
	    
	  row = snp_index[snp.pos-1];
	  if(row == SNP_INDEX_NONE) {
	    /* ignore SNPs that are present in haps file but missing from
	     * genotype file. 
	     */
	    
	    if(!missing_geno_probs) {
	      my_warn("%s:%d: SNP %s as position %ld is present in "
		      "impute2_haps file but not in impute2 file\n",
		      __FILE__, __LINE__, snp.name, snp.pos);
	    }

	    missing_geno_probs += 1;
	  } else {
	    write_h5matrix_row(haplotype_info, row, haplotypes);
	    n_haplotype_row += 1;
	  }

	  line_num += 1;

	  if((line_num % 10000) == 0) {
	    fprintf(stderr, ".");
	  }
	}
	
	fprintf(stderr, "\n");
	if(missing_geno_probs) {
	  /* there were missing entries in geno_probs file */
	  my_warn("%s:%d: %ld SNPs were present in "
		  "impute2_haps file but not in impute2 file\n",
		  __FILE__, __LINE__, missing_geno_probs);
	}

	if(n_haplotype_row < file_info.n_row) {
	  /* there were missing entries in haplotype file */
	  my_warn("%s:%d: %ld SNPs were present in "
		  "impute file but not in impute2_haps file\n",
		  __FILE__, __LINE__, file_info.n_row - n_haplotype_row);
	}
	if(n_haplotype_row > file_info.n_row) {
	  my_warn("%s:%d: there were duplicate SNP entries in "
		  "impute2_haps file\n", __FILE__, __LINE__);
	}
	
	
	gzclose(gzf);

	/* cleanup haplotype memory for this chrom */
	my_free(haplotypes);
	close_h5matrix(haplotype_info);
    }

    fprintf(stderr, "\n");
    
    /* write snp_index data */
    if(args->snp_index_file) {
      write_snp_index(snp_index_info, snp_index);
      close_h5vector(snp_index_info);
    }

    /* cleanup snp_index and snp_tab memory for this chromosome */
    my_free(snp_index);
    
    if(snp_tab) {
      snp_tab_free(snp_tab);
    }
  }
  
  
  impute_info_free(impute_info);
}






void parse_vcf(Arguments *args, Chromosome *all_chroms, int n_chrom,
	       H5MatrixInfo *gprob_info, H5MatrixInfo *haplotype_info,
	       H5VectorInfo *snp_index_info, hid_t snp_tab_h5file) {
  VCFInfo *vcf;
  FileInfo file_info;
  SNPTab *snp_tab;
  long i, j;
  float *geno_probs;
  char *haplotypes;
  long *snp_index;
  SNP snp;
  hsize_t row;
  gzFile gzf;
  Chromosome *chrom;
  SampleTab *samp_tab;
  
  vcf = vcf_info_new();

  /* loop over input files, there should be one for each chromosome */
  for(i = 0; i < args->n_input_files; i++) {
    /* guess chromosome file filename */
    chrom = chrom_guess_from_file(args->input_files[i],
				  all_chroms, n_chrom);
    if(chrom == NULL) {
      my_err("%s:%d: could not guess chromosome from filename '%s'\n",
	     __FILE__, __LINE__, args->input_files[i]);
    }

    fprintf(stderr, "chromosome: %s, length: %ldbp\n",
	    chrom->name, chrom->len);

    fprintf(stderr, "reading from file %s\n", args->input_files[i]);
    gzf = util_must_gzopen(args->input_files[i], "rb");

    set_file_info(gzf, args->input_files[i], args, &file_info, vcf, NULL);
    
    /* initialize output files and memory to hold genotypes,
     * haplotypes, etc
     */
    if(args->geno_prob_file) {
      if(vcf->n_sample == 0) {
	my_warn("VCF file has 0 samples, cannot write genotype probabilities");
	geno_probs = NULL;
      } else {
	geno_probs = my_malloc(file_info.n_geno_prob_col *
			       sizeof(float));
	init_h5matrix(gprob_info, file_info.n_row,
		      file_info.n_geno_prob_col,
		      GENO_PROB_DATATYPE, chrom->name);

	if((vcf->n_sample*3) != file_info.n_geno_prob_col) {
	  my_warn("%s:%d number of samples*3 (%d*3=%d) does not match "
		  "number of genotype columns for chromosome %s "
		  "(%d)\n", __FILE__, __LINE__, vcf->n_sample,
		  vcf->n_sample*3, chrom->name, file_info.n_geno_prob_col);
	}
      
	/* create table of sample names for genotype probs */
	samp_tab = sample_tab_from_names(gprob_info->h5file, chrom->name,
					 vcf->sample_names, vcf->n_sample);
	sample_tab_free(samp_tab);
      }
    } else {
      geno_probs = NULL;
    }
    if(args->haplotype_file) {
      if(vcf->n_sample == 0) {
	my_warn("VCF file has 0 samples, cannot write haplotypes");
	haplotypes = NULL;
      } else {
	haplotypes = my_malloc(file_info.n_haplotype_col * sizeof(char));
	init_h5matrix(haplotype_info, file_info.n_row,
		      file_info.n_haplotype_col,
		      HAPLOTYPE_DATATYPE, chrom->name);

	if((vcf->n_sample*2) != file_info.n_haplotype_col) {
	  my_warn("%s:%d number of samples (%d*2=%d) does not match "
		  "number of haplotype columns for chromosome %s "
		  "(%d)\n", __FILE__, __LINE__, vcf->n_sample,
		  vcf->n_sample*2, chrom->name, file_info.n_haplotype_col);
	}
      
	/* create table of sample names for haplotypes */
	samp_tab = sample_tab_from_names(haplotype_info->h5file, chrom->name,
					 vcf->sample_names, vcf->n_sample);
	sample_tab_free(samp_tab);
      }
    } else {
      haplotypes = NULL;
    }
    if(args->snp_index_file) {
      /* SNP index vector is same length as chromosome,
       * initialize to SNP_INDEX_NONE (-1)
       */
      snp_index = my_malloc(chrom->len * sizeof(long));
      for(j = 0; j < chrom->len; j++) {
	snp_index[j] = SNP_INDEX_NONE;
      }
      init_h5vector(snp_index_info, chrom->len,
		    SNP_INDEX_DATATYPE, chrom->name);
    } else {
      snp_index = NULL;
    }
    if(args->snp_tab_file) {
      snp_tab = snp_tab_new(snp_tab_h5file, chrom->name,
			    file_info.n_row);
    } else {
      snp_tab = NULL;
    }

    row = 0;

    fprintf(stderr, "parsing file and writing to HDF5 files\n");

    while(vcf_read_line(gzf, vcf, &snp,
			geno_probs, haplotypes) != -1) {
      
      if(geno_probs) {
	write_h5matrix_row(gprob_info, row, geno_probs);
      }
      if(haplotypes) {
	write_h5matrix_row(haplotype_info, row, haplotypes);
      }

      /*  set snp_index element at this chromosome position
       * to point to row in matrices / SNP table 
       */
      if(snp.pos > chrom->len || snp.pos < 1) {
	my_err("%s:%d: SNP %s position (%ld) is outside of "
	       "chromomosome %s range: 1-%ld", __FILE__, __LINE__,
	       snp.name, snp.pos, chrom->name, chrom->len);
      }

      if(snp_index) {
	/* set value in snp_index array to point to current row */
	snp_index[snp.pos-1] = row;
      }

      /* append row to SNP table */
      if(snp_tab) {
	snp_tab_append_row(snp_tab, &snp);
      }

      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }
    fprintf(stderr, "\n");

    if(row != file_info.n_row) {
      my_warn("%s:%d: expected %ld data rows, but only read %ld\n",
	      __FILE__, __LINE__, file_info.n_row, row);
    }

    /* write snp_index data */
    if(snp_index) {
      write_snp_index(snp_index_info, snp_index);
    }
        
    /* clean up memory etc. for this chromosome */
    if(geno_probs) {
      my_free(geno_probs);
      close_h5matrix(gprob_info);
    }
    if(haplotypes) {
      my_free(haplotypes);
      close_h5matrix(haplotype_info);
    }
    if(snp_index) {
      my_free(snp_index);
      close_h5vector(snp_index_info);
    }
    if(snp_tab) {
      snp_tab_free(snp_tab);
    }
    gzclose(gzf);
  }

  vcf_info_free(vcf);
}





int main(int argc, char **argv) {
  Arguments args;
  Chromosome *all_chroms;
  int n_chrom;
  H5MatrixInfo gprob_info, haplotype_info;
  H5VectorInfo snp_index_info;
  hid_t snp_tab_h5file;

  parse_args(&args, argc, argv);

  /* read chromosomes */
  all_chroms = chrom_read_file(args.chrom_file, &n_chrom);

  fprintf(stderr, "long alleles will be truncated to %dbp\n", SNP_MAX_ALLELE);

  /* create new HDF5 file(s) */
  if(args.geno_prob_file) {
    gprob_info.h5file = create_h5file(args.geno_prob_file);
    fprintf(stderr, "writing genotype probabilities to: %s\n",
	    args.geno_prob_file);
  }
  if(args.haplotype_file) {
    haplotype_info.h5file = create_h5file(args.haplotype_file);
    fprintf(stderr, "writing haplotypes to: %s\n",
	    args.haplotype_file);
  }
  if(args.snp_index_file) {
    snp_index_info.h5file = create_h5file(args.snp_index_file);
    fprintf(stderr, "writing SNP index to: %s\n", args.snp_index_file);
  }
  if(args.snp_tab_file) {
    snp_tab_h5file = create_h5file(args.snp_tab_file);
    fprintf(stderr, "writing SNP table to: %s\n", args.snp_tab_file);
  }
    
  
  if(args.format == FORMAT_VCF) {
    parse_vcf(&args, all_chroms, n_chrom,
	      &gprob_info, &haplotype_info, &snp_index_info,
	      snp_tab_h5file);
  }
  else if(args.format == FORMAT_IMPUTE) {
    parse_impute(&args, all_chroms, n_chrom,
		 &gprob_info, &haplotype_info, &snp_index_info,
		 snp_tab_h5file);
  }
  else {
    my_err("%s:%d: unknown format\n", __FILE__, __LINE__);
  }
    
  chrom_array_free(all_chroms, n_chrom);

  /* TODO: fix small mem leak: sample names never freed */
  
  /* close HDF5 files */
  if(args.geno_prob_file) {
    H5Fclose(gprob_info.h5file);
  }
  if(args.haplotype_file) {
    H5Fclose(haplotype_info.h5file);
  }
  if(args.snp_index_file) {
    H5Fclose(snp_index_info.h5file);
  }
  if(args.snp_tab_file) {
    H5Fclose(snp_tab_h5file);
  }

  fprintf(stderr, "done\n");
  
  return 0;
}
