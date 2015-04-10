
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <hdf5.h>
#include <stdint.h>

#include "vcf.h"
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
#define HAPLOTYPE_DATATYPE H5T_NATIVE_SHORT
#define SNP_INDEX_DATATYPE H5T_NATIVE_LONG

#define SNP_INDEX_NONE -1

typedef struct {
  /* chromInfo.txt.gz file that chromosomes are read from */
  char *chrom_file;

  /* HDF5 file that genotype probabilities are written to */
  char *geno_prob_file;
  char *haplotype_file;
  char *snp_index_file;
  char *snp_tab_file;
  char **vcf_files;
  
  int n_vcf_files;
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


void usage(char **argv) {
  fprintf(stderr, "\nusage: %s OPTIONS VCF1 [VCF2 ...]\n"
	  "\n"
	  "Description:\n"
	  "  This program converts VCF files containing genetic\n"
	  "  polymorphism data into HDF5 files. HDF5 files are\n"
	  "  platform-independent binary files that provide very\n"
	  "  efficient programmatic access to data in C, Python (via\n"
	  "  PyTables) and several other programming languages\n"
	  "\n"
	  "VCF Input Files:\n"
	  "     A separate VCF input file must be provided for each\n"
	  "     chromosome. The filename must contain the name of the\n"
	  "     chromosome. Chromosome names should match those in the\n"
	  "     CHROM_FILE, which is provided by the --chrom option.\n"
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
	  "Output Options:\n"
	  "  --geno_prob GENO_PROB_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write genotype probabilities to.\n"
	  "     This option can only be used if the input VCF files provide\n"
	  "     genotype likelihoods, (GL in the FORMAT specifier).\n"
	  "\n"
	  "  --haplotype HAPLOTYPE_OUTPUT_FILE [optional]\n"
	  "     Path to HDF5 file to write haplotypes to.\n"
	  "     This option can only be used if the input VCF files provide\n"
	  "     genotypes (GT in the VCF FORMAT specifier). The genotypes\n"
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
	  "  vcf2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \\\n"
	  "         --haplotype haplotypes.h5 \\\n"
	  "         --snp_index snp_index.h5 \\\n"
	  "         --snp_tab   snp_tab.h5 \\\n"
	  "         data/1000G/ALL.chr*.vcf.gz\n"
	  "\n"
	  "  # read 1000 genomes VCF files that contain genotype likelihoods\n"
	  "  # and write genotype probabilties to HDF5 file\n"
	  "  vcf2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \\\n"
	  "         --geno_prob geno_probs.h5 \\\n"
	  "         1000G/supporting/genotype_likelihoods/shapeit2/ALL.chr*.gl.vcf.gz\n"
	  "\n"
	  "\n", argv[0], argv[0], argv[0]);
  exit(-1);
}


void parse_args(Arguments *args, int argc, char **argv) {
  int c, i;
  
   static struct option loptions[] = {
     {"chrom",     required_argument, 0, 'c'},
     {"geno_prob", required_argument, 0, 'p'},
     {"haplotype", required_argument, 0, 'h'},
     {"snp_index", required_argument, 0, 'i'},
     {"snp_tab", required_argument, 0, 't'},
     {0,0,0,0}
   };
   args->chrom_file = NULL;
   args->geno_prob_file = NULL;
   args->haplotype_file = NULL;
   args->snp_index_file = NULL;
   args->snp_tab_file = NULL;

   while(1) {
     c = getopt_long(argc, argv, "c:p:h:i:t:", loptions, NULL);
     
     if(c == -1) {
       break;
     }     
     switch (c) {
       /** TODO: fix this very small mem leak: filenames never freed */
       case 'c': args->chrom_file = util_str_dup(optarg); break;
       case 'p': args->geno_prob_file = util_str_dup(optarg); break;
       case 'h': args->haplotype_file = util_str_dup(optarg); break;
       case 'i': args->snp_index_file = util_str_dup(optarg); break;
       case 't': args->snp_tab_file = util_str_dup(optarg); break;
       default: usage(argv); break;
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
   if(!args->geno_prob_file && !args->haplotype_file &&
      !args->snp_index_file && !args->snp_tab_file) {
     fprintf(stderr, "at least one of --geno_prob, --haplotype, "
	     "--snp_index should be specified\n");
     usage(argv);
   }
}



/**
 * Returns appropriate chromosome for filename by looking for match
 * with chromosome name in filename. Returns NULL if no match found.
 */
Chromosome *guess_chrom(char *filename, Chromosome *chroms,
			int n_chrom) {
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


void init_h5vector(H5VectorInfo *info, hsize_t len, hid_t datatype,
		   const char *name) {
  int rank = 1;
  herr_t status;
  hsize_t dim[1];
  hsize_t chunk[1];

  info->len = len;
  info->datatype = datatype;

  dim[0] = len;
  chunk[0] = VEC_CHUNK;
    
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
  info->chunk[0] = ROW_CHUNK;
  info->chunk[1] = COL_CHUNK;
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



int main(int argc, char **argv) {
  VCFInfo vcf;
  Arguments args;
  Chromosome *all_chroms, *chrom;
  int n_chrom;
  H5MatrixInfo gprob_info, haplotype_info;
  H5VectorInfo snp_index_info;
  SNPTab *snp_tab;
  hsize_t n_geno_prob_col, n_haplotype_col;
  hsize_t n_row, row;
  hid_t snp_tab_h5file;
  gzFile gzf;
  float *geno_probs;
  char *haplotypes;
  long *snp_index;
  size_t n_lines;
  long i, j;

  parse_args(&args, argc, argv);

  /* read chromosomes */
  all_chroms = chrom_read_file(args.chrom_file, &n_chrom);
  
  /* create new HDF5 file(s) */
  if(args.geno_prob_file) {
    gprob_info.h5file = create_h5file(args.geno_prob_file);
    fprintf(stderr, "writing genotype probabilities to: %s\n",
	    args.geno_prob_file);
  }
  if(args.haplotype_file) {
    haplotype_info.h5file = create_h5file(args.haplotype_file);
    fprintf(stderr, "writing haplotypes to: %s\n", args.haplotype_file);
  }
  if(args.snp_index_file) {
    snp_index_info.h5file = create_h5file(args.snp_index_file);
    fprintf(stderr, "writing SNP index to: %s\n", args.snp_index_file);
  }
  if(args.snp_tab_file) {
    snp_tab_h5file = create_h5file(args.snp_tab_file);
    fprintf(stderr, "writing SNP table to: %s\n", args.snp_tab_file);
  }
    
  for(i = 0; i < args.n_vcf_files; i++) {
    chrom = guess_chrom(args.vcf_files[i], all_chroms, n_chrom);
    if(chrom == NULL) {
      my_err("%s:%d: could not guess chromosome from filename '%s'\n",
	     __FILE__, __LINE__, args.vcf_files[i]);
    }
    fprintf(stderr, "chromosome: %s, length: %ldbp\n",
	    chrom->name, chrom->len);
    
    fprintf(stderr, "reading from VCF file %s\n", args.vcf_files[i]);
    gzf = util_must_gzopen(args.vcf_files[i], "rb");
    
    fprintf(stderr, "counting lines in file\n");
    n_lines = util_gzcount_lines(gzf);
    fprintf(stderr, "  total lines: %ld\n", n_lines);

    fprintf(stderr, "reading VCF header\n");
    vcf_read_header(gzf, &vcf);
    fprintf(stderr, "  VCF header lines: %ld\n", vcf.n_header_lines);
    fprintf(stderr, "  VCF samples: %ld\n", vcf.n_samples);

    
    n_row = n_lines - vcf.n_header_lines;
    
    n_geno_prob_col = vcf.n_samples * 3;
    n_haplotype_col = vcf.n_samples * 2;

    if(args.geno_prob_file) {
      geno_probs = my_malloc(n_geno_prob_col * sizeof(float)); 
      init_h5matrix(&gprob_info, n_row, n_geno_prob_col,
		    GENO_PROB_DATATYPE, chrom->name);
    } else {
      geno_probs = NULL;
    }
    if(args.haplotype_file) {
      haplotypes = my_malloc(n_haplotype_col * sizeof(char));
      init_h5matrix(&haplotype_info, n_row, n_haplotype_col,
		    HAPLOTYPE_DATATYPE, chrom->name);
    } else {
      haplotypes = NULL;
    }
    if(args.snp_index_file) {
      /* SNP index vector is same length as chromosome,
       * initialize to SNP_INDEX_NONE (-1)
       */
      snp_index = my_malloc(chrom->len * sizeof(long));
      for(j = 0; j < chrom->len; j++) {
	snp_index[j] = SNP_INDEX_NONE;
      }
      init_h5vector(&snp_index_info, chrom->len,
		    SNP_INDEX_DATATYPE, chrom->name);
    } else {
      snp_index = NULL;
    }
    if(args.snp_tab_file) {
      snp_tab = snp_tab_new(snp_tab_h5file, chrom->name);
    } else {
      snp_tab = NULL;
    }    
    
    row = 0;

    fprintf(stderr, "parsing VCF and writing to HDF5 files\n");
    
    while(vcf_read_line(gzf, &vcf, geno_probs, haplotypes) != -1) {
      
       /* fprintf(stderr, "chrom: %s, id: %s, pos: %ld, ref: %s, alt: %s " */
       /* 	      "qual: %s, filter: %s, info: %s format: %s\n", */
       /* 	      vcf.chrom, vcf.id, vcf.pos, vcf.ref_allele, vcf.alt_allele, */
       /* 	      vcf.qual, vcf.filter, vcf.info, vcf.format); */
      if(geno_probs) {
	write_h5matrix_row(&gprob_info, row, geno_probs);
      }
      if(haplotypes) {
	write_h5matrix_row(&haplotype_info, row, haplotypes);
      }

      /*  set snp_index element at this chromosome position
       * to point to row in matrices / SNP table 
       */
      if(vcf.snp.pos > chrom->len || vcf.snp.pos < 1) {
	my_err("%s:%d: SNP position (%ld) is outside of "
	       "chromomosome %s range:1-%ld", __FILE__, __LINE__,
	       vcf.snp.pos, chrom->len);
      }

      if(snp_index) {
	snp_index[vcf.snp.pos-1] = row;
      }

      /* append row to SNP table */
      if(snp_tab) {
	snp_tab_append_row(snp_tab, &vcf.snp);
      }
      
      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }
    fprintf(stderr, "\n");

    if(row != n_row) {
      my_warn("%s:%d: expected %ld data rows, but only read %ld\n",
	      __FILE__, __LINE__, n_row, row);
    }

    /* write snp_index data */
    if(snp_index) {
      write_snp_index(&snp_index_info, snp_index);
    }
        
    /* clean up memory etc. for this chromosome */
    if(geno_probs) {
      my_free(geno_probs);
      close_h5matrix(&gprob_info);
    }
    if(haplotypes) {
      my_free(haplotypes);
      close_h5matrix(&haplotype_info);
    }
    if(snp_index) {
      my_free(snp_index);
      close_h5vector(&snp_index_info);
    }
    if(snp_tab) {
      snp_tab_free(snp_tab);
    }
    gzclose(gzf);
  }
  
  chrom_array_free(all_chroms, n_chrom);

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
