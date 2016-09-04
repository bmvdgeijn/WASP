
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <hdf5.h>
#include <stdint.h>

#include "seq.h"
#include "util.h"
#include "memutil.h"
#include "chrom.h"


/* 
 * These values control size of chunks that are used for 
 * data transfer and compression. Changing them
 * affects effciency of reading/writing time and compression.
 * Not sure what optimal values are...
 */
#define VEC_CHUNK 100000

#define SEQ_DATATYPE H5T_NATIVE_CHAR


typedef struct {
  /* chromInfo.txt.gz file that chromosomes are read from */
  char *chrom_file;

  /* HDF5 file that sequence is written to */
  char *seq_file;
  char **input_files;
  
  int n_input_files;
} Arguments;



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
  fprintf(stderr, "\nusage: %s OPTIONS FASTA_FILE1 [FASTA_FILE2 ...]\n"
	  "\n"
	  "Description:\n"
	  "  This program converts FASTA files into HDF5 files of\n"
	  "  genome sequence\n"
	  "\n"
	  "Input Files:\n"
	  "     A separate FASTA file must be provided\n"
	  "     for each chromosome. The filename must contain the\n"
	  "     name of the chromosome. Chromosome names should\n"
	  "     match those in the CHROM_FILE, which is provided by\n"
	  "     the --chrom option.\n"
	  "\n"
	  "Input Options:\n"
	  "  --chrom CHROM_FILE [required]\n"
	  "     Path to chromInfo.txt file (may be gzipped) with\n"
	  "     list of chromosomes for the relevant genome\n"
	  "     assembly. Each line in file should contain\n"
	  "     tab-separated chromosome name and chromosome length\n"
	  "     (in basepairs). chromInfo.txt files can be\n"
	  "     downloaded from the UCSC genome browser. For\n"
	  "     example, a chromInfo.txt.gz file for hg19 can be\n"
	  "     downloaded from\n"
	  "     http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/\n"
	  "\n"
	  "Output Options:\n"
	  "  --seq SEQ_OUTPUT_FILE [required]\n"
	  "     Path to output HDF5 file to write sequence to\n"
	  "\n"
	  "Example:\n"
	  "  # read FASTA files and write HDF5 file\n"
	  "  fasta2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \\\n"
	  "         --seq seq.h5 \\\n"
	  "         data/ucsc/hg19/chr*.fasta.gz\n"
	  "\n"
	  "\n", argv[0]);
}


void parse_args(Arguments *args, int argc, char **argv) {
  int c;
  
   static struct option loptions[] = {
     {"chrom",     required_argument, 0, 'c'},
     {"seq",    required_argument, 0, 's'},
     {0,0,0,0}
   };
   args->chrom_file = NULL;
   args->seq_file = NULL;

   while(1) {
     c = getopt_long(argc, argv, "c:s:", loptions, NULL);
     
     if(c == -1) {
       break;
     }     
     switch (c) {
       /** TODO: fix this very small mem leak: filenames never freed */
       case 'c': args->chrom_file = util_str_dup(optarg); break;
       case 's': args->seq_file = util_str_dup(optarg); break;
       default: usage(argv); break;
     }
   }
      
   if(optind < argc) {
     args->n_input_files = argc - optind;
     args->input_files = &argv[optind];
   } else {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "must specify at least one FASTA input file\n");
     exit(-1);
   }
   if(!args->chrom_file) {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "--chrom argument is missing\n");
     exit(-1);
   }
   if(!args->seq_file) {
     usage(argv);
     fprintf(stderr, "Error: ");
     fprintf(stderr, "--seq argument is missing\n");
     exit(-1);
   }
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
    
  /* define a dataspace which can hold vector of data etc... */
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



void write_seq(H5VectorInfo *info, char *seq_str) {
  herr_t status;
  
  status = H5Dwrite(info->dataset, info->datatype, H5S_ALL,
		    info->file_dataspace, H5P_DEFAULT, seq_str);

  if(status < 0) {
    my_err("%s:%rd: failed to write data", __FILE__, __LINE__);
  }
}



void close_h5vector(H5VectorInfo *info) {
  H5Sclose(info->file_dataspace);
  H5Dclose(info->dataset);
}




void parse_fasta(Arguments *args, H5VectorInfo *seq_vec_info) {
  Seq *seq;
  Chromosome *chrom, *all_chroms;
  char *seq_str;
  int i, n_chrom;

  seq = seq_new();
  
  /* read chromosomes */
  all_chroms = chrom_read_file(args->chrom_file, &n_chrom);

  for(i = 0; i < args->n_input_files; i++) {
    chrom = chrom_guess_from_file(args->input_files[i],
				  all_chroms, n_chrom);
        
    if(chrom == NULL) {
      my_err("%s:%d: could not guess chromosome from filename "
	     "%s\n", __FILE__, __LINE__, args->input_files[i]);
    }

    fprintf(stderr, "chromosome: %s, length: %ldbp\n",
	    chrom->name, chrom->len);

    /* seq sequence from fasta file */
    seq_read_fasta_from_file(seq, args->input_files[i]);
    
    if(chrom->len != seq->len) {
      my_err("%s:%d: chromosome length %ld does not match sequence "
	     "length %ld", __FILE__, __LINE__, chrom->len, seq->len);
    }

    init_h5vector(seq_vec_info, chrom->len, SEQ_DATATYPE,
		  chrom->name);

    seq_str = seq_get_seqstr(seq);
    
    /* write sequence to HDF5 */
    fprintf(stderr, "writing to HDF5 file\n");
    write_seq(seq_vec_info, seq_str);

    my_free(seq_str);

    close_h5vector(seq_vec_info);
  }

  seq_free(seq);
  
  chrom_array_free(all_chroms, n_chrom);
}



int main(int argc, char **argv) {
  Arguments args;
  H5VectorInfo seq_vec_info;
  
  parse_args(&args, argc, argv);
  
  /* create new HDF5 file */
  seq_vec_info.h5file = create_h5file(args.seq_file);
  fprintf(stderr, "writing sequence to: %s\n", args.seq_file);

  /* read sequences from FASTA files and write to HDF5 file */
  parse_fasta(&args, &seq_vec_info);
  
  /* close HDF5 file */
  H5Fclose(seq_vec_info.h5file);

  fprintf(stderr, "done\n");
  
  return 0;
}
