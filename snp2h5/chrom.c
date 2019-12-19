
#include <limits.h>
#include <string.h>

#include "memutil.h"
#include "chrom.h"
#include "util.h"



/**
 * Returns appropriate chromosome for filename by looking for match
 * with chromosome name in filename. Returns NULL if no match found.
 */
Chromosome *chrom_guess_from_file(char *path,
				  Chromosome *chroms,
				  int n_chrom) {
  int i, j, longest_match, n1, n2;
  char *filename;
  Chromosome *match_chrom;

  longest_match = 0;
  match_chrom = NULL;

  
  /* find last occurance of '/' or '\' in file
   * assume that everything before this is in directories
   * leading up to filename
   */
  filename = rindex(path, '/');
  if(filename == NULL) {
    filename = rindex(path, '\\');
    if(filename == NULL) {
      filename = path;
    }
  }

  fprintf(stderr, "guessing chromosome name from filename %s\n",
	  filename);
  
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

  if(match_chrom) {
    fprintf(stderr, "best matching chromosome: %s\n", match_chrom->name);
  } else {
    fprintf(stderr, "could not find matching chromosome\n");
  }
  
  return match_chrom;
}



/**
 * Returns a (deep) copy of the provided chromosome
 */
Chromosome *chrom_copy(const Chromosome *chrom) {
  Chromosome *new_chrom;

  new_chrom = my_new(Chromosome, 1);
  new_chrom->id = chrom->id;

  if(chrom->name) {
    new_chrom->name = util_str_dup(chrom->name);
  } else {
    new_chrom->name = NULL;
  }

  if(chrom->assembly) {
    new_chrom->assembly = util_str_dup(chrom->assembly);
  } else {
    new_chrom->assembly = NULL;
  }

  new_chrom->len = chrom->len;
  
  return new_chrom;
}



Chromosome *chrom_read_gzfile(const char *filename, int *n_chrom) {
  char buf[LINE_MAX], name_buf[LINE_MAX];
  Chromosome *chroms;
  gzFile gzf;
  int n, i;

  gzf = util_must_gzopen(filename, "r");
  *n_chrom = util_gzcount_lines(gzf);
  
  if(*n_chrom < 1) {
    my_err("%s:%d: chromosome file '%s' is empty\n", 
	   __FILE__, __LINE__, filename);
  }

  chroms = my_new(Chromosome, *n_chrom);
  for(i = 0; i < *n_chrom; i++) {
    if(!gzgets(gzf, buf, sizeof(buf))) {
      my_err("%s:%d: expected %d lines in file, but only read %d\n", 
	     __FILE__, __LINE__, *n_chrom, i);
    }

    n = sscanf(buf, "%s %ld", name_buf, &chroms[i].len);
    if(n < 2) {
      my_err("%s:%d: line did not have at least 2 tokens:\n'%s'",
	     __FILE__, __LINE__, buf);
    }
    chroms[i].name = util_str_dup(name_buf);
    chroms[i].assembly = NULL;
    chroms[i].id = i;
    
    if(chroms[i].len < 1) {
      my_err("%s:%d: chrom length (%ld) should be >= 1",
	     __FILE__, __LINE__, chroms[i].len);
    }
  }

  gzclose(gzf);

  return chroms;
}



/**
 * Reads an array chromosomes from a file containing a name and length on
 * each line, separated by a white space character.
 */
Chromosome *chrom_read_file(const char *filename, int *n_chrom) {
  char buf[LINE_MAX], name_buf[LINE_MAX];
  Chromosome *chroms;
  FILE *f;
  int n, i;

  if(util_has_gz_ext(filename)) {
    return chrom_read_gzfile(filename, n_chrom);
  }

  f = util_must_fopen(filename, "r");
  *n_chrom = util_fcount_lines(f);
  
  if(*n_chrom < 1) {
    my_err("%s:%d: chromosome file '%s' is empty\n", 
	   __FILE__, __LINE__, filename);
  }

  chroms = my_new(Chromosome, *n_chrom);
  for(i = 0; i < *n_chrom; i++) {
    if(!fgets(buf, sizeof(buf), f)) {
      my_err("%s:%d: expected %d lines in file, but only read %d\n", 
	     __FILE__, __LINE__, *n_chrom, i);
    }

    n = sscanf(buf, "%s %ld", name_buf, &chroms[i].len);
    if(n < 2) {
      my_err("%s:%d: line did not have at least 2 tokens:\n'%s'",
	     __FILE__, __LINE__, buf);
    }
    chroms[i].name = util_str_dup(name_buf);
    chroms[i].assembly = NULL;
    chroms[i].id = i;
    
    if(chroms[i].len < 1) {
      my_err("%s:%d: chrom length (%ld) should be >= 1",
	     __FILE__, __LINE__, chroms[i].len);
    }
  }

  fclose(f);

  return chroms;
}



/**
 * Frees memory allocated for an array of chromosomes
 */
void chrom_array_free(Chromosome *chroms, int n_chrom) {
  int i;

  for(i = 0; i < n_chrom; i++) {
    if(chroms[i].name) {
      my_free(chroms[i].name);
    }
    if(chroms[i].assembly) {
      my_free(chroms[i].assembly);
    }
  }

  my_free(chroms);
}


/**
 * Frees memory that was allocated for a single chromosome
 */
void chrom_free(Chromosome *chrom) {
  if(chrom->name != NULL) {
    my_free(chrom->name);
  }

  if(chrom->assembly != NULL) {
    my_free(chrom->assembly);
  }

  my_free(chrom);
}

