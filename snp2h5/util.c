
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <zlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>

#include "memutil.h"
#include "util.h"
#include "err.h"


/**
 * Returns TRUE if the provided filename ends with ".gz", returns
 * FALSE otherwise.
 */
int util_has_gz_ext(const char *filename) {
  size_t len;

  len = strlen(filename);

  if((len > 3) && (strncmp(&filename[len-3], ".gz",3)==0)) {
    return TRUE;
  }
  return FALSE;
}



/**
 * Helper function, reads an entire file into memory and returns the
 * result as a null-terminated string. The returned string should be
 * freed when it is no longer needed.
 */
char *util_read_entire_file(char *filename) {
  FILE *fh;
  struct stat fs;
  char *buf;

  fh = fopen(filename, "r");

  if(fh == NULL) {
    my_err("%s:%d: Could not open tree file '%s'", __FILE__, 
	    __LINE__, filename);
  }

  /* determine file size */
  fstat(fileno(fh), &fs);

  /* read entire file at once */
  buf = my_new(char, fs.st_size + 1);
  if(fread(buf, fs.st_size, 1, fh) == 0) {
    my_err("%s:%d: Could not read entire file '%s'", __FILE__,
	    __LINE__, filename);
  }

  buf[fs.st_size] = '\0';

  return buf;
}



/**
 * Counts the number of newline characters in the file pointed to by
 * filehandle. The filehandle is rewound to the beginning of the file
 * before and after the count of newlines.
 */
long util_fcount_lines(FILE *fh) {
  long line_count;
  char c;

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  line_count = 0;
  while((c = getc(fh)) != EOF) {
    if(c == '\n') {
      line_count++;
    }
  }

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}





/**
 * write a subset of lines which are flagged as TRUE in the provided
 * array (that should have a length corresponding to the number lines
 * in the file)
 */
void util_fwrite_line_subset(FILE *input_fh, FILE *output_fh,
			     long n_lines, unsigned char *line_flags) {
  long i, len;
  char buf[UTIL_FGETS_BUF_SZ];
  
  i = 0;
  while(fgets(buf, UTIL_FGETS_BUF_SZ, input_fh) != NULL) {
    len = strlen(buf);
    if(buf[len-1] != '\n') {
      my_err("%s:%d: line %ld is too long\n", __FILE__, __LINE__, i+1);
    }

    if(i >= n_lines) {
      my_err("%s:%d: line has more lines than expected (%ld)", 
	     __FILE__, __LINE__, n_lines);
    }
    
    if(line_flags[i]) {
      /* write line to output file */
      fprintf(output_fh, "%s", buf);
    }
	
    i += 1;
  }
}




/**
 * Counts number of newline characters in file with provided path
 */
long util_count_lines(const char *filename) {
  char c;
  long line_count;
  gzFile gzf;

  gzf = util_must_gzopen(filename, "rb");

  line_count = 0;
  while((c = gzgetc(gzf)) != EOF) {
    if(c == '\n') {
      line_count += 1;
    }
  }
  
  gzclose(gzf);

  return line_count;
}


/**
 * Counts number of newline characters in file, but stops counting
 * if max_lines is reached.
 */
long util_count_lines_max(const char *filename, long max_lines) {
  char c;
  long line_count;
  gzFile gzf;

  gzf = util_must_gzopen(filename, "rb");

  line_count = 0;
  while((c = gzgetc(gzf)) != EOF) {
    if(c == '\n') {
      line_count += 1;

      if(line_count >= max_lines) {
	return max_lines;
      }
    }
  }
  
  gzclose(gzf);

  return line_count;
}




/**
 * Counts the number of '\n' characters in the provided gzipped file.
 * The gzFile is rewound to the beginning of the file before and after the count of
 * newlines.
 */
long util_gzcount_lines(gzFile gzf) {
  char c;
  long line_count;

  if(gzrewind(gzf) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  line_count = 0;

  while((c = gzgetc(gzf)) != EOF) {
    if(c == '\n') {
      line_count += 1;
    }
  }


  if(gzrewind(gzf) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}




/**
 * Counts the number of lines in a file that begin with a provided
 * string. The filehandle is rewound to the beginning of the file
 * before and after the count of newlines.
 */
long util_fcount_lines_match(FILE *fh, const char *starts_with) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  int started_with, at_line_start;
  int match_len, len;

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  match_len = strlen(starts_with);
  if(match_len == 0) {
    return 0;
  }

  if(match_len > UTIL_FGETS_BUF_SZ) {
    my_err("%s:%d: length of string to match must be"
	    "<= %d bytes", __FILE__, __LINE__, UTIL_FGETS_BUF_SZ);
  }

  line_count = 0;
  started_with  = FALSE;
  at_line_start = TRUE;

  while(fgets(buf, UTIL_FGETS_BUF_SZ, fh) != NULL) {
    if(at_line_start) {
      /* we are at the beginning of a line, does it match the string? */
      if(strncmp(starts_with, buf, match_len) == 0) {
	started_with = TRUE;
      }
    }

    len = strlen(buf);
    if(buf[len-1] == '\n') {
      /* we are at the end of the line, did the beginning match? */
      if(started_with) {
	line_count++;
	started_with = FALSE;
      }
      at_line_start = TRUE;
    } else {
      at_line_start = FALSE;
    }
  }

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}



/**
 * Counts the number of lines in a gzipped file that begin with a
 * provided string. The gzFile is rewound to the beginning of the
 * file before and after the count of newlines.
 */
long util_gzcount_lines_match(gzFile gzf, const char *starts_with) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  int started_with, at_line_start;
  int match_len, len;

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  match_len = strlen(starts_with);
  if(match_len == 0) {
    return 0;
  }

  if(match_len > UTIL_FGETS_BUF_SZ) {
    my_err("%s:%d: length of string to match must be"
	    "<= %d bytes", __FILE__, __LINE__, UTIL_FGETS_BUF_SZ);
  }

  line_count = 0;
  started_with  = FALSE;
  at_line_start = TRUE;

  while(gzgets(gzf, buf, UTIL_FGETS_BUF_SZ) != NULL) {
    if(at_line_start) {
      /* we are at the beginning of a line, does it match the string? */
      if(strncmp(starts_with, buf, match_len) == 0) {
	started_with = TRUE;
      }
    }

    len = strlen(buf);
    if(buf[len-1] == '\n') {
      /* we are at the end of the line, did the beginning match? */
      if(started_with) {
	line_count++;
	started_with = FALSE;
      }
      at_line_start = TRUE;
    } else {
      at_line_start = FALSE;
    }
  }

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}



/**
 * returns TRUE if file with provided path exists, FALSE otherwise
 */
int util_file_exists(const char *path) {
  struct stat s;
  return stat(path, &s) == 0;
}


/**
 * concatenates a variable number of strings together into a single
 * string and returns the result. The last argument provided must be
 * NULL, to indicate the end of the variable arguments.  The returned
 * string should be freed when it is no longer needed.
 *
 * example:
 *   
 *
 *  char *my_str = util_str_concat("hello", " ", "world!", NULL);
 *  fprintf(stderr, "str: %s\n", my_str);
 *  my_free(my_str);
 *
 */
char *util_str_concat(const char *str1, ...) {
  size_t ttl_len;
  char *next_str, *new_str;
  va_list args;
  
  if(str1 == NULL) {
    return NULL;
  }

  /* first pass through variable args: determine total size of str */
  va_start(args, str1);

  ttl_len = strlen(str1);

  while(1) {
    next_str = va_arg(args, char *);
    if(next_str == NULL) {
      /* end of args list */
      break;
    }
    ttl_len += strlen(next_str);
  }

  va_end(args);


  /* second pass through args: paste strings together */
  new_str = my_new(char, ttl_len + 1);
  
  strcpy(new_str, str1);
  ttl_len = strlen(str1);

  va_start(args, str1);
  while(1) {
    next_str = va_arg(args, char *);
    if(next_str == NULL) {
      break;
    }

    /* append next string to end of current str */
    strcpy(&new_str[ttl_len], next_str);

    ttl_len += strlen(next_str);
  }
  va_end(args);

  return new_str;
}




/**
 * Reads a single line from provided the gzfile handle and returns it 
 * as null-terminated string (the '\n' is removed). The memory for the line
 * is dynamically allocated and should be freed.
 * Returns NULL if at EOF.
 */
char *util_gzgets_line(gzFile gzf) {
  char *line, c;
  size_t size, i;

  size = 65536;
  line = my_malloc(size);

  i = 0;
  while((c = gzgetc(gzf)) != -1) {
    if(i >= size) {
      /* buffer is full, double its size */
      size = size*2;
      /* fprintf(stderr, "expanding mem to %ld bytes\n", size); */
      line = my_realloc(line, size);
    }
    
    if(c == '\n') {
      /* end of line, terminate  */
      line[i] = '\0';
      return line;
    }

    line[i] = c;
    i++;
  }

  if(i == 0) {
    /* at EOF */
    return NULL;
  }

  
  /* did not hit a '\n' before EOF, add terminating '\0' to line 
   * before returning it
   */
  if(i+1 >= size) {
    line = my_realloc(line, i+1);
  }
  line[i+1] = '\0';
  
  return line;
}


/**
 * Reads a single line from provided the gzfile handle into the provided
 * buffer. The buffer memory is reallocated as required and the buffer size
 * is set appropriately. Newlines chars are replaced with terminating '\0'. 
 * Returns the length of the string read from file (typically one less than 
 * number of bytes read because newline replaced with '\0') or -1 if at EOF.
 */
size_t util_gzgetline(gzFile gzf, char **lineptr, size_t *size) {
  char c;
  size_t i;

  i = 0;
  while((c = gzgetc(gzf)) != -1) {
    if(i >= *size) {
      /* buffer is full, double its size */
      *size = *size * 2;
      /* fprintf(stderr, "expanding mem to %ld bytes\n", *size); */
      *lineptr = my_realloc(*lineptr, *size);
    }
    
    if(c == '\n') {
      /* end of line, terminate  */
      (*lineptr)[i] = '\0';
      return i;
    }

    (*lineptr)[i] = c;
    i++;
  }

  if(i == 0) {
    /* at EOF */
    return -1;
  }

  
  /* did not hit a '\n' before EOF, add terminating '\0' to line 
   * before returning it
   */
  if(i+1 >= *size) {
    *lineptr = my_realloc(*lineptr, i+1);
    *size = i+1;
  }
  (*lineptr)[i+1] = '\0';
  
  return i;
}

				    

/**
 * Like fgets, but reads an entire line, without the
 * need to provide a buffer of a fixed size.
 * The returned line should be freed when it is no
 * longer needed.
 * The newline character is replaced by a '\0'
 */
char *util_fgets_line(FILE *fh) {
  char *line, c;
  size_t size, i;

  size = 1024;
  line = my_malloc(size);

  i = 0;
  while((c = fgetc(fh)) != -1) {
    if(i >= size) {
      /* buffer is full, double its size */
      size = size*2;
      line = my_realloc(line, size);
    }
    
    if(c == '\n') {
      /* end of line, terminate  */
      line[i] = '\0';
      return line;
    }

    line[i] = c;
    i++;
  }

  if(i == 0) {
    /* at EOF */
    return NULL;
  }

  
  /* did not hit a '\n' before EOF, add terminating '\0' to line 
   * before returning it
   */
  if(i+1 >= size) {
    line = my_realloc(line, i+1);
  }
  line[i+1] = '\0';
  
  return line;
}
  



/**
 * Makes a copy of a string and returns it. The new string should
 * be freed when it is no longer needed
 */
char *util_str_dup(const char *str) {
  char *new_str;
  new_str = my_new(char, strlen(str)+1);
  strcpy(new_str, str);

  return new_str;
}


/**
 * Makes a copy of a string up to n characters long and returns it.
 * If the string is longer than n characters, then a string of n
 * characters (plus a '\0') is returned. The new string should be
 * freed when it is no longer needed
 */
char *util_str_ndup(const char *str, size_t n) {
  char *new_str;
  size_t len;

  len = 0;
  while(len < n && str[len] != '\0') {
    len++;
  }

  new_str = my_new(char, len+1);
  strncpy(new_str, str, len);
  new_str[len] = '\0';

  return new_str;
}


/** 
 * more practical version of strncpy that:
 * [1] returns length of dest string
 * [2] always null terminates by setting dest[n-1] to '\0'
 * 
 * This function still gives no indication that src string is
 * truncated, however.
 */

size_t util_strncpy(char *dest, const char *src, size_t n) {
  size_t i = 0;
  if(n > 0) {
    while((i < n-1) && (src[i] != '\0')) {
      dest[i] = src[i];

      i++;
    }
    /* padd to end with \0 */
    memset(&dest[i], '\0', n-i);
    /* dest[i] = '\0'; */
  }
	 
  return i;
}


/**
 * Does an in-place uppercase of all characters in a string
 */
void util_str_uc(char *str) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    str[i] = toupper(str[i]);
    i++;
  }

  return;
}



/**
 * Does an in-place lowercase of all characters in a string
 */
void util_str_lc(char *str) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    str[i] = tolower(str[i]);
    i++;
  }

  return;
}


/*
 * Does an in-place replacement of one character by another
 * in a string.
 */
void util_str_replace(char *str, char from, char to) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    if(str[i] == from) {
      str[i] = to;
    }
    i++;
  }
}


/*
 * Does an in-place removal of all instances of
 * a specified character in a string.
 */
void util_str_remove_char(char *str, const char c) {
  size_t i,j;

  i = j = 0;
  while(str[i] != '\0') {
    if(str[i] != c) {
      str[j] = str[i];
      j++;
    }
    i++;
  }
  str[j] = '\0';
}






/**
 * Converts a long integer to a string with commas delimited
 * the hundreds thousands, millions etc. The returned string should be
 * freed when it is no longer needed.
 */
char *util_long_to_comma_str(const long x) {
  char *str, *comma_str;
  int len,i,j,k;
  int max_commas;

  str = my_new(char, UTIL_STR_BUF_SZ);
  max_commas = (UTIL_STR_BUF_SZ / 3);
  
  len = snprintf(str, UTIL_STR_BUF_SZ-max_commas, "%lu", x);

  if(len == 0) {
    my_err("%s:%d: could not write number %lu to string", 
	    __FILE__, __LINE__, x);
  }
  
  max_commas = len / 3;
  comma_str = my_new(char, len + max_commas + 1);

  /* copy the string in reverse, adding commas every three */
  j = 0;
  k = 0;
  for(i = len-1; i >= 0; i--) {
    comma_str[j] = str[i];
    k++;
    
    if(i > 0 && ((k %3) == 0)) {
      /* insert comma */
      j++;
      comma_str[j] = ',';
    }
    j++;
  }

  comma_str[j] = '\0';
  my_free(str);

  util_breverse(comma_str, j);

  return comma_str;
}



/**
 * Reverses the contents of a provided array. Arguments
 * are similar to those for qsort():
 *   base - ptr to the beginning of the array
 *   nmemb - number of elements in array
 *   size - size of each element in bytes
 */
void util_reverse(void *base, const size_t nmemb, const size_t size) {
  void *tmp;
  unsigned char *fwd_ptr, *rev_ptr;

  /* allocate mem for holding swap */
  tmp = malloc(size);
  if(tmp == NULL) {
    my_err("%s:%d: could not allocate mem", __FILE__, __LINE__);
  }
  
  /* initiailize pointers so that they point to beginning and end of array */
  fwd_ptr = (unsigned char *)base;
  rev_ptr = &((unsigned char *)base)[nmemb*size - size];

  /* move ptrs towards each other from either end of array */
  while(fwd_ptr < rev_ptr) {
    /* swap memory at each pointer */
    tmp = memcpy(tmp, fwd_ptr, size);
    memcpy(fwd_ptr, rev_ptr, size);
    memcpy(rev_ptr, tmp, size);

    /* advance ptrs */
    fwd_ptr += size;
    rev_ptr -= size;
  }

  free(tmp);
}


/**
 * Reverses a byte array of a given length
 */
void util_breverse(void *base, const size_t nmemb) {
  util_reverse(base, nmemb, sizeof(unsigned char));
}

/*
 * Does an in-place reversal of the order of the characters within a
 * string. This function calls strlen to determine the length of the
 * string. If this is already known, or this function will be called
 * repeatedly it is more efficient to call util_breverse instead.
 *
 */
void util_str_reverse(char *str) {
  util_breverse(str, strlen(str));
}


/*
 * Performs an in-place removal of whitespace characters from the
 * provided string.
 */
void util_str_remove_whitespace(char *str) {
  long i,j;

  i = j = 0;
  while(str[i] != '\0') {
    if(!isspace((unsigned char)str[i])) {
      str[j] = str[i];
      j++;
    }
    i++;
  }

  str[j] = '\0';
}


/*
 * Performs in-place removal of leading whitespace
 * characters from the provided string.
 */
void util_str_lstrip(char *str) {
  long leading_ws, len;
  
  leading_ws = 0;
  while(str[leading_ws] != '\0' && isspace((unsigned char)str[leading_ws])) {
    leading_ws++;
  }

  if(leading_ws > 0) {
    len = strlen(str);
    memmove(str, &str[leading_ws], len-leading_ws+1);
  }
}


/*
 * Performs in-place removal of trailing whitespace
 * characters from the provided string.
 */
void util_str_rstrip(char *str) {
  long len,i;

  len = strlen(str);
  
  i = len-1;
  while(i >= 0 && isspace((unsigned char)str[i])) {
    str[i] = '\0';
    i--;
  }
}


/**
 * Performs in-place removal of leading and trailing whitespace
 * from provided string.
 */
void util_str_strip(char *str) {
  util_str_lstrip(str);
  util_str_rstrip(str);
}


/**
 * Returns TRUE if the start of str exactly matches the characters in
 * start
 */
int util_str_starts_with(const char *str, const char *start) {
  size_t i;

  i = 0;

  while(start[i] != '\0') {
    if(str[i] != start[i]) {
      return FALSE;
    }
    i++;
  }

  return TRUE;
}





/**
 * Returns TRUE if the end of str exactly matches the characters in
 * end, returns FALSE otherwise.
 */
int util_str_ends_with(const char *str, const char *end) {
  size_t i, j;

  i = strlen(str);
  j = strlen(end);

  if(j > i) {
    /* string to match is longer than provided string */
    return FALSE;
  }

  while(j > 0) {
    if(str[i-1] != end[j-1]) {
      /* character does not match */
      return FALSE;
    }
    i--;
    j--;
  }

  return TRUE;
}


/**
 * Splits a string into n_tok or fewer tokens
 * using space and tab as delimitors.  Returns the
 * number of tokens read. The provided string is modified by replacing
 * delimitors with '\0'.
 */
int util_str_split(char *str, char **tokens, const size_t n_tok) {  
  size_t i = 0;

  if(n_tok < 1) {
    my_err("%s:%d: number of tokens must be at least 1\n",
	   __FILE__, __LINE__);
  }

  while(((tokens[i] = strsep(&str, " \t")) != NULL)) {
    if(tokens[i][0] != '\0') {
      i++;

      if(i == n_tok) {
	break;
      }
    }
  }
  
  return i;
}


/**
 * Comparison function used to compare two doubles Useful for sorting.
 * Returns -1 if x < y, returns 0 if x == y, returns 1 if x > y
 */
int util_dbl_cmp(const void *x, const void *y) {
  if(*(double *)x < *(double *)y) {
    return -1;
  }
  if(*(double *)x > *(double *)y) {
    return 1;
  }
  return 0;
}


/**
 * Reads size bytes from a gzipped file or prints an error and aborts.
 */
void util_must_gzread(gzFile gzf, void *buf, size_t size) {
  int errnum;

  if(gzread(gzf, buf, size) != size) {
    if(gzerror(gzf, &errnum)) {
      my_err("error reading %lld bytes: %s", (long long)size,
	     gzerror(gzf, &errnum));
    } else {
      my_err("end of file reading %lld bytes", (long long)size);
    }
  }
}

/**
 * Reads size bytes from a file or prints an error and aborts. 
 * Based on Jim Kent's mustRead function
 */
void util_must_fread(FILE *file, void *buf, size_t size) {
  if(size != 0 && fread(buf, size, 1, file) != 1) {
    if(ferror(file)) {
      my_err("error reading %lld bytes: %s", (long long)size, 
	     strerror(ferror(file)));
    } else {
      my_err("end of file reading %lld bytes", (long long)size);
    }
  }
}




/**
 * Writes size bytes to a file or prints an error and aborts.
 */
void util_must_fwrite(FILE *file, void *buf, size_t size) {
  if(size != 0 && fwrite(buf, size, 1, file) != 1) {
    if(ferror(file)) {
      my_err("error writing %lld bytes: %s", (long long)size, 
            strerror(ferror(file)));
    } else {
      my_err("error writing %lld bytes\n", (long long)size);
    }
  }
}


/**
 * Writes size bytes to a gzipped file or prints an error and aborts.
 */
void util_must_gzwrite(gzFile gzf, void *buf, size_t size) {
  int errnum;

  if(gzwrite(gzf, buf, size) != size) {
    if(gzerror(gzf, &errnum)) {
      my_err("error writing %lld bytes: %s", (long long)size,
	     gzerror(gzf, &errnum));
    } else {
      my_err("error writing %lld bytes\n", (long long)size);
    }
  }
}





/** 
 * Wrapper around fopen that prints an error and aborts on
 * failure to open the file
 */
FILE *util_must_fopen(const char *path, const char *mode) {
  FILE *f;

  errno = 0;
  f = fopen(path, mode);
  if(f == NULL) {
    my_err("%s:%d: error opening file '%s' in mode '%s': %s", 
	   __FILE__, __LINE__, path, mode, strerror(errno));
  }

  return f;
}



/** 
 * Wrapper around fopen that prints an error and aborts on
 * failure to open the file
 */
gzFile util_must_gzopen(const char *path, const char *mode) {
  gzFile f;

  if(strcmp(ZLIB_VERSION, zlibVersion()) != 0) {
    my_warn("zlib compilation (%s) and runtime (%s) versions do not match.\n",
	    ZLIB_VERSION, zlibVersion());
  }
  
  if((strcmp(mode, "wb") == 0) && !util_str_ends_with(path, ".gz")) {
    my_warn("%s:%d: file '%s' does not end with '.gz' but "
	    "opening with gzopen anyway\n", __FILE__, __LINE__,
	    path);
  }
  
  errno = 0;
  
  f = gzopen(path, mode);
  if(f == NULL) {
    if(errno) {
      my_err("%s:%d: error opening file '%s' in mode '%s': %s", 
	     __FILE__, __LINE__, path, mode, strerror(errno));
    } else {
      my_err("%s:%d: error opening file '%s' in mode '%s'", 
	     __FILE__, __LINE__, path, mode);
    }
  }

  return f;
}


/**
 * opens a gzipped file in write mode, but aborts if the file already
 * exists.
 */
gzFile util_check_gzopen(const char *path) {
  if(util_file_exists(path)) {
    my_err("%s:%d: output file already exists: %s",
	   __FILE__, __LINE__, path);
  }
  return util_must_gzopen(path, "wb");
}



/**
 * Wrapper around strtol that aborts and writes an error message on
 * failure
 */
long util_parse_long(const char *str) {
  long val;

  errno = 0;
  val = strtol(str, NULL, 10);
  if(val == 0 && errno) {
    perror("strtol");
    my_err("%s:%d: failed to parse long string '%s'", __FILE__, __LINE__, str);
  }

  return val;
}


/**
 * returns true if the provided string is one of 'na' or 'nan', ignoring case
 */
static int is_nan_str(const char *str) {
  int i;
  char uc_str[4];

  /* convert first three chars to uppercase string */
  i = 0;
  while((str[i] != '\0') && (!isspace(str[i])) && (i < 3)) {
    uc_str[i] = toupper(str[i]);
    i++;
  }
  uc_str[i] = '\0';

  /* is this a NAN? */
  if(i == 3) {
    return (strncmp(uc_str, "NAN", 3) == 0);
  }

  /* is this a NA? */
  return (strncmp(uc_str, "NA", 2) == 0);
}

/**
 * Wrapper around strtod that aborts and writes an error message on
 * failure
 */
double util_parse_double(const char *str) {
  double val;

  if(is_nan_str(str)) {
    return NAN;
  }

  errno = 0;
  val = strtod(str, NULL);
  if(errno) {
    perror("strtod");
    my_err("%s:%d: failed to parse double string '%s'", __FILE__, __LINE__,
	   str);
  }

  return val;
}
