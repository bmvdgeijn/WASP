
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "err.h"
#include "memutil.h"

/**
 * Wrapper for malloc that prints an error message and exits
 * if memory could not be allocated
 */
void *my_malloc(size_t n_bytes) {
  void *mem;
 
  mem = malloc(n_bytes);
  if(!mem) {
    my_err("%s:%d: could not allocate %ld bytes of memory\n",
	  __FILE__, __LINE__, n_bytes);
  }

  return mem;
}


/**
 * Wrapper for malloc that allocates memory and 
 * then sets every byte to 0 (NULL).
 */
void *my_malloc0(size_t n_bytes) {
  void *mem;

  mem = my_malloc(n_bytes);
  memset(mem, 0, n_bytes);
  return mem;
}


/**
 * Wrapper for realloc that prints an error message and exits
 * if memory could not be allocated
 */
void *my_realloc(void *ptr, size_t n_bytes) {
  ptr = realloc(ptr, n_bytes);
  if(!ptr) {
    my_err("%s:%d: could not allocate %ld bytes of memory\n",
	  __FILE__, __LINE__, n_bytes);
  }

  return ptr;
}



/**
 * Wrapper around free that checks that provided
 * ptr is not NULL before calling free
 */
void __MY_FREE(void *ptr, const char *filename, const int line_num) {
  if(ptr == NULL) {
    my_err("%s:%d: attempt to free NULL ptr at: %s:%d\n",
	  __FILE__, __LINE__, filename, line_num);
  }
  free(ptr);
}



