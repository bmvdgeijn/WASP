
#ifndef __MEMUTIL_H__

#define __MEMUTIL_H__

#include <stdlib.h>
#include "err.h"


void *my_malloc(size_t n_bytes);
void *my_malloc0(size_t n_bytes);
void *my_realloc(void *ptr, size_t n_bytes);
void __MY_FREE(void *ptr, const char *filename, const int line_num);

#define my_new(struct_type, n) \
  ((struct_type *)my_malloc(sizeof(struct_type) * n))

#define my_new0(struct_type, n) \
  ((struct_type *)my_malloc0(sizeof(struct_type) * n))

#define my_free(ptr) __MY_FREE(ptr, __FILE__, __LINE__); (ptr) = (void *)NULL

#endif
