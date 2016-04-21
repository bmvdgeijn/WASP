#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "err.h"

void my_err(const char *format, ...) {
  int i;

  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, format, args);
  
  i = strlen(format);
  if(i < 1 || format[i-1] != '\n') {
    /* add new line */
    fprintf(stderr, "\n");
  }

  va_end(args);
  exit(2);
}




void my_warn(const char *format, ...) {
  int i;

  va_list args;
  va_start(args, format);
  fprintf(stderr, "WARNING: ");
  vfprintf(stderr, format, args);
  
  i = strlen(format);
  if(i < 1 || format[i-1] != '\n') {
    /* add new line */
    fprintf(stderr, "\n");
  }

  va_end(args);
}


void my_verbose(const char *format, ...) {
  #ifdef VERBOSE
  va_list args;
  va_start(args, format);
  vafprintf(stderr, format, args)
  va_end(args);
  #endif
}
