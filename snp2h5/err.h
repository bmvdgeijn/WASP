#ifndef __ERR_H__
#define __ERR_H__

#include <stdio.h>
#include <stdarg.h>

void my_err(const char *format, ...);
void my_warn(const char *format, ...);
void my_verbose(const char *format, ...);

#endif
