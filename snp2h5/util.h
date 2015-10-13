#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#include <zlib.h>


#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif 

#define UTIL_FGETS_BUF_SZ 16348
#define UTIL_STR_BUF_SZ 20


#define util_fread_one(file, var) util_must_fread((file), &(var), sizeof(var))
#define util_fwrite_one(file, var) util_must_fwrite((file), &(var), sizeof(var))

#define util_gzread_one(file, var) util_must_gzread((file), &(var), sizeof(var))
#define util_gzwrite_one(file, var) util_must_gzwrite((file), &(var), sizeof(var))


int util_has_gz_ext(const char *filename);

char *util_read_entire_file(char *filename);
long util_count_lines(const char *filename);
long util_fcount_lines(FILE *fh);
long util_gzcount_lines(gzFile gzf);
long util_fcount_lines_match(FILE *fh, const char *starts_with);
long util_gzcount_lines_match(gzFile gzf, const char *starts_with);

char *util_fgets_line(FILE *fh);
char *util_gzgets_line(gzFile gzf);
size_t util_gzgetline(gzFile gzf, char **lineptr, size_t *size);

FILE *util_must_fopen(const char *path, const char *mode);
gzFile util_must_gzopen(const char *path, const char *mode);
gzFile util_check_gzopen(const char *path);
void util_must_fread(FILE *file, void *buf, size_t size);
void util_must_fwrite(FILE *file, void *buf, size_t size);
void util_must_gzwrite(gzFile gzf, void *buf, size_t size);
void util_must_gzread(gzFile gzf, void *buf, size_t size);
int util_file_exists(const char *path);

char *util_long_to_comma_str(const long x);

void util_reverse(void *base, const size_t nmemb, const size_t size);
void util_breverse(void *base, const size_t nmemb);

char *util_str_concat(const char *str1, ...);
char *util_str_dup(const char *str);
char *util_str_ndup(const char *str, size_t n);
size_t util_strncpy(char *dest, const char *src, size_t n);
void util_str_replace(char *str, char from, char to);
void util_str_reverse(char *str);
void util_str_uc(char *str);
void util_str_lc(char *str);
void util_str_remove_char(char *str, const char c);
void util_str_remove_whitespace(char *str);
void util_str_lstrip(char *str);
void util_str_rstrip(char *str);
void util_str_strip(char *str);
int util_str_starts_with(const char *str, const char *start);
int util_str_ends_with(const char *str, const char *end);
int util_str_split(char *str, char **tokens, const size_t n_tok);

int util_dbl_cmp(const void *x, const void *y);

long util_parse_long(const char *str);
double util_parse_double(const char *str);


#endif
