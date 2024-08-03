#ifndef COUNTWORDS
#define COUNTWORDS

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#define FILE_ERROR 1
#define MEMORY_ERROR 2

extern int error_type;
extern char *error_message;

long **count_short_words(char *, int);
long **count_bipart_words(char *, int, int, int);

#endif
