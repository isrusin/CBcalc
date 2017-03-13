#ifndef COUNTWORDS
#define COUNTWORDS

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

int error_type;
const int FILE_ERROR;
const int MEMORY_ERROR;
char *error_message;

long **count_short_words(char *, int);
long **count_bipart_words(char *, int, int, int);

#endif
