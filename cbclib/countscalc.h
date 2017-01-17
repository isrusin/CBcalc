#ifndef COUNTWORDS
#define COUNTWORDS

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

int error_type;
const int FILE_ERROR;
const int MEMMORY_ERROR;
char *error_message;


typedef struct {
    int len;
    unsigned long mask;
    unsigned long val;
    int index;
} site_t;

typedef struct {
    /*  Site     len  mask    shift  dmask   ushift
        x(N)xxx  3    111111  6      111111  4
        xxx(N)x  3    111111  2      11      0
        xx(N)xx  2    1111    4      1111    0
    */
    int size;
    int len;
    unsigned int mask, dmask;
    int shift, ushift;
    unsigned int *arr;
    int index;
} bipart_t;

int translate(char);
int skip(gzFile);
int initialize_short(gzFile, site_t *);
int initialize_bipart(gzFile, bipart_t *);
void countup_short(gzFile, site_t *, long *);
void countup_bipart(gzFile, bipart_t *, long *);
void countup_short_res(site_t *, long **);
void countup_bipart_res(bipart_t *, int, long **);
void countup_subsites(int, int, long **);
bipart_t make_bpsite(int, int, int, int, int);
long **count_short_words(char *, int);
long **count_bipart_words(char *, int, int, int);

#endif
