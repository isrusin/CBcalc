#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>

#include "count_words.h"

extern char* error_message;
extern int error_type;
extern const int FILE_ERROR;
extern const int MEMMORY_ERROR;

void site_to_string(unsigned long, int, char *);
void print_counts(FILE *, int, int, int, long **);

int main(int argc, char **argv){
    if(argc != 6){
        printf("usage: ./count_short_words input.fasta[.gz] output.cnt "
               "word_lengtn gap_position gap_length\n");
        return 1;
    }
    const int len = atoi(argv[3]);
    if(len < 1 || len > 14){
        printf("bad word length, only [1..14] is allowed\n");
        return 1;
    }
    const int pos = atoi(argv[4]);
    if(pos < 1 || pos >= len){
        printf("bad gap position\n");
        return 1;
    }
    const int gap = atoi(argv[5]);
    if(len < 1 || len > 14){
        printf("bad gap length, only [1..14] is allowed\n");
        return 1;
    }
    long **counts = count_bipart_words(argv[1], len, pos, gap);
    if(!counts){
        if(error_type == FILE_ERROR)
            printf("Input file error: %s\n", error_message);
        else if(error_type == MEMMORY_ERROR)
            printf("Memmory error: %s\n", error_message);
        else
            printf("Unknown error\n");
        return 1;
    }
    FILE *cnt;
    cnt = fopen(argv[2], "wb");
    assert(cnt);
    print_counts(cnt, len, pos, gap, counts);
    fclose(cnt);
    int i;
    for(i = 0; i < len - pos; i ++){
        free(counts[i]);
    }
    free(counts);
}

void site_to_string(unsigned long site, int len, char *dest){
    const char nucls[] = {'A', 'C', 'G', 'T'};
    int i = len - 1;
    while(i >= 0){
        dest[i] = nucls[site % 4];
        site >>= 2;
        i --;
    }
}

void print_counts(FILE *cnt, int len, int pos, int gap, long **countsp){
    unsigned long max = 1ul << 2 * (pos + 1);
    char dest[len + 1];
    long *counts;
    int i;
    for(i = pos + 1; i <= len; i ++){
        dest[i] = '\0';
        counts = countsp[len - i];
        fprintf(cnt, "#word_length=%d, gap_position=%d, gap_length=%d, "
                "total_number=%d\n", i, pos, gap, counts[max]);
        unsigned long site;
        for(site = 0ul; site < max; site ++){
            site_to_string(site, i, dest);
            fprintf(cnt, "%s\t%d\n", dest, counts[site]);
        }
        max <<= 2;
    }
}
