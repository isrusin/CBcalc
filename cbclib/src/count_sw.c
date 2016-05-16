#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>

#include "count_words.h"

void site_to_string(unsigned long, int, char *);
void print_counts(FILE *, int, long **);

int main(int argc, char **argv){
    if(argc != 4){
        printf("usage: ./count_short_words input.fasta[.gz] output.cnt "
               "word_length\n");
        return 1;
    }
    const int len = atoi(argv[3]);
    if(len < 1 || len > 14){
        printf("bad word length, only [1..14] is allowed\n");
        return 1;
    }
    long **counts = count_short_words(argv[1], len);
    FILE *cnt;
    cnt = fopen(argv[2], "wb");
    assert(cnt);
    print_counts(cnt, len, counts);
    fclose(cnt);
    int i;
    for(i = 0; i < len; i ++){
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

void print_counts(FILE *cnt, int len, long **countsp){
    unsigned long max = 1ul << 2;
    char dest[len + 1];
    long *counts;
    int i;
    for(i = 1; i <= len; i ++){
        dest[i] = '\0';
        counts = countsp[len - i];
        fprintf(cnt, "#word_length=%d, total_number=%d\n", i, counts[max]);
        unsigned long site;
        for(site = 0ul; site < max; site ++){
            site_to_string(site, i, dest);
            fprintf(cnt, "%s\t%d\n", dest, counts[site]);
        }
        max <<= 2;
    }
}
