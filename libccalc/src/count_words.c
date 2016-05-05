#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>

#include "count_words.h"

int translate(char nucl){
	switch(nucl){
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		case '\n': case '\t': case ' ': case '-': return -2;
		default: return -1;
	}
}

int skip(gzFile fasta){
	int nucl;
	while((nucl = gzgetc(fasta)) != -1){
		if(nucl == '>'){
			while((nucl = gzgetc(fasta)) != '\n')
				if(nucl == -1)
					return 0;
			continue;
		}
		if(translate(nucl) >= 0){
			gzungetc(nucl, fasta);
			return 1;
		}
	}
	return 0;
}

int initialize_short(gzFile fasta, site_t *sp){
	sp->index = 0;
	int nucl, tnucl;
	sp->val = 0ul;
	while(sp->index < sp->len){
		if((nucl = gzgetc(fasta)) == -1)
			return 0;
		tnucl = translate(nucl);
		if(tnucl < -1)
			continue;
		if(tnucl >= 0){
			sp->val = (sp->val << 2) + tnucl;
			sp->index ++;
		}else{
			gzungetc(nucl, fasta);
			return 0;
		}
	}
	return 1;
}

int initialize_bipart(gzFile fasta, bipart_t *sp){
	unsigned int site = 0u;
	int i, nucl, tnucl;
	sp->index = 0;
	for(i = 0; i < sp->len - 1;){
		if((nucl = gzgetc(fasta)) == -1)
			return 0;
		tnucl = translate(nucl);
		if(tnucl < -1)
			continue;
		if(tnucl >= 0){
			site = (site << 2) + tnucl;
			i ++;
		}else{
			gzungetc(nucl, fasta);
			return 0;
		}
	}
	while(sp->index < sp->size){
		if((nucl = gzgetc(fasta)) == -1)
			return 0;
		tnucl = translate(nucl);
		if(tnucl < -1)
			continue;
		if(tnucl >= 0){
			site = ((site << 2) + tnucl) & sp->mask;
			sp->arr[sp->index] = site;
			sp->index ++;
		}else{
			gzungetc(nucl, fasta);
			return 0;
		}
	}
	return 1;
}

void countup_short(gzFile fasta, site_t *sp, long *counts){
	int nucl, tnucl;
	unsigned long num_index = sp->mask + 1ul;
	while((nucl = gzgetc(fasta)) != -1){
		tnucl = translate(nucl);
		if(tnucl < -1)
			continue;
		if(tnucl >= 0){
			counts[num_index] ++;
			counts[sp->val] ++;
			sp->val = ((sp->val << 2) + tnucl) & sp->mask;
		}else{
			gzungetc(nucl, fasta);
			break;
		}
	}
}

void countup_bipart(gzFile fasta, bipart_t *sp, long *counts){
	int nucl, tnucl;
	unsigned int uhalf, dhalf;
	unsigned long site;
	unsigned long num_index = 1ul << (2*sp->len - sp->ushift + sp->shift);
	int index = sp->size;
	while((nucl = gzgetc(fasta)) != -1){
		tnucl = translate(nucl);
		if(tnucl < -1)
			continue;
		if(tnucl >= 0){
			counts[num_index] ++;
			dhalf = sp->arr[index - 1];
			if(index == sp->size)
				index = 0;
			uhalf = sp->arr[index] >> sp->ushift;
			site = (uhalf << sp->shift) + (dhalf & sp->dmask);
			counts[site] ++;
			sp->arr[index] = ((dhalf << 2) + tnucl) & sp->mask;
			index ++;
		}else{
			gzungetc(nucl, fasta);
			break;
		}
	}
	sp->index = index;
}

void countup_short_res(site_t *sp, long **countsp){
	unsigned long mask = sp->mask;
	while(sp->index < sp->len){
		mask >>= 2;
		countsp ++;
		sp->index ++;
	}
	unsigned long num_index;
	while(mask){
		num_index = mask + 1ul;
		(*countsp)[num_index] ++;
		(*countsp)[sp->val & mask] ++;
		countsp ++;
		mask >>= 2;
	}
}

void countup_bipart_res(bipart_t *sp, int real_uindex, long **countsp){
	unsigned long mask = sp->dmask;
	int shift = sp->shift;
	int uindex = sp->index % sp->size;
	int dindex = (sp->index + sp->size - 1) % sp->size;
	while(uindex != real_uindex){
		mask >>= 2;
		shift -= 2;
		countsp ++;
		uindex ++;
	}
	unsigned long site;
	unsigned int uhalf;
	unsigned int dhalf = sp->arr[dindex];
	unsigned long num_index = (1ul + sp->mask >> sp->ushift) << shift;
	while(mask){
		(*countsp)[num_index] ++;
		num_index >>= 2;
		uindex %= sp->size;
		uhalf = sp->arr[uindex] >> sp->ushift;
		site = (uhalf << shift) + (dhalf & mask);
		(*countsp)[site] ++;
		countsp ++;
		uindex ++;
		mask >>= 2;
		shift -= 2;
	}
}

void countup_subsites(int len, int pos, long **countsp){
	long *src = countsp[0];
	long *dst = countsp[1];
	int depth = 0;
	int max_depth = len - pos - 1;
	unsigned long site;
	unsigned long max = 1ul << len * 2;
	while(depth < max_depth){
		for(site = 0ul; site < max; site ++)
			dst[site >> 2] += src[site];
		dst[max >> 2] += src[max];
		max >>= 2;
		depth ++;
		src = countsp[depth];
		dst = countsp[depth + 1];
	}
}

bipart_t make_bpsite(int gap, int ulen, int dlen, int minlen, int maxlen){
	bipart_t bpsite;
	bpsite.size = gap + minlen + 1;
	bpsite.len = maxlen;
	bpsite.mask = (1u << 2 * maxlen) - 1u;
	bpsite.shift = dlen * 2;
	bpsite.dmask = (1u << bpsite.shift) - 1u;
	bpsite.ushift = (maxlen - ulen) * 2;
	bpsite.arr = (unsigned int *)calloc(bpsite.size, sizeof(unsigned int));
	bpsite.index = 0;
	return bpsite;
}

long **count_short_words(char *filename, int len){
	
	gzFile fasta;
	fasta = gzopen(filename, "rb");
	assert(fasta);
	
	long **counts;
	counts = (long **)calloc(len, sizeof(long *));
	int i;
	unsigned long num = (1ul << len * 2);
	for(i = 0; i < len; i ++){
		counts[i] = (long *)calloc(num + 1, sizeof(long));// +1 for total
		num >>= 2;
	}
	assert(counts);
	
	site_t site;
	site.len = len;
	site.mask = (1ul << len * 2) - 1ul;
	site.index = 0;
	while(1){
		if(skip(fasta) == 0)
			break;
		if(initialize_short(fasta, &site)){
			countup_short(fasta, &site, counts[0]);
		}
		countup_short_res(&site, counts);
	}
	countup_subsites(len, 0, counts);
	gzclose(fasta);
	return counts;
}

long **count_bipart_words(char *filename, int len, int pos, int gap){
	
	gzFile fasta;
	fasta = gzopen(filename, "rb");
	assert(fasta);
	
	int ulen = pos;
	int dlen = len - ulen;
	int minlen = ulen < dlen ? ulen : dlen;
	int maxlen = len - minlen;
	
	long **counts;
	counts = (long **)calloc(dlen, sizeof(long *));
	int i;
	unsigned long num = (1ul << len * 2);
	for(i = 0; i < dlen; i ++){
		counts[i] = (long *)calloc(num + 1, sizeof(long));// +1 for total
		num >>= 2;
	}
	assert(counts);
	
	bipart_t bpsite;
	bpsite = make_bpsite(gap, ulen, dlen, minlen, maxlen);
	while(1){
		if(skip(fasta) == 0)
			break;
		int uindex = bpsite.size;
		if(initialize_bipart(fasta, &bpsite)){
			countup_bipart(fasta, &bpsite, counts[0]);
			uindex = bpsite.index % bpsite.size;
		}
		countup_bipart_res(&bpsite, uindex, counts);
	}
	free(bpsite.arr);
	countup_subsites(len, ulen, counts);
	gzclose(fasta);
	return counts;
}
