#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

#define VERSION "0.0.2"

int main(int argc,char** argv) {
    fprintf(stderr, "Usage is: %s inputBam N outputBam1 outputBam2 ... ouotputBamN (version=" VERSION ")\n", argv[0]);
    fprintf(stderr, "  N is the number of partitions of inputBam into outputBams. \n");
    //fprintf(stderr, "Usage is: %s inputBam prob1 outputBam1 prob2 outputBam2 (version=" VERSION ")\n", argv[0]);
    //fprintf(stderr, "  prob1 and prob2 are the probabilities that the reads in inputBam are in outputBam1 and outputBam2, respectively. \n");
    
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;
    
    if (argc < 5) return -1;
    in = sam_open(argv[1], "r");
    if (in==NULL) return -1;
    if ((header = sam_hdr_read(in)) == 0) return -1;
    b = bam_init1();
    
    int nbams = atoi(argv[2]);
	samFile** outs = calloc(nbams, sizeof(samFile*));
	for (int i = 0; i < nbams; i++) {
		outs[i] = sam_open(argv[3+i], "wb2");
		assert(outs[i] != NULL);
		int sam_write_ret = sam_hdr_write(outs[i], header);
		assert(sam_write_ret >= 0);
	}
	
    while (sam_read1(in, header, b) >= 0) {
        const uint32_t subsam_seed = 256UL-1UL;
        /* https://github.com/samtools/samtools/blob/52bd699e0a75b6d39504d1f1bdb38dacdc903921/sam_view.c#L93
         * uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
         * if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
         * */
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ subsam_seed);
        double randfrac = (double)(k&0xffffff) / 0x1000000;
        
		for (int i = 0; i < nbams; i++) {
			if (randfrac < (double)(i+1)/(double)nbams) {
				int sam_write_ret = sam_write1(outs[i], header, b);
				assert(sam_write_ret);
				break;
			}
		}
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
	
	for (int i = 0; i < nbams; i++) {
		sam_close(outs[i]);
	}
	free(outs);

    return 0;
}

