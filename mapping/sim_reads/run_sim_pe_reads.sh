#!/bin/bash

WASP=$HOME/proj/WASP

python sim_pe_reads.py --seq $WASP/test_data/seq.h5 \
       --n_reads 10000 \
       --read_len 36 \
       --hap_file $WASP/example_data/genotypes/chr22.hg19.haplotype.txt.gz \
       --out_fastq1 $WASP/example_data/sim_pe_reads1.fastq.gz \
       --out_fastq2 $WASP/example_data/sim_pe_reads2.fastq.gz

