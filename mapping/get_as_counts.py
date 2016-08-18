

import sys
import argparse
import numpy as np
import pysam

import util
import snptable



def write_results(out_f, chrom_name, snp_tab, ref_matches,
                  alt_matches, oth_matches):
    
    for i in range(snp_tab.n_snp):
        out_f.write("%s %d %s %s %d %d %d\n" % (chrom_name, snp_tab.snp_pos[i],
                                             snp_tab.snp_allele1[i],
                                             snp_tab.snp_allele2[i],
                                             ref_matches[i],
                                             alt_matches[i],
                                             oth_matches[i]))


def write_header(out_f):
    out_f.write("CHROM SNP.POS REF.ALLELE ALT.ALLELE OTHER.COUNT "
                "ALT.COUNT OTH.COUNT\n")
    


def parse_options():
    parser = argparse.ArgumentParser(description="This script outputs "
                                     "allele-specific counts for SNPs, using "
                                     "reads from the provided BAM file. "
                                     "Currently indels are not output and "
                                     "chromosomes with no mapped reads "
                                     "are skipped. Output "
                                     "is written to stdout, with a single "
                                     "header row and the following "
                                     "columns: <chromosome> <snp_position> "
                                     "<ref_allele> <alt_allele> "
                                     "<ref_allele_count> <alt_allele_count> "
                                     "<other_count>. Reads that overlap "
                                     "multiple SNPs will be counted multiple "
                                     "times in the output (this behavior "
                                     "differs from the "
                                     "extract_haplotype_read_counts.py "
                                     "script).")

    parser.add_argument("bam_filename",
                        help="Sorted BAM file containing reads")
    
    parser.add_argument("snp_dir",
                        help=("Directory containing SNPs. "
                              "This directory should contain one file per "
                              "chromosome named like chr<#>.snps.txt.gz. "
                              "Each file should contain 3 columns: position "
                              "RefAllele AltAllele"))


    return parser.parse_args()
                        
    

def main(bam_filename, snp_dir):
    out_f = sys.stdout
    
    bam = pysam.Samfile(bam_filename)
        
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = snptable.SNPTable()
    read_pair_cache = {}

    # keep track of number of ref matches, non-ref matches, and other
    # for each SNP
    snp_ref_match = None
    snp_alt_match = None
    snp_other_match = None
    
    for read in bam:
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome

            if cur_chrom:
                # write out results from last chromosome
                write_results(out_f, cur_chrom, snp_tab, snp_ref_match,
                              snp_alt_match, snp_oth_match)
            
            cur_chrom = bam.getrname(read.tid)
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)
            snp_filename = "%s/%s.snps.txt.gz" % (snp_dir, cur_chrom)
            
            # clear SNP table and results, read next SNP filename
            snp_tab.read_file(snp_filename)
            snp_ref_match = np.zeros(snp_tab.n_snp, dtype=np.int16)
            snp_alt_match = np.zeros(snp_tab.n_snp, dtype=np.int16)
            snp_oth_match = np.zeros(snp_tab.n_snp, dtype=np.int16)

            sys.stderr.write("read %d SNPs\n" % snp_tab.n_snp)
                 
        if read.is_secondary:
            # this is a secondary alignment (i.e. read was aligned more than
            # once and this has align score that <= best score)
            continue

        # loop over all SNP that overlap this read
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

        for snp_i, read_pos in zip(snp_idx, snp_read_pos):
            snp_pos = snp_tab.snp_pos[snp_i]
            ref_allele = snp_tab.snp_allele1[snp_i]
            alt_allele = snp_tab.snp_allele2[snp_i]
            
            if ref_allele == read.query[read_pos-1]:
                snp_ref_match[snp_i] += 1
            elif alt_allele == read.query[read_pos-1]:
                snp_alt_match[snp_i] += 1
            else:
                snp_oth_match[snp_i] += 1

    if cur_chrom:
        # write results for final chromosome
        write_results(out_f, cur_chrom, snp_tab, snp_ref_match,
                      snp_alt_match, snp_oth_match)


    


if __name__ == "__main__":
    options = parse_options()

    main(options.bam_filename, options.snp_dir)

    
