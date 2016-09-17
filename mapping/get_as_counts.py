import sys
import argparse
import numpy as np
import pysam

import util
import snptable

import tables

import os

def write_results(out_f, chrom_name, snp_tab, ref_matches,
                  alt_matches, oth_matches, geno_sample):

    haps = None
    has_haps = False
    
    if geno_sample:
        # get index for this sample in the haplotype table
        samp_idx_dict = dict(zip(snp_tab.samples,
                                 range(len(snp_tab.samples))))

        if geno_sample in samp_idx_dict:
            idx = samp_idx_dict[geno_sample]
            geno_hap_idx = np.array([idx*2, idx*2+1], dtype=np.int)
            haps = snp_tab.haplotypes[:,geno_hap_idx]
            has_haps = True
            sys.stderr.write("geno_hap_idx: %s\n" % repr(geno_hap_idx))
        else:
            sys.stderr.write("WARNING: sample %s is not present for "
                             "chromosome %s\n" % (geno_sample, chrom_name))
            haps = None
            has_haps = False

    for i in range(snp_tab.n_snp):
        if has_haps:
            geno_str = "%d|%d" % (haps[i, 0], haps[i, 1])
        else:
            geno_str = "NA"
        out_f.write("%s %d %s %s %s %d %d %d\n" %
                    (chrom_name, snp_tab.snp_pos[i],
                     snp_tab.snp_allele1[i], snp_tab.snp_allele2[i],
                     geno_str, ref_matches[i], alt_matches[i],
                     oth_matches[i]))


def write_header(out_f):
    out_f.write("CHROM SNP.POS REF.ALLELE ALT.ALLELE GENOTYPE REF.COUNT "
                "ALT.COUNT OTHER.COUNT\n")


    
def parse_samples(samples_str):
    """Gets list of samples from --samples argument. This may be 
    a comma-delimited string or a path to a file. If a file is provided 
    then the first column of the file is assumed to be the sample name"""

    if samples_str is None:
        return None
        
    # first check if this is a path to a file
    if os.path.exists(samples_str) and not os.path.isdir(samples_str):
        samples = []

        if util.is_gzipped(samples_str):
            f = gzip.open(samples_str)
        else:
            f = open(samples_str)

        for line in f:
            # assume first token in line is sample name
            samples.append(line.split()[0])

        sys.stderr.write("read %d sample names from file '%s'\n" %
                         (len(samples), samples_str))
                    
        f.close()
    else:    
        # otherwise assume comma-delimited string
        if ("," not in samples_str and len(samples_str) > 15) \
           or ("/" in samples_str):
            sys.stderr.write("WARNING: --samples argument (%s) "
                             "does not look like sample name "
                             "but is not path to valid file. "
                             "Assuming it is a sample name anyway."
                             % samples_str)

        samples = samples_str.split(",")
        sys.stderr.write("SAMPLES: %s\n"% repr(samples))


    return samples


    


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
                                     "<ref_allele> <alt_allele> <genotype> "
                                     "<ref_allele_count> <alt_allele_count> "
                                     "<other_count>. Reads that overlap "
                                     "multiple SNPs will be counted multiple "
                                     "times in the output (this behavior "
                                     "differs from the "
                                     "extract_haplotype_read_counts.py "
                                     "script).")


    parser.add_argument("--snp_dir", action='store', 
                        help=("Directory containing SNP text files "
                              "This directory should contain one file per "
                              "chromosome named like chr<#>.snps.txt.gz. "
                              "Each file should contain 3 columns: position "
                              "RefAllele AltAllele"),
                        default=None)
        

    parser.add_argument("--snp_tab",
                        help="Path to HDF5 file to read SNP information "
                        "from. Each row of SNP table contains SNP name "
                        "(rs_id), position, allele1, allele2.",
                        metavar="SNP_TABLE_H5_FILE",
                        default=None)
    
    parser.add_argument("--snp_index",
                        help="Path to HDF5 file containing SNP index. The "
                        "SNP index is used to convert the genomic position "
                        "of a SNP to its corresponding row in the haplotype "
                        "and snp_tab HDF5 files.",
                        metavar="SNP_INDEX_H5_FILE",
                        default=None)
    
    parser.add_argument("--haplotype",
                        help="Path to HDF5 file to read phased haplotypes "
                        "from. When generating alternative reads "
                        "use known haplotypes from this file rather "
                        "than all possible allelic combinations.",
                        metavar="HAPLOTYPE_H5_FILE",
                        default=None)

    parser.add_argument("--samples",
                        help="Use only haplotypes and SNPs that are "
                        "polymorphic in these samples. "
                        "SAMPLES can either be a comma-delimited string "
                        "of sample names or a path to a file with one sample "
                        "name per line (file is assumed to be "
                        "whitespace-delimited and first column is assumed to "
                        "be sample name). Sample names should match those "
                        "present in the haplotype HDF5 file. Samples are "
                        "ignored if no haplotype file is provided.",
                        metavar="SAMPLES", default=None)


    parser.add_argument("--genotype_sample",
                        metavar="GENO_SAMPLE",
                        help="output genotypes for sample with name "
                        "GENO_SAMPLE alongside allele-specific counts. "
                        "GENO_SAMPLE must match one "
                        "of the names present in the haplotype HDF5 file. "
                        "If the --samples argument is provided then "
                        "GENO_SAMPLE must also be one of the specified "
                        "samples. If --genotype_sample is "
                        "not provided or the GENO_SAMPLE does not match any "
                        "of the samples in haplotype file then NA is "
                        "output for genotype.", default=None)
        
    parser.add_argument("bam_filename", action='store',
                        help="Coordinate-sorted input BAM file "
                        "containing mapped reads.")


    options = parser.parse_args()
    
    if options.snp_dir:
        if(options.snp_tab or options.snp_index or options.haplotype):
            parser.error("expected --snp_dir OR (--snp_tab, --snp_index and "
                         "--haplotype) arguments but not both")
    else:
        if not (options.snp_tab and options.snp_index and options.haplotype):
            parser.error("either --snp_dir OR (--snp_tab, "
                         "--snp_index AND --haplotype) arguments must be "
                         "provided")
     
    return options
                        
    

def main(bam_filename, snp_dir=None, snp_tab_filename=None,
         snp_index_filename=None, haplotype_filename=None, samples=None,
         geno_sample=None):

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

    
    if geno_sample and not haplotype_filename:
        sys.stderr.write("WARNING: cannot obtain genotypes for sample "
                         "%s without --haplotype argument\n")
        geno_sample = None

    sys.stderr.write("GENOTYPE_SAMPLE: %s\n" % geno_sample)

    if snp_tab_filename:
        if (not snp_index_filename) or (not haplotype_filename):
            raise ValueError("--snp_index and --haplotype must be provided "
                             "if --snp_tab is provided")
        snp_tab_h5 = tables.openFile(snp_tab_filename, "r")
        snp_index_h5 = tables.openFile(snp_index_filename, "r")
        hap_h5 = tables.openFile(haplotype_filename, "r")
    else:
        snp_tab_h5 = None
        snp_index_h5 = None
        hap_h5 = None
        
    for read in bam:
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome

            if cur_chrom:
                # write out results from last chromosome
                write_results(out_f, cur_chrom, snp_tab, snp_ref_match,
                              snp_alt_match, snp_oth_match, geno_sample)
            
            cur_chrom = bam.getrname(read.tid)
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)

            # read SNPs for next chromomsome
            if snp_tab_h5:
                # read SNPs from HDF5 files, reduce to set that are
                # polymorphic in specified samples
                snp_tab.read_h5(snp_tab_h5, snp_index_h5, hap_h5,
                                cur_chrom, samples=samples)
            elif snp_dir:
                # read SNPs from text file
                snp_filename = "%s/%s.snps.txt.gz" % (snp_dir, cur_chrom)
                snp_tab.read_file(snp_filename)
            else:
                raise ValueError("--snp_dir OR (--snp_tab, --snp_index, "
                                 "and --hap_h5) must be defined")

            sys.stderr.write("read %d SNPs\n" % snp_tab.n_snp)
            
            # clear SNP table and results             
            snp_ref_match = np.zeros(snp_tab.n_snp, dtype=np.int16)
            snp_alt_match = np.zeros(snp_tab.n_snp, dtype=np.int16)
            snp_oth_match = np.zeros(snp_tab.n_snp, dtype=np.int16)
                
                 
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
            
            if ref_allele == read.query_sequence[read_pos-1]:
                snp_ref_match[snp_i] += 1
            elif alt_allele == read.query_sequence[read_pos-1]:
                snp_alt_match[snp_i] += 1
            else:
                snp_oth_match[snp_i] += 1

    if cur_chrom:
        # write results for final chromosome
        write_results(out_f, cur_chrom, snp_tab, snp_ref_match,
                      snp_alt_match, snp_oth_match, geno_sample)


    


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()
    samples = parse_samples(options.samples)

    
    main(options.bam_filename, 
         snp_dir=options.snp_dir,
         snp_tab_filename=options.snp_tab,
         snp_index_filename=options.snp_index,
         haplotype_filename=options.haplotype,
         samples=samples, geno_sample=options.genotype_sample)
    

    
