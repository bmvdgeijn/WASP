import sys
import os
import subprocess
import gzip
import argparse
import numpy as np

import pysam

import snptable



MAX_SEQS_DEFAULT = 10
            

                
                
                
                


        
class DataFiles(object):
    """Object to hold names and filehandles for all input / output 
    datafiles"""
    
    def __init__(self, bam_filename, snp_dir, is_sorted, is_paired):
        self.is_paired = is_paired
        self.bam_filename = bam_filename

        self.snp_dir = snp_dir
        
        name_split = self.bam_filename.split(".")
        if len(name_split) > 1:
           self.prefix = ".".join(name_split[:-1])
        else:
            self.prefix = name_split[0]
        
        # TODO: could allow names of output files to be specified
        # on command line rather than appending name to prefix

        if not is_sorted:
            sort_bam(self.bam_filename, self.prefix)
            self.bam_sort_filename = self.prefix + ".sort.bam"
        else:
            self.bam_sort_filename = self.bam_filename

        self.keep_filename = self.prefix + ".keep.bam"
        self.remap_filename = self.prefix + ".to.remap.bam"

        self.fastq1_filename = None
        self.fastq2_filename = None
        self.fastq1 = None
        self.fastq2 = None

        if self.is_paired:
            self.fastq1_filename = prefix + ".remap.fq1.gz"
            self.fastq2_filename = prefix + ".remap.fq2.gz"
            self.fastq1 = gzip.open(self.fastq1_filename, "wb")
            self.fastq2 = gzip.open(self.fastq2_filename, "wb")
        else:
            self.fastq1_filename = self.prefix + ".remap.fq.gz"
            self.fastq1 = gzip.open(self.fastq1_filename, "wb")

        self.input_bam = pysam.Samfile(self.bam_sort_filename, "rb")
        self.keep_bam = pysam.Samfile(self.keep_filename, "wb",
                                      template=self.input_bam)
        self.remap_bam = pysam.Samfile(self.remap_filename, "wb",
                                       template=self.input_bam)
    


        

    



def parse_options():
    parser=argparse.ArgumentParser()

    parser.add_argument("--is_paired_end", "-p", action='store_true',
                        dest='is_paired_end', 
                        default=False,
                        help=("Indicates that reads are paired-end "
                              "(default is single)."))
    
    parser.add_argument("--is_sorted", "-s", action='store_true',
                        dest='is_sorted', 
                        default=False,
                        help=('Indicates that the input bam file'
                              ' is coordinate sorted (default '
                              'is False).'))
    
    parser.add_argument("--max_seqs", type=int, default=MAX_SEQS_DEFAULT,
                        help="The maximum number of sequences with different "
                        "allelic combinations to consider remapping "
                        "(default=%d). Reads with more allelic combinations "
                        "than MAX_SEQs are discarded" % MAX_SEQS_DEFAULT)
                        
    parser.add_argument("bam_filename", action='store',
                        help="Coordinate sorted bam file containing reads.")

    parser.add_argument("snp_dir", action='store', 
                        help=("Directory containing SNPs "
                              "This directory should contain one file per "
                              "chromosome named like chr<#>.snps.txt.gz. "
                              "Each file should contain 3 columns: position "
                              "RefAllele AltAllele"))
                        
    return parser.parse_args()



        
def sort_bam(input_bam, output_prefix):
    """Calls samtools sort on input_bam filename and writes to
    output_bam. Takes into account that the command line arguments 
    for samtools sort have changed between versions."""

    output_bam = output_prefix + ".sort.bam"
    
    # first try new way of using samtools sort
    failed = False
    cmd = "samtools sort -o " + output_bam + " " + input_bam
    sys.stderr.write("running command: %s\n" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        sys.stderr.write("samtools sort command failed:\n%s\n" %
                         str(e))
        failed = True
    if not os.path.exists(output_bam):
        sys.stderr.write("output file %s does not exist\n" % output_bam)
        failed = True
        
    if failed:
        # OLD way of calling samtools (changed in newer versions)
        sys.stderr.write("samtools sort command failed, trying old samtools "
                         "syntax\n")
        
        cmd = "samtools sort " + input_bam + " " + output_prefix
        sys.stderr.write("running command: %s\n" % cmd)

        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            sys.stderr.write("samtools sort command failed:\n%s\n" %
                             str(e))
            exit(1)
        
        if not os.path.exists(paths.sorted_output_bam):
            raise IOError("Failed to create sorted BAM file '%s'" %
                          paths.sorted_output_bam)





def write_read(read, snp_tab, snp_idx, read_idx):
    snp_allele1 = [' '] * read.qlen
    snp_allele2 = [' '] * read.qlen

    for (s_idx, r_idx) in zip(snp_idx, read_idx):
        a1 = snp_tab.snp_allele1[s_idx]
        a2 = snp_tab.snp_allele2[s_idx]

        snp_allele1[r_idx] = a1
        snp_allele2[r_idx] = a2

    sys.stderr.write("READ: %s\n" % read.query)
    sys.stderr.write("A1:   %s\n" % "".join(snp_allele1))
    sys.stderr.write("A2:   %s\n" % "".join(snp_allele2))
    
    
        


def count_ref_alt_matches(read, snp_tab, snp_idx, read_idx):
    ref_alleles = snp_tab.snp_allele1[snp_idx]
    alt_alleles = snp_tab.snp_allele2[snp_idx]
    
    ref_count = 0
    alt_count = 0
    other_count = 0

    for i in range(len(snp_idx)):
        if ref_alleles[i] == read.query[read_idx[i]]:
            # read matches reference allele
            ref_count += 1
        elif alt_alleles[i] == read.query[read_idx[i]]:
            # read matches non-reference allele
            alt_count += 1
        else:
            # read matches neither ref nor other
            other_count += 1

    return ref_count, alt_count, other_count
    


def generate_reads(read_seq, ref_alleles, alt_alleles, read_idx, i):
    """Recursively generate set of reads with all possible combinations
    of alleles (i.e. 2^n combinations where n is the number of snps overlapping
    the reads)
    """
    # create new version of this read with both reference and
    # alternative versions of allele at this index
    idx = read_idx[i]
    ref_read = read_seq[:idx] + ref_alleles[i] + read_seq[idx+1:]
    alt_read = read_seq[:idx] + alt_alleles[i] + read_seq[idx+1:]

    if i == len(read_idx)-1:
        # this was the last SNP
        return [ref_read, alt_read]

    # continue recursively with other SNPs overlapping this read
    reads1 = generate_reads(ref_read, ref_alleles, alt_alleles, read_idx, i+1)
    reads2 = generate_reads(alt_read, ref_alleles, alt_alleles, read_idx, i+1)

    return reads1 + reads2
                


def write_fastq(fastq_file, orig_read, new_seqs):
    for new_seq in new_seqs:
        # TODO: give each fastq a new name giving:
        # 1 - the original name of the read
        # 2 - the coordinate that it should map to
        # 3 - the number of the read
        # 4 - the total number of reads being remapped
        name = orig_read.qname
        fastq_file.write("@%s\n%s\n+%s\n%s\n" %
                         (name, new_seq, name, orig_read.qual))
                         
    

    
def filter_reads(files, max_seqs=MAX_SEQS_DEFAULT, max_snps=5):
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = snptable.SNPTable()

    n_ref_match = 0
    n_alt_match = 0
    n_other = 0
    
    
    for read in files.input_bam:
        # TODO: handle paired end reads!!!!
        
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome
            cur_chrom = files.input_bam.getrname(read.tid)

            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)
            snp_filename = "%s/%s.snps.txt.gz" % (files.snp_dir, cur_chrom)

            snp_tab.read_file(snp_filename)

        # check if read overlaps SNPs or indels
        snp_idx, read_idx = snp_tab.get_overlapping_snps(read)

        if len(snp_idx) > 0:
            ref_alleles = snp_tab.snp_allele1[snp_idx]
            alt_alleles = snp_tab.snp_allele2[snp_idx]

            ref, alt, other = count_ref_alt_matches(read, snp_tab,
                                                    snp_idx, read_idx)
            n_ref_match += ref
            n_alt_match += alt
            n_other += other
            
            
            # TODO: limit recursion, by throwing out read
            #       if it overlaps too many SNPs
            read_seqs = generate_reads(read.query, ref_alleles, alt_alleles,
                                       read_idx, 0)
            
            # make set of unique reads, we don't want to remap
            # duplicates, or the read that matches original
            unique_reads = set(read_seqs)
            if read.query in unique_reads:
                unique_reads.remove(read.query)
             
            # write read to fastq file for remapping
            write_fastq(files.fastq1, read, unique_reads)

            # write read to 'to remap' BAM
            # this is probably not necessary with new implmentation
            # but kept for consistency with previous version of script
            files.remap_bam.write(read)

        else:
            # no SNPs overlap read, write to keep file
            files.keep_bam.write(read)

    sys.stderr.write("read SNP ref matches: %d\n" % n_ref_match)
    sys.stderr.write("read SNP alt matches: %d\n" % n_alt_match)
    sys.stderr.write("read SNP mismatches: %d\n" % n_other)

    total = n_other + n_ref_match + n_alt_match
    if total > 0:
        mismatch_pct = 100.0 * float(n_other) / total

        if mismatch_pct > 20.0:
            sys.stderr.write("WARNING: many read SNP overlaps do not match "
                             "either allele (%.1f%%). SNP coordinates "
                             "in input file xmay be incorrect.\n" % mismatch_pct)
    
    
                     
        



def main(bam_filenames, snp_dir, is_paired_end=False,
         is_sorted=False, max_seqs=MAX_SEQS_DEFAULT):

    files = DataFiles(bam_filenames, snp_dir, is_sorted, is_paired_end)
    filter_reads(files, max_seqs=max_seqs)
    
    

    



if __name__ == '__main__':
    options = parse_options()

    main(options.bam_filename, options.snp_dir,
         is_paired_end=options.is_paired_end, is_sorted=options.is_sorted,
         max_seqs=options.max_seqs)
    
