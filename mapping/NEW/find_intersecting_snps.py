import sys
import os
import subprocess
import gzip
import string
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

        self.fastq_single_filename = None
        self.fastq1_filename = None
        self.fastq2_filename = None
        self.fastq1 = None
        self.fastq2 = None
        self.fastq_single = None

        if self.is_paired:
            self.fastq1_filename = self.prefix + ".remap.fq1.gz"
            self.fastq2_filename = self.prefix + ".remap.fq2.gz"
            self.fastq1 = gzip.open(self.fastq1_filename, "wb")
            self.fastq2 = gzip.open(self.fastq2_filename, "wb")
            self.fastq_single_filename = self.prefix + ".remap.single.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wb")
        else:
            self.fastq_single_filename = self.prefix + ".remap.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wb")

        self.input_bam = pysam.Samfile(self.bam_sort_filename, "rb")
        self.keep_bam = pysam.Samfile(self.keep_filename, "wb",
                                      template=self.input_bam)
        self.remap_bam = pysam.Samfile(self.remap_filename, "wb",
                                       template=self.input_bam)
    


        
class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self):
        # number of read matches to reference allele
        self.ref_count = 0
        # number of read matches to alternative allele
        self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        self.other_count = 0

        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0

        # paired reads map to different chromosomes
        self.discard_different_chromosome = 0

        # number of reads discarded because overlap an indel
        self.discard_indel = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0
        
        # number of single reads kept
        self.keep_single = 0
        # number of read pairs kept
        self.keep_pair = 0

        # number of single reads that need remapping
        self.remap_single = 0
        # number of read pairs kept
        self.remap_pair = 0
        

    def write(self, file_handle):
        sys.stderr.write("DISCARD reads:\n"
                         "  improper pair: %d\n"
                         "  different chromosome: %d\n"
                         "  indel: %d\n"
                         "  secondary alignment: %d\n"
                         "  excess allelic combinations: %d\n"
                         "KEEP reads:\n"
                         "  single-end: %d\n"
                         "  pairs: %d\n"
                         "REMAP reads:\n"
                         "  single-end: %d\n"
                         "  pairs: %d\n" %
                         (self.discard_improper_pair,
                          self.discard_different_chromosome,
                          self.discard_indel,
                          self.discard_secondary,
                          self.discard_excess_reads,
                          self.keep_single,
                          self.keep_pair,
                          self.remap_single,
                          self.remap_pair))

        file_handle.write("read SNP ref matches: %d\n" % self.ref_count)
        file_handle.write("read SNP alt matches: %d\n" % self.alt_count)
        file_handle.write("read SNP mismatches: %d\n" % self.other_count)
        
        total = self.ref_count + self.alt_count + self.other_count
        if total > 0:
            mismatch_pct = 100.0 * float(self.other_count) / total
            if mismatch_pct > 10.0:
                sys.stderr.write("WARNING: many read SNP overlaps do not match "
                                 "either allele (%.1f%%). SNP coordinates "
                                 "in input file xmay be incorrect.\n" %
                                 mismatch_pct)
    



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





DNA_COMP = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global DNA_COMP

    if DNA_COMP is None:
        DNA_COMP = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(DNA_COMP)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]



        
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





def write_read(read, snp_tab, snp_idx, read_pos):
    snp_allele1 = [' '] * read.qlen
    snp_allele2 = [' '] * read.qlen

    for (s_idx, r_idx) in zip(snp_idx, read_pos):
        a1 = snp_tab.snp_allele1[s_idx]
        a2 = snp_tab.snp_allele2[s_idx]

        snp_allele1[r_pos-1] = a1
        snp_allele2[r_pos-1] = a2

    sys.stderr.write("READ: %s\n" % read.query)
    sys.stderr.write("A1:   %s\n" % "".join(snp_allele1))
    sys.stderr.write("A2:   %s\n" % "".join(snp_allele2))
    
    
        


def count_ref_alt_matches(read, read_stats, snp_tab, snp_idx, read_pos):
    ref_alleles = snp_tab.snp_allele1[snp_idx]
    alt_alleles = snp_tab.snp_allele2[snp_idx]
    
    for i in range(len(snp_idx)):
        if ref_alleles[i] == read.query[read_pos[i]-1]:
            # read matches reference allele
            read_stats.ref_count += 1
        elif alt_alleles[i] == read.query[read_pos[i]-1]:
            # read matches non-reference allele
            read_stats.alt_count += 1
        else:
            # read matches neither ref nor other
            read_stats.other_count += 1
    


def generate_reads(read_seq, ref_alleles, alt_alleles, read_pos, i):
    """Recursively generate set of reads with all possible combinations
    of alleles (i.e. 2^n combinations where n is the number of snps overlapping
    the reads)
    """
    # create new version of this read with both reference and
    # alternative versions of allele at this index
    idx = read_pos[i]-1
    ref_read = read_seq[:idx] + ref_alleles[i] + read_seq[idx+1:]
    alt_read = read_seq[:idx] + alt_alleles[i] + read_seq[idx+1:]

    if i == len(read_pos)-1:
        # this was the last SNP
        return [ref_read, alt_read]

    # continue recursively with other SNPs overlapping this read
    reads1 = generate_reads(ref_read, ref_alleles, alt_alleles, read_pos, i+1)
    reads2 = generate_reads(alt_read, ref_alleles, alt_alleles, read_pos, i+1)

    return reads1 + reads2
                


def write_fastq(fastq_file, orig_read, new_seqs):
    n_seq = len(new_seqs)
    i = 1
    for new_seq in new_seqs:
        # Give each fastq a new name giving:
        # 1 - the original name of the read
        # 2 - the coordinate that it should map to
        # 3 - the number of the read
        # 4 - the total number of reads being remapped
        name = "%s.%d.%d.%d" % (orig_read.qname, orig_read.pos+1, i, n_seq)
                                       
        fastq_file.write("@%s\n%s\n+%s\n%s\n" %
                         (name, new_seq, name, orig_read.qual))

        i += 1

        
def write_pair_fastq(fastq_file1, fastq_file2, orig_read1, orig_read2,
                     new_pairs):

    n_pair = len(new_pairs)
    i = 1
    for pair in new_pairs:
        # give each fastq record a new name giving:
        # 1 - the original name of the read
        # 2 - the coordinates the two ends of the pair should map to
        # 3 - the number of the read
        # 4 - the total number of reads being remapped

        pos_str = "%d-%d" % (min(orig_read1.pos+1, orig_read2.pos+1),
                             max(orig_read1.pos+1, orig_read2.pos+1))
        
        name = "%s.%s.%d.%d" % (orig_read1.qname, pos_str, i+1, n_pair)
        
        fastq_file1.write("@%s\n%s\n+%s\n%s\n" %
                          (name, pair[0], name, orig_read1.qual))

        rev_seq = revcomp(pair[1])
        fastq_file2.write("@%s\n%s\n+%s\n%s\n" %
                          (name, rev_seq, name, orig_read2.qual))

        i += 1
                         
    

    
def filter_reads(files, max_seqs=MAX_SEQS_DEFAULT, max_snps=5):
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = snptable.SNPTable()
    read_stats = ReadStats()
    read_pair_cache = {}
        
    for read in files.input_bam:
        # TODO: handle paired end reads!!!!
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome
            cur_chrom = files.input_bam.getrname(read.tid)

            if len(read_pair_cache) != 0:
                sys.stderr.write("WARNING: failed to find pairs for %d "
                                 "reads on this chromosome\n" %
                                 len(read_pair_cache))
            read_pair_cache = {}
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)
            snp_filename = "%s/%s.snps.txt.gz" % (files.snp_dir, cur_chrom)

            snp_tab.read_file(snp_filename)

        if read.is_secondary:
            # this is a secondary alignment (i.e. read was aligned more than
            # once and this has align score that <= best score)
            read_stats.discard_secondary += 1
            continue

        if read.is_paired:
            if read.next_reference_name is None:
                # other side of pair not mapped...
                process_single_read(read, read_stats, files, snp_tab, max_seqs)
            elif(read.next_reference_name == cur_chrom or
                 read.next_reference_name == "="):
                # other pair mapped to same chrom
                if not read.is_proper_pair:
                    read_stats.discard_improper_pair += 1
                    continue

                if read.qname in read_pair_cache:
                    # we already saw prev pair, retrieve from cache
                    read1 = read_pair_cache[read.qname]
                    read2 = read
                    del read_pair_cache[read.qname]

                    if read2.next_reference_start != read1.reference_start:
                        sys.stderr.write("WARNING: read pair positions "
                                         "do not match for pair %s\n" %
                                         read.qname)
                    else:
                        process_paired_read(read1, read2, read_stats,
                                            files, snp_tab, max_seqs)
                else:
                    # we need to wait for next pair
                    read_pair_cache[read.qname] = read
            else:
                # other side of pair mapped to different
                # chromosome, discard this read
                read_stats.discard_different_chromosome += 1

        else:
            process_single_read(read, read_stats, files, snp_tab, max_seqs)
                
    
    read_stats.write(sys.stderr)
                     

def process_paired_read(read1, read2, read_stats, files, snp_tab, max_seqs):
    """Checks if either end of read pair overlaps SNPs or indels
    and writes read pair (or generated read pairs) to appropriate 
    output files"""

    new_reads = []    
    for read in (read1, read2):
        # check if either read overlaps SNPs or indels
        # check if read overlaps SNPs or indels
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)
        
        if len(indel_idx) > 0:
            # for now discard this read pair, we want to improve this to handle
            # the indel reads appropriately
            read_stats.discard_indel += 2
            # TODO: add option to handle indels instead of throwing out reads
            return

        if len(snp_idx) > 0:
            ref_alleles = snp_tab.snp_allele1[snp_idx]
            alt_alleles = snp_tab.snp_allele2[snp_idx]

            count_ref_alt_matches(read, read_stats, snp_tab, snp_idx,
                                  snp_read_pos)

            # TODO: limit recursion here by discarding reads that
            # overlap too many SNPs            
            read_seqs = generate_reads(read.query, ref_alleles, alt_alleles,
                                       snp_read_pos, 0)
            
            new_reads.append(read_seqs)
        else:
            # no SNPs or indels overlap this read
            new_reads.append([])
            
    if len(new_reads[0]) == 0 and len(new_reads[1]) == 0:
        # neither read overlapped SNPs or indels
        files.keep_bam.write(read1)
        files.keep_bam.write(read2)
        read_stats.keep_pair += 1
    else:
        # add original version of both sides of pair
        new_reads[0].append(read1.query)
        new_reads[1].append(read2.query)

        if len(new_reads[0]) + len(new_reads[1]) > max_seqs:
            # quit now before generating a lot of read pairs
            read_stats.discard_excess_reads += 2
            return 

        # collect all unique combinations of read pairs
        unique_pairs = set([])
        n_unique_pairs = 0
        for new_read1 in new_reads[0]:
            for new_read2 in new_reads[1]:
                pair = (new_read1, new_read2)
                if pair in unique_pairs:
                    pass
                else:
                    n_unique_pairs += 1
                    if n_unique_pairs > max_seqs:
                        read_stats.discard_excess_reads += 2
                        return
                    unique_pairs.add(pair)

        # remove original read pair, if present
        orig_pair = (read1.query, read2.query)
                                 
        if orig_pair in unique_pairs:
            unique_pairs.remove(orig_pair)
            
        # write read pair to fastqs for remapping
        write_pair_fastq(files.fastq1, files.fastq2, read1, read2, unique_pairs)

        # Write read to 'remap' BAM for consistency with previous
        # implementation of script. Probably not needed and will result in
        # BAM that is not coordinate sorted. Possibly remove this...
        files.remap_bam.write(read1)
        files.remap_bam.write(read2)
        read_stats.keep_pair += 1
        

        
    

    

def process_single_read(read, read_stats, files, snp_tab, max_seqs):
    """Check if a single read overlaps SNPs or indels, and writes
    this read (or generated read pairs) to appropriate output files"""
                
    # check if read overlaps SNPs or indels
    snp_idx, snp_read_pos, \
        indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

    if len(indel_idx) > 0:
        # for now discard this read, we want to improve this to handle
        # the indel reads appropriately
        read_stats.discard_indel += 1
        # TODO: add option to handle indels instead of throwing out reads
        return

    if len(snp_idx) > 0:
        ref_alleles = snp_tab.snp_allele1[snp_idx]
        alt_alleles = snp_tab.snp_allele2[snp_idx]

        count_ref_alt_matches(read, read_stats, snp_tab, snp_idx,
                              snp_read_pos)

        # TODO: limit recursion, by throwing out read
        #       if it overlaps too many SNPs
        read_seqs = generate_reads(read.query, ref_alleles, alt_alleles,
                                   snp_read_pos, 0)

        # make set of unique reads, we don't want to remap
        # duplicates, or the read that matches original
        unique_reads = set(read_seqs)
        if read.query in unique_reads:
            unique_reads.remove(read.query)

        if len(unique_reads) < max_seqs:
            # write read to fastq file for remapping
            write_fastq(files.fastq_single, read, unique_reads)

            # write read to 'to remap' BAM
            # this is probably not necessary with new implmentation
            # but kept for consistency with previous version of script
            files.remap_bam.write(read)
        else:
            # discard read
            read_stats.discard_excess_reads += 1
            return

    else:
        # no SNPs overlap read, write to keep file
        files.keep_bam.write(read)
        read_stats.keep_single += 1
            




def main(bam_filenames, snp_dir, is_paired_end=False,
         is_sorted=False, max_seqs=MAX_SEQS_DEFAULT):

    files = DataFiles(bam_filenames, snp_dir, is_sorted, is_paired_end)
    filter_reads(files, max_seqs=max_seqs)
    
    

    



if __name__ == '__main__':
    options = parse_options()

    main(options.bam_filename, options.snp_dir,
         is_paired_end=options.is_paired_end, is_sorted=options.is_sorted,
         max_seqs=options.max_seqs)
    
