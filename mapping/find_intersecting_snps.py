import sys
import os
import gzip
import argparse
import numpy as np

import pysam

import util
import snptable

import tables

MAX_SEQS_DEFAULT = 64
MAX_SNPS_DEFAULT = 6


class DataFiles(object):
    """Object to hold names and filehandles for all input / output 
    datafiles"""
    
    def __init__(self, bam_filename, is_sorted, is_paired,
                 output_dir=None, snp_dir=None,
                 snp_tab_filename=None, snp_index_filename=None,
                 haplotype_filename=None, samples=None):
        # flag indicating whether reads are paired-end
        self.is_paired = is_paired
        
        # prefix for output files
        self.prefix = None

        # name of input BAM filename
        self.bam_filename = bam_filename        
        # name of sorted input bam_filename
        # (new file is created if input file is not
        #  already sorted)
        self.bam_sort_filename = None
        # pysam file handle for input BAM
        self.input_bam = None

        # name of output keep and to.remap BAM files
        self.keep_filename = None
        self.remap_filename = None

        # pysam file handles for output BAM filenames
        self.keep_bam = None
        self.remap_bam = None

                
        # name of output fastq files
        self.fastq_single_filename = None
        self.fastq1_filename = None
        self.fastq2_filename = None
        self.fastq1 = None
        self.fastq2 = None
        self.fastq_single = None

        # name of directory to read SNPs from
        self.snp_dir = snp_dir

        # paths to HDF5 files to read SNP info from
        self.snp_tab_filename = snp_tab_filename
        self.snp_index_filename = snp_index_filename
        self.haplotype_filename = haplotype_filename

        if self.snp_tab_filename:
            self.snp_tab_h5 = tables.openFile(snp_tab_filename, "r")
            self.snp_index_h5 = tables.openFile(snp_index_filename, "r")
            self.hap_h5 = tables.openFile(haplotype_filename, "r")
        else:
            self.snp_tab_h5 = None
            self.snp_index_h5 = None
            self.hap_h5 = None

            
        # separate input directory and bam filename
        tokens = self.bam_filename.split("/")
        bam_dir = "/".join(tokens[:-1])
        filename = tokens[-1]

        if output_dir is None:
            # if no output dir specified, use same directory as input
            # bam file
            output_dir = bam_dir
        else:
            if output_dir.endswith("/"):
                # strip trailing '/' from output dir name
                output_dir = output_dir[:-1]
                
        name_split = filename.split(".")
        if len(name_split) > 1:
           self.prefix = output_dir + "/" + ".".join(name_split[:-1])
        else:
            self.prefix = output_dir + "/" + name_split[0]
            
        # TODO: could allow names of output files to be specified
        # on command line rather than appending name to prefix
        sys.stderr.write("prefix: %s\n" % self.prefix)
        
        if not is_sorted:
            util.sort_bam(self.bam_filename, self.prefix)
            self.bam_sort_filename = self.prefix + ".sort.bam"
        else:
            self.bam_sort_filename = self.bam_filename

        self.keep_filename = self.prefix + ".keep.bam"
        self.remap_filename = self.prefix + ".to.remap.bam"

        sys.stderr.write("reading reads from:\n  %s\n" %
                         self.bam_sort_filename)
        
        sys.stderr.write("writing output files to:\n")
        
        if self.is_paired:
            self.fastq1_filename = self.prefix + ".remap.fq1.gz"
            self.fastq2_filename = self.prefix + ".remap.fq2.gz"
            self.fastq1 = gzip.open(self.fastq1_filename, "wb")
            self.fastq2 = gzip.open(self.fastq2_filename, "wb")
            self.fastq_single_filename = self.prefix + ".remap.single.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wb")
            sys.stderr.write("  %s\n  %s\n  %s\n" %
                             (self.fastq1_filename,
                              self.fastq2_filename,
                              self.fastq_single_filename))
            
        else:
            self.fastq_single_filename = self.prefix + ".remap.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wb")
            sys.stderr.write("  %s\n" % (self.fastq_single_filename))

        self.input_bam = pysam.Samfile(self.bam_sort_filename, "rb")
        self.keep_bam = pysam.Samfile(self.keep_filename, "wb",
                                      template=self.input_bam)
        self.remap_bam = pysam.Samfile(self.remap_filename, "wb",
                                       template=self.input_bam)
        sys.stderr.write("  %s\n  %s\n" % (self.keep_filename,
                                           self.remap_filename))


    
        
    def close(self):
        """close open filehandles"""
        filehandles = [self.keep_bam, self.remap_bam, self.fastq1,
                       self.fastq2, self.fastq_single,
                       self.snp_tab_h5, self.snp_index_h5,
                       self.hap_h5]

        for fh in filehandles:
            if fh:
                fh.close()

        
class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self):
        # number of read matches to reference allele
        self.ref_count = 0
        # number of read matches to alternative allele
        self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        self.other_count = 0

        # number of reads discarded becaused not mapped
        self.discard_unmapped = 0
        
        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0

        # number of reads discarded because mate unmapped
        self.discard_mate_unmapped = 0

        # paired reads map to different chromosomes
        self.discard_different_chromosome = 0

        # number of reads discarded because overlap an indel
        self.discard_indel = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # number of reads discarded because of too many overlapping SNPs
        self.discard_excess_snps = 0
        
        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0
        
        # reads where we expected to see other pair, but it was missing
        # possibly due to read-pairs with different names
        self.discard_missing_pair = 0
        
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
                         "  unmapped: %d\n"
                         "  mate unmapped: %d\n"
                         "  improper pair: %d\n"
                         "  different chromosome: %d\n"
                         "  indel: %d\n"
                         "  secondary alignment: %d\n"
                         "  excess overlapping snps: %d\n"
                         "  excess allelic combinations: %d\n"
                         "  missing pairs (e.g. mismatched read names): %d\n"
                         "KEEP reads:\n"
                         "  single-end: %d\n"
                         "  pairs: %d\n"
                         "REMAP reads:\n"
                         "  single-end: %d\n"
                         "  pairs: %d\n" %
                         (self.discard_unmapped,
                          self.discard_mate_unmapped,
                          self.discard_improper_pair,
                          self.discard_different_chromosome,
                          self.discard_indel,
                          self.discard_secondary,
                          self.discard_excess_snps,
                          self.discard_excess_reads,
                          self.discard_missing_pair,
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
                                 "in input file may be incorrect.\n" %
                                 mismatch_pct)
    



def parse_options():
    
    parser = argparse.ArgumentParser(description="Looks for SNPs and indels "
                                     "overlapping reads. If a read overlaps "
                                     "SNPs, alternative versions of the read "
                                     "containing different alleles are created "
                                     "and written to files for remapping. "
                                     "Reads that do not overlap SNPs or indels "
                                     "are written to a 'keep' BAM file."
                                     "Reads that overlap indels are presently "
                                     "discarded.")
                                   

    parser.add_argument("--is_paired_end", "-p", action='store_true',
                        dest='is_paired_end', 
                        default=False,
                        help=("Indicates that reads are paired-end "
                              "(default is single)."))
    
    parser.add_argument("--is_sorted", "-s", action='store_true',
                        dest='is_sorted', 
                        default=False,
                        help=('Indicates that the input BAM file'
                              ' is coordinate-sorted (default '
                              'is False).'))
    
    parser.add_argument("--max_seqs", type=int, default=MAX_SEQS_DEFAULT,
                        help="The maximum number of sequences with different "
                        "allelic combinations to consider remapping "
                        "(default=%d). Read pairs with more allelic "
                        "combinations than MAX_SEQs are discarded" %
                        MAX_SEQS_DEFAULT)

    parser.add_argument("--max_snps", type=int, default=MAX_SNPS_DEFAULT,
                        help="The maximum number of SNPs allowed to overlap "
                        "a read before discarding the read. Allowing higher "
                        "numbers will decrease speed and increase memory "
                        "usage (default=%d)."
                         % MAX_SNPS_DEFAULT)
    
    parser.add_argument("--output_dir", default=None,
                        help="Directory to write output files to. If not "
                        "specified, output files are written to the "
                        "same directory as the input BAM file.")

    parser.add_argument("--snp_dir", action='store', 
                        help="Directory containing SNP text files "
                        "This directory should contain one file per "
                        "chromosome named like chr<#>.snps.txt.gz. "
                        "Each file should contain 3 columns: position "
                        "RefAllele AltAllele. This option should "
                        "only be used if --snp_tab, --snp_index, "
                        "and --haplotype arguments are not used."
                        " If this argument is provided, all possible "
                        "allelic combinations are used (rather "
                        "than set of observed haplotypes).",
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
                        "name per line (file is assumed to be whitespace-"
                        "delimited and first column is assumed to be sample "
                        "name). Sample names should match those present in the "
                        "--haplotype file. Samples are ignored if no haplotype "
                        "file is provided.",
                        metavar="SAMPLES")
                        
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
    
    if options.samples and not options.haplotype:
        # warn because no way to use samples if haplotype file not specified
        sys.stderr.write("WARNING: ignoring --samples argument "
                         "because --haplotype argument not provided")

    return options









def write_read(read, snp_tab, snp_idx, read_pos):
    snp_allele1 = [' '] * read.qlen
    snp_allele2 = [' '] * read.qlen

    for (s_idx, r_idx) in zip(snp_idx, read_pos):
        a1 = snp_tab.snp_allele1[s_idx]
        a2 = snp_tab.snp_allele2[s_idx]

        snp_allele1[r_pos-1] = a1
        snp_allele2[r_pos-1] = a2

    sys.stderr.write("READ: %s\n" % read.query_sequence)
    sys.stderr.write("A1:   %s\n" % "".join(snp_allele1))
    sys.stderr.write("A2:   %s\n" % "".join(snp_allele2))
    
    
        


def count_ref_alt_matches(read, read_stats, snp_tab, snp_idx, read_pos):
    ref_alleles = snp_tab.snp_allele1[snp_idx]
    alt_alleles = snp_tab.snp_allele2[snp_idx]
    
    for i in range(len(snp_idx)):
        if ref_alleles[i] == read.query_sequence[read_pos[i]-1]:
            # read matches reference allele
            read_stats.ref_count += 1
        elif alt_alleles[i] == read.query_sequence[read_pos[i]-1]:
            # read matches non-reference allele
            read_stats.alt_count += 1
        else:
            # read matches neither ref nor other
            read_stats.other_count += 1
            


def get_unique_haplotypes(haplotypes, snp_idx):
    """returns list of vectors of unique haplotypes for this set of SNPs"""
    haps = haplotypes[snp_idx,:].T
    
    # create view of data that joins all elements of column
    # into single void datatype
    h = np.ascontiguousarray(haps).view(np.dtype((np.void, haps.dtype.itemsize * haps.shape[1])))

    # get index of unique columns
    _, idx = np.unique(h, return_index=True)

    return haps[idx,:]
    



            
def generate_haplo_reads(read_seq, snp_idx, read_pos, ref_alleles, alt_alleles,
                         haplo_tab):
    haps = get_unique_haplotypes(haplo_tab, snp_idx)

    # sys.stderr.write("UNIQUE haplotypes: %s\n"
    #                  "read_pos: %s\n"
    #                 % (repr(haps), read_pos))
    
    read_len = len(read_seq)

    new_read_list = []

    # loop over haplotypes
    for hap in haps:
        new_read = []
        cur_pos = 1

        missing_data = False

        # loop over the SNPs to get alleles that make up this haplotype
        for i in range(len(hap)):
            if read_pos[i] > cur_pos:
                # add segment of read
                new_read.append(read_seq[cur_pos-1:read_pos[i]-1])
            # add segment for appropriate allele
            if hap[i] == 0:
                # reference allele
                new_read.append(ref_alleles[i])
            elif hap[i] == 1:
                # alternate allele
                new_read.append(alt_alleles[i])
            else:
                # haplotype has unknown genotype or phasing so skip it...
                # not sure if this is the best thing to do, could instead
                # assume ambiguity of this allele and generate reads with
                # both possible alleles
                missing_data = True
                break
            
            cur_pos = read_pos[i] + 1

        if read_len >= cur_pos:
            # add final segment
            new_read.append(read_seq[cur_pos-1:read_len])
            
        if not missing_data:
            new_seq = "".join(new_read)

            # sanity check: read should be same length as original
            if len(new_seq) != read_len:
                raise ValueError("Expected read len to be %d but "
                                 "got %d.\n"
                                 "ref_alleles: %s\n"
                                 "alt_alleles: %s\n"
                                 "read_pos: %s\n"
                                 "snp_idx: %s\n"
                                 "haps: %s\n" 
                                 % (read_len, len(new_seq),
                                    repr(ref_alleles), repr(alt_alleles),
                                    repr(read_pos), repr(snp_idx),
                                    repr(haps)))

            new_read_list.append("".join(new_seq))

    return new_read_list

    

            
def generate_reads(read_seq, read_pos, ref_alleles, alt_alleles, i):
    """Recursively generate set of reads with all possible combinations
    of alleles (i.e. 2^n combinations where n is the number of snps overlapping
    the reads)
    """
    # TODO: this would use a lot less memory if re-implemented
    # to not use recursion
    
    # create new version of this read with both reference and
    # alternative versions of allele at this index
    idx = read_pos[i]-1
    ref_read = read_seq[:idx] + ref_alleles[i] + read_seq[idx+1:]
    alt_read = read_seq[:idx] + alt_alleles[i] + read_seq[idx+1:]

    if i == len(read_pos)-1:
        # this was the last SNP
        return [ref_read, alt_read]

    # continue recursively with other SNPs overlapping this read
    reads1 = generate_reads(ref_read, read_pos, ref_alleles, alt_alleles, i+1)
    reads2 = generate_reads(alt_read, read_pos, ref_alleles, alt_alleles,  i+1)

    return reads1 + reads2
                


def write_fastq(fastq_file, orig_read, new_seqs):
    n_seq = len(new_seqs)
    i = 1
    for new_seq in new_seqs:
        # Give each read a new name giving:
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
        
        name = "%s.%s.%d.%d" % (orig_read1.qname, pos_str, i, n_pair)
        
        fastq_file1.write("@%s\n%s\n+%s\n%s\n" %
                          (name, pair[0], name, orig_read1.qual))

        rev_seq = util.revcomp(pair[1])
        fastq_file2.write("@%s\n%s\n+%s\n%s\n" %
                          (name, rev_seq, name, orig_read2.qual))

        i += 1
                         


        
    
def filter_reads(files, max_seqs=MAX_SEQS_DEFAULT, max_snps=MAX_SNPS_DEFAULT,
                 samples=None):
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = snptable.SNPTable()
    read_stats = ReadStats()
    read_pair_cache = {}
    cache_size = 0
    read_count = 0
    
    for read in files.input_bam:
        read_count += 1
        # if (read_count % 100000) == 0:
        #     sys.stderr.write("\nread_count: %d\n" % read_count)
        #     sys.stderr.write("cache_size: %d\n" % cache_size)

        # TODO: need to change this to use new pysam API calls
        # but need to check pysam version for backward compatibility
        if read.tid == -1:
            # unmapped read
            read_stats.discard_unmapped += 1
            continue
        
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome
            cur_chrom = files.input_bam.getrname(read.tid)

            if len(read_pair_cache) != 0:
                sys.stderr.write("WARNING: failed to find pairs for %d "
                                 "reads on this chromosome\n" %
                                 len(read_pair_cache))
                read_stats.discard_missing_pair += len(read_pair_cache)
            read_pair_cache = {}
            cache_size = 0
            read_count = 0
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)

            # use HDF5 files if they are provided, otherwise use text
            # files from SNP dir
            if files.snp_tab_h5:
                sys.stderr.write("reading SNPs from file '%s'\n" %
                                 files.snp_tab_h5.filename)
                snp_tab.read_h5(files.snp_tab_h5, files.snp_index_h5,
                                files.hap_h5, cur_chrom, samples)
            else:
                snp_filename = "%s/%s.snps.txt.gz" % (files.snp_dir, cur_chrom)
                sys.stderr.write("reading SNPs from file '%s'\n" % snp_filename)
                snp_tab.read_file(snp_filename)
            
            sys.stderr.write("processing reads\n")

        if read.is_secondary:
            # this is a secondary alignment (i.e. read was aligned more than
            # once and this has align score that <= best score)
            read_stats.discard_secondary += 1
            continue

        if read.is_paired:
            if read.mate_is_unmapped:
                # other side of pair not mapped
                # we could process as single... but these not likely
                # useful so discard
                # process_single_read(read, read_stats, files,
                #                     snp_tab, max_seqs, max_snps)
                read_stats.discard_mate_unmapped += 1
            elif(read.next_reference_name == cur_chrom or
                 read.next_reference_name == "="):
                # other pair mapped to same chrom

                # sys.stderr.write("flag: %s" % read.flag)
                if not read.is_proper_pair:
                    # sys.stderr.write(' => improper\n')
                    read_stats.discard_improper_pair += 1
                    continue
                # sys.stderr.write(' => proper\n')

                if read.qname in read_pair_cache:
                    # we already saw prev pair, retrieve from cache
                    read1 = read_pair_cache[read.qname]
                    read2 = read
                    del read_pair_cache[read.qname]
                    cache_size -= 1

                    if read2.next_reference_start != read1.reference_start:
                        sys.stderr.write("WARNING: read pair positions "
                                         "do not match for pair %s\n" %
                                         read.qname)
                    else:
                        process_paired_read(read1, read2, read_stats,
                                            files, snp_tab, max_seqs,
                                            max_snps)
                else:
                    # we need to wait for next pair
                    read_pair_cache[read.qname] = read

                    cache_size += 1

                    
            else:
                # other side of pair mapped to different
                # chromosome, discard this read
                read_stats.discard_different_chromosome += 1

        else:
            process_single_read(read, read_stats, files, snp_tab,
                                max_seqs, max_snps)

    if len(read_pair_cache) != 0:
        sys.stderr.write("WARNING: failed to find pairs for %d "
                         "reads on this chromosome\n" %
                         len(read_pair_cache))
        read_stats.discard_missing_pair += len(read_pair_cache)
    
    read_stats.write(sys.stderr)
                     

def process_paired_read(read1, read2, read_stats, files,
                        snp_tab, max_seqs, max_snps):
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

            # limit recursion here by discarding reads that
            # overlap too many SNPs
            if len(snp_read_pos) > max_snps:
                read_stats.discard_excess_snps += 1
                return

            if files.hap_h5:
                # generate reads using observed set of haplotypes
                read_seqs = generate_haplo_reads(read.query_sequence,
                                                 snp_idx,
                                                 snp_read_pos,
                                                 ref_alleles, alt_alleles,
                                                 snp_tab.haplotypes)
            else:
                # generate all possible allelic combinations of reads
                read_seqs = generate_reads(read.query_sequence, snp_read_pos,
                                           ref_alleles, alt_alleles, 0)
            
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
        new_reads[0].append(read1.query_sequence)
        new_reads[1].append(read2.query_sequence)

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
        orig_pair = (read1.query_sequence, read2.query_sequence)
                                 
        if orig_pair in unique_pairs:
            unique_pairs.remove(orig_pair)
            
        # write read pair to fastqs for remapping
        write_pair_fastq(files.fastq1, files.fastq2, read1, read2,
                         unique_pairs)

        # Write read to 'remap' BAM for consistency with previous
        # implementation of script. Probably not needed and will result in
        # BAM that is not coordinate sorted. Possibly remove this...
        files.remap_bam.write(read1)
        files.remap_bam.write(read2)
        read_stats.remap_pair += 1
        

        
    

    

def process_single_read(read, read_stats, files, snp_tab, max_seqs,
                        max_snps):
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

        # limit recursion here by discarding reads that
        # overlap too many SNPs
        if len(snp_read_pos) > max_snps:
            read_stats.discard_excess_snps += 1
            return

        if files.hap_h5:
            read_seqs = generate_haplo_reads(read.query_sequence, snp_idx,
                                             snp_read_pos,
                                             ref_alleles, alt_alleles,
                                             snp_tab.haplotypes)
        else:
            read_seqs = generate_reads(read.query_sequence,  snp_read_pos,
                                       ref_alleles, alt_alleles, 0)

        # make set of unique reads, we don't want to remap
        # duplicates, or the read that matches original
        unique_reads = set(read_seqs)
        if read.query_sequence in unique_reads:
            unique_reads.remove(read.query_sequence)
        
        if len(unique_reads) == 0:
            # only read generated matches original read,
            # so keep original
            files.keep_bam.write(read)
            read_stats.keep_single += 1
        elif len(unique_reads) < max_seqs:
            # write read to fastq file for remapping
            write_fastq(files.fastq_single, read, unique_reads)

            # write read to 'to remap' BAM
            # this is probably not necessary with new implmentation
            # but kept for consistency with previous version of script
            files.remap_bam.write(read)
            read_stats.remap_single += 1
        else:
            # discard read
            read_stats.discard_excess_reads += 1
            return

    else:
        # no SNPs overlap read, write to keep file
        files.keep_bam.write(read)
        read_stats.keep_single += 1
            



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


        
def main(bam_filenames, is_paired_end=False,
         is_sorted=False, max_seqs=MAX_SEQS_DEFAULT,
         max_snps=MAX_SNPS_DEFAULT, output_dir=None,
         snp_dir=None, snp_tab_filename=None,
         snp_index_filename=None,
         haplotype_filename=None, samples=None):

    files = DataFiles(bam_filenames,  is_sorted, is_paired_end,
                      output_dir=output_dir,
                      snp_dir=snp_dir,
                      snp_tab_filename=snp_tab_filename,
                      snp_index_filename=snp_index_filename,
                      haplotype_filename=haplotype_filename)
    
    filter_reads(files, max_seqs=max_seqs, max_snps=max_snps,
                 samples=samples)

    files.close()
    
    

if __name__ == '__main__':
    sys.stderr.write("command line: %s\n" % " ".join(sys.argv))
    sys.stderr.write("python version: %s\n" % sys.version)
    sys.stderr.write("pysam version: %s\n" % pysam.__version__)

    util.check_pysam_version()
        
    options = parse_options()
    samples = parse_samples(options.samples)
    
    main(options.bam_filename,
         is_paired_end=options.is_paired_end, is_sorted=options.is_sorted,
         max_seqs=options.max_seqs, max_snps=options.max_snps,
         output_dir=options.output_dir,
         snp_dir=options.snp_dir,
         snp_tab_filename=options.snp_tab,
         snp_index_filename=options.snp_index,
         haplotype_filename=options.haplotype,
         samples=samples)
         
    
