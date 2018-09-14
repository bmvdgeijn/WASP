import sys
import os
import gzip
import argparse
import numpy as np
from itertools import product, groupby

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
            self.snp_tab_h5 = tables.open_file(snp_tab_filename, "r")
            self.snp_index_h5 = tables.open_file(snp_index_filename, "r")
            self.hap_h5 = tables.open_file(haplotype_filename, "r")
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
        
        # create output dir if does not exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

            
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
            self.fastq1 = gzip.open(self.fastq1_filename, "wt")
            self.fastq2 = gzip.open(self.fastq2_filename, "wt")
            self.fastq_single_filename = self.prefix + ".remap.single.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wt")
            sys.stderr.write("  %s\n  %s\n  %s\n" %
                             (self.fastq1_filename,
                              self.fastq2_filename,
                              self.fastq_single_filename))
            
        else:
            self.fastq_single_filename = self.prefix + ".remap.fq.gz"
            self.fastq_single = gzip.open(self.fastq_single_filename, "wt")
            sys.stderr.write("  %s\n" % (self.fastq_single_filename))

        self.input_bam = pysam.Samfile(self.bam_sort_filename, "r")
        self.keep_bam = pysam.Samfile(self.keep_filename, "w",
                                      template=self.input_bam)
        self.remap_bam = pysam.Samfile(self.remap_filename, "w",
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

        # number of chimeric reads discarded
        self.discard_supplementary = 0

        # number of reads discarded because of too many overlapping SNPs
        self.discard_excess_snps = 0
        
        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0

        # when read pairs share SNP locations but have different alleles there
        self.discard_discordant_shared_snp = 0
        
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
                         "  supplementary alignment: %d\n"
                         "  excess overlapping snps: %d\n"
                         "  excess allelic combinations: %d\n"
                         "  read pairs with discordant shared SNPs: %d\n"
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
                          self.discard_supplementary,
                          self.discard_excess_snps,
                          self.discard_excess_reads,
                          self.discard_discordant_shared_snp,
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
        ref = ref_alleles[i].decode("utf-8")
        alt = alt_alleles[i].decode("utf-8")
        
        if ref == read.query_sequence[read_pos[i]-1]:
            # read matches reference allele
            read_stats.ref_count += 1
        elif alt == read.query_sequence[read_pos[i]-1]:
            # read matches non-reference allele
            read_stats.alt_count += 1
        else:
            # read matches neither ref nor other
            read_stats.other_count += 1
            

def get_unique_haplotypes(haplotypes, phasing, snp_idx):
    """
    returns list of vectors of unique haplotypes for this set of SNPs
    all possible combinations of ref/alt are calculated at unphased sites
    """
    haps = haplotypes[snp_idx,:].T
    
    # get phasing of SNPs for each individual as bool array
    # True = unphased and False = phased
    if phasing is not None:
        phasing = np.logical_not(phasing[snp_idx, :].T.astype(bool))
    else:
        # assume all SNPs are unphased
        phasing = np.full((int(haps.shape[0]/2), haps.shape[1]), True)

    # if a haplotype has unphased SNPs, generate all possible allelic
    # combinations and add each combination as a new haplotype
    new_haps = []
    # iterate through each individual
    for i in range(len(phasing)):
        # get the haplotype data for this individual
        hap_pair = haps[i*2:(i*2)+2]
        # get bool index of cols in hap_pair that contain hets
        hets = np.not_equal(hap_pair[0], hap_pair[1])
        # get bool index of cols in hap_pair that are both unphased and hets
        phase = np.logical_and(phasing[i], hets)
        # get all combinations of indices at unphased cols
        # then, index into hap_pair with each combination
        for j in product([0, 1], repeat=sum(phase)):
            # index into hap_pair using all ref genes at genotyped columns
            ref = np.repeat(0, hap_pair.shape[1])
            np.place(ref, phase, j)
            new_haps.append(hap_pair[ref, range(len(ref))])
            # index into hap_pair using all alt genes at genotyped columns
            alt = np.repeat(1, hap_pair.shape[1])
            np.place(alt, phase, j)
            new_haps.append(hap_pair[alt, range(len(alt))])
    # add new haps to old haps
    haps = np.concatenate((haps, new_haps))


    # create view of data that joins all elements of column
    # into single void datatype
    h = np.ascontiguousarray(haps).view(np.dtype((np.void, haps.dtype.itemsize * haps.shape[1])))

    # get index of unique columns
    _, idx = np.unique(h, return_index=True)


    return haps[idx,:]


            
def generate_haplo_reads(read_seq, snp_idx, read_pos, ref_alleles, alt_alleles,
                         haplo_tab, phase_tab):
    """
      read_seq - a string representing the the sequence of the read in question
      snp_index - a list of indices of SNPs that this read overlaps
      read_pos - a list of positions in read_seq that overlap SNPs
      ref_alleles - a np array of reference alleles with
                    indices corresponding to snp_index
      alt_alleles - a np array of alternate alleles with
                    indices corresponding to snp_index
      haplo_tab - a pytables node with haplotypes from haplotype.h5
    """
    haps = get_unique_haplotypes(haplo_tab, phase_tab, snp_idx)

    # sys.stderr.write("UNIQUE haplotypes: %s\n"
    #                  "read_pos: %s\n"
    #                 % (repr(haps), read_pos))
    
    read_len = len(read_seq)

    new_read_list = set()

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
                new_read.append(ref_alleles[i].decode("utf-8"))
            elif hap[i] == 1:
                # alternate allele
                new_read.append(alt_alleles[i].decode("utf-8"))
            else:
                # haplotype has unknown genotype so skip it...
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

            new_read_list.add("".join(new_seq))

    return new_read_list

    

            
def generate_reads(read_seq, read_pos, ref_alleles, alt_alleles):
    """Generate set of reads with all possible combinations
    of alleles (i.e. 2^n combinations where n is the number of snps overlapping
    the reads)
    """
    reads = [read_seq]
    # iterate through all snp locations
    for i in range(len(read_pos)):
        idx = read_pos[i]-1
        # for each read we've already created...
        for j in range(len(reads)):
            read = reads[j]
            # create a new version of this read with both reference
            # and alternative versions of the allele at this index
            reads.append(
              read[:idx] + ref_alleles[i].decode("utf-8") + read[idx+1:]
            )
            reads.append(
              read[:idx] + alt_alleles[i].decode("utf-8") + read[idx+1:]
            )
    return set(reads)


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

        if read.is_supplementary:
            # this is a supplementary alignment (ie chimeric and not the representative alignment)
            read_stats.discard_supplementary += 1
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


def slice_read(read, indices):
    """slice a read by an array of indices"""
    return "".join(np.array(list(read))[indices])


def group_reads_by_snps(reads, snp_read_pos):
    """
    group the reads by strings containing the combinations of ref/alt alleles
    among the reads at the shared_snps. return a list of sets of reads - one
    for each group
    """
    # group the reads by the snp string and create a list to hold the groups
    return [
        set(reads) for hap, reads in
        groupby(
          # note that groupby needs the data to be sorted by the same key func
          sorted(reads, key=lambda read: slice_read(read, snp_read_pos)),
          key=lambda read: slice_read(read, snp_read_pos)
        )
    ]


def read_pair_combos(old_reads, new_reads, max_seqs, snp_idx, snp_read_pos):
    """
    Collects all unique combinations of read pairs. Handles the possibility of
    shared SNPs among the pairs (ie doesn't treat them as independent).
    Returns False before more than max_seqs pairs are created or None
    when the original read pair has discordant alleles at shared SNPs.
    Input:
        old_reads - a tuple of length 2, containing the pair of original reads
        new_reads - a list of two sets, each containing the reads generated
                    from old_reads for remapping
        snp_index - a list of two lists of the indices of SNPs that overlap
                    with old_reads
        snp_read_pos - a list of two lists of the positions in old_reads where
                       SNPs are located
    Output:
        unique_pairs - a set of tuples, each representing a unique pair of
                       new_reads
    """
    # get the indices of the shared SNPs in old_reads
    for i in range(len(snp_read_pos)):
        # get the indices of the SNP indices that are in both reads
        idx_idxs = np.nonzero(np.in1d(snp_idx[i], snp_idx[(i+1) % 2]))[0]
        # now, use the indices in idx_idxs to get the relevant snp positions
        # and convert positions to indices
        snp_read_pos[i] = np.array(snp_read_pos[i], dtype=int)[idx_idxs] - 1
    # check: are there discordant alleles at the shared SNPs?
    # if so, discard these reads
    if (
        slice_read(old_reads[0], snp_read_pos[0])
        != slice_read(old_reads[1], snp_read_pos[1])
    ):
        return None
    # group reads by the alleles they have at shared SNPs
    for i in range(len(new_reads)):
        new_reads[i] = group_reads_by_snps(
            new_reads[i], snp_read_pos[i]
        )
    unique_pairs = set()
    # calculate unique combinations of read pairs only among reads that
    # have the same alleles at shared SNPs (ie if they're in the correct group)
    for group in range(len(new_reads[0])):
        for pair in product(new_reads[0][group], new_reads[1][group]):
            if len(unique_pairs) <= max_seqs:
                unique_pairs.add(pair)
            else:
                return False
    return unique_pairs


def process_paired_read(read1, read2, read_stats, files,
                        snp_tab, max_seqs, max_snps):
    """Checks if either end of read pair overlaps SNPs or indels
    and writes read pair (or generated read pairs) to appropriate
    output files"""

    new_reads = []
    pair_snp_idx = []
    pair_snp_read_pos = []

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
                                                 snp_tab.haplotypes,
                                                 snp_tab.phase)
            else:
                # generate all possible allelic combinations of reads
                read_seqs = generate_reads(read.query_sequence, snp_read_pos,
                                           ref_alleles, alt_alleles)
            
            new_reads.append(read_seqs)
            pair_snp_idx.append(snp_idx)
            pair_snp_read_pos.append(snp_read_pos)
        else:
            # no SNPs or indels overlap this read
            new_reads.append(set())
            pair_snp_idx.append([])
            pair_snp_read_pos.append([])

    if len(new_reads[0]) == 0 and len(new_reads[1]) == 0:
        # neither read overlapped SNPs or indels
        files.keep_bam.write(read1)
        files.keep_bam.write(read2)
        read_stats.keep_pair += 1
    else:
        # add original version of both sides of pair
        new_reads[0].add(read1.query_sequence)
        new_reads[1].add(read2.query_sequence)

        if len(new_reads[0]) + len(new_reads[1]) > max_seqs:
            # quit now before generating a lot of read pairs
            read_stats.discard_excess_reads += 2
            return

        # get all unique combinations of read pairs
        unique_pairs = read_pair_combos(
            (read1.query_sequence, read2.query_sequence), new_reads,
            max_seqs, pair_snp_idx, pair_snp_read_pos
        )
        # if unique_pairs is None or False we should discard these reads
        if unique_pairs is None:
            read_stats.discard_discordant_shared_snp += 1
            return
        elif not unique_pairs:
            read_stats.discard_excess_reads += 2
            return

        # remove original read pair, if present
        unique_pairs.discard((read1.query_sequence, read2.query_sequence))
            
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
                                             snp_tab.haplotypes,
                                             snp_tab.phase)
        else:
            read_seqs = generate_reads(read.query_sequence,  snp_read_pos,
                                       ref_alleles, alt_alleles)

        # we don't want the read that matches the original
        read_seqs.discard(read.query_sequence)
        
        if len(read_seqs) == 0:
            # only read generated matches original read,
            # so keep original
            files.keep_bam.write(read)
            read_stats.keep_single += 1
        elif len(read_seqs) < max_seqs:
            # write read to fastq file for remapping
            write_fastq(files.fastq_single, read, read_seqs)

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
            f = gzip.open(samples_str, "rt")
        else:
            f = open(samples_str, "rt")

        for line in f:
            # assume first token in line is sample name
            samples.append(line.split()[0])

        sys.stderr.write("read %d sample names from file '%s'\n" %
                         (len(samples), samples_str))
                    
        f.close()
    else:    
        # otherwise assume comma-delimited string
        if ("/" in samples_str or "\\" in samples_str):
            sys.stderr.write("WARNING: --samples argument (%s) "
                             "does not look like list of sample names "
                             "(contains '/' or '\\') but is not path to "
                             "valid file. Assuming it is list of sample "
                             "names anyway." % samples_str)

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
    sys.stderr.write("pytables version: %s\n" % tables.__version__)

    util.check_pysam_version()
    util.check_pytables_version()
        
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
         
    
