import sys
import os
import subprocess
import gzip
import argparse
import numpy as np

import pysam


NUCLEOTIDES = set(['A', 'C', 'T', 'G'])
SNP_UNDEF = -1

MAX_SEQS_DEFAULT = 10

# codes for CIGAR string
BAM_CMATCH     = 0   # M - match/mismatch to ref M
BAM_CINS       = 1   # I - insertion in read relative to ref
BAM_CDEL       = 2   # D - deletion in read relative to ref
BAM_CREF_SKIP  = 3   # N - skipped region from reference (e.g. intron)
BAM_CSOFT_CLIP = 4   # S - soft clipping (clipped sequence present in seq)
BAM_CHARD_CLIP = 5   # H - hard clipping (clipped sequence NOT present in seq)
BAM_CPAD       = 6   # P - padding (silent deletion from padded reference)
BAM_CEQUAL     = 7   # = - sequence match
BAM_CDIFF      = 8   # X - sequence mismatch

class SNPTable(object):
    def __init__(self):
        self.clear()

    def clear(self):
        # snp_index and indel_index are arrays of length
        # max(snp_pos, indel_pos) that provide lookup
        # into snp_pos, snp_allele1, etc. by chromosome position.
        # For example, if the first and second snps on the chromosome are
        # at positions 1234, 1455 then elements 1233 and 1444 of the
        # snp_index array will be 0 and 1 (and can be used to lookup
        # info for the SNP in snp_pos, snp_allele1, snp_allele2 arrays)
        self.snp_index = np.array([], dtype=np.int32)
        self.snp_pos = np.array([], dtype=np.int32)
        self.snp_allele1 = np.array([], dtype="|S1")
        self.snp_allele2 = np.array([], dtype="|S1")
        self.indel_index = np.array([], dtype=np.int32)
        self.indel_pos = np.array([], dtype=np.int32)
        self.indel_allele1 = []
        self.indel_allele2 = []
        
    
    def read_file(self, filename):
        """read in SNPs and indels from text input file"""
        sys.stderr.write("reading SNPs from file '%s'\n" % filename)

        try:
            if filename.endswith(".gz"):
                f = gzip.open(filename)
            else:
                f = open(filename, "r")
        except IOError:
            sys.stderr.write("WARNING: unable to read from file '%s', "
                             "assuming no SNPs for this chromosome\n" %
                             filename)
            self.clear()
            return
        
        snp_pos_list = []
        snp_allele1_list = []
        snp_allele2_list = []
        max_pos = 0
        indel_pos_list = []
        indel_allele1_list = []
        indel_allele2_list = []

        for line in f:
            words = line.split()

            if(len(words) != 3):
                raise ValueError("expected 3 values per SNP file line got %d:\n"
                                 "%s\n" % (line, len(words)))

            pos = int(words[0])
            a1 = words[1].upper()
            a2 = words[2].upper()

            if pos <= 0:
                raise ValueError("expected SNP position to be >= 1:\n%s\n" %
                                 line)

            if pos > max_pos:
                max_pos = pos

            if (len(a1) == 1) and (len(a2) == 1):
                if a1 in NUCLEOTIDES and a2 in NUCLEOTIDES:
                    # this is a SNP
                    snp_pos_list.append(pos)
                    snp_allele1_list.append(a1)
                    snp_allele2_list.append(a2)
                else:
                    if ("-" in a1) or ("-" in a2):
                        # 1bp indel
                        # reads overlapping indels are thrown out, although
                        # we will likely handle indels better soon
                        indel_pos_list.append(pos)
                        indel_allele1_list.append(a1)
                        indel_allele2_list.append(a2)
                    else:
                        sys.stderr.write("WARNING: unexpected character "
                                         "in SNP:\n%s\n" % line)

        f.close()

        # convert lists to numpy arrays, which allow for faster
        # lookups and use less memory
        self.snp_pos = np.array(snp_pos_list, dtype=np.int32)
        del snp_pos_list
        self.snp_allele1 = np.array(snp_allele1_list, dtype="|S1")
        del snp_allele1_list
        self.snp_allele2 = np.array(snp_allele2_list, dtype="|S1")
        del snp_allele2_list

        self.indel_pos = np.array(indel_pos_list, dtype=np.int32)
        del indel_pos_list

        # make another array that makes it easy to lookup SNPs by their position
        # on the chromosome
        self.snp_index = np.empty(max_pos, dtype=np.int32)
        self.snp_index[:] = SNP_UNDEF
        self.snp_index[self.snp_pos-1] = np.arange(self.snp_pos.shape[0])

        self.indel_index = np.empty(max_pos, dtype=np.int32)
        self.indel_index[:] = SNP_UNDEF
        self.indel_index[self.indel_pos-1] = np.arange(self.indel_pos.shape[0])
        self.indel_allele1 = indel_allele1_list
        self.indel_allele2 = indel_allele2_list

    
    def get_overlapping_snps(self, read):
        """returns list containing indices of overlapping SNPs, as well 
        as index in read sequence of the overlapping SNPs"""
        
        # read.cigar is a list of tuples. Each tuple has two entries. The first
        # entry specifies the character in the cigar and the second entry
        # specifies the length of that character. The values are
        # M       BAM_CMATCH      0
        # I       BAM_CINS        1
        # D       BAM_CDEL        2
        # N       BAM_CREF_SKIP   3
        # S       BAM_CSOFT_CLIP  4
        # H       BAM_CHARD_CLIP  5
        # P       BAM_CPAD        6
        # =       BAM_CEQUAL      7
        # X       BAM_CDIFF       8
        # E.g. (0, 5) means 5 matches, and (4, 2) means a soft clip of 2bp

        read_start = 0
        read_end = 0
        genome_start = read.pos + 1
        genome_end = read.pos + 1
        
        snp_idx_list = []
        read_idx_list = []
        
        for cigar in read.cigar:
            op = cigar[0] # CIGAR 'operation'
            op_len  = cigar[1] # length of operation
            
            if (op == BAM_CMATCH) or (op == BAM_CEQUAL) or (op == BAM_CDIFF):
                # match or mismatch to reference
                read_start = read_end + 1
                read_end = read_start + op_len - 1
                genome_start = genome_start + 1
                genome_end = genome_start + op_len - 1

                # TODO: check for SNP in this genome segment
                s = genome_start - 1
                e = min(genome_end, self.snp_index.shape[0])
                s_idx = self.snp_index[s:e]
                offsets = np.where(s_idx != SNP_UNDEF)[0]
                
                if offsets.shape[0] > 0:
                    # there are overlapping SNPs
                    snp_idx_list.extend(s_idx[offsets])
                    # get the offset of the SNPs into the read
                    read_idx = offsets + read_start - 1
                    read_idx_list.extend(read_idx)

            elif op == BAM_CINS:
                # insert in read relative to reference
                read_start = read_end + 1
                read_end = read_start + op_len - 1

                # genome sequence does not advance, no possibility
                # for read to overlap SNP, since these bases do
                # not exist in reference

                # TODO: handle indels

            elif op == BAM_CDEL:
                # deletion in read relative to reference
                genome_start = genome_start + 1
                genome_end   = genome_start + op_len - 1

                # read sequence does not advance, no possibility
                # for read to overlap SNP, since these bases do
                # not exist in read

                # TODO handle indels

            elif op == BAM_CREF_SKIP:
                # section of skipped reference, such as intron
                genome_end = genome_end + op_len
                genome_start = genome_end

                # do nothing with SNPs/indels in this region
                # since they are skipped
                
            elif op == BAM_CSOFT_CLIP:
                # this part of read skipped
                read_start = read_end + 1
                read_end = read_start + op_len - 1

                # This is like insert, but at the beginning of the read.

                # TODO: handle indels?

            elif seq_type == BAM_CHARD_CLIP:
                # these bases not included in read or genome
                pass

            elif seq_type == BAM_CPAD:
                # like an insert, likely only used in multiple-sequence
                # alignment where inserts may be of different lengths
                # in different seqs
                read_start += read_end + 1
                read_end = read_start + op_len - 1

                # TODO: handle indels?

            else:
                raise ValueError("unknown CIGAR code %d" % op)


        if read_end != read.qlen:
            raise ValueError("length of read segments in CIGAR %d "
                             "does not add up to query length (%d)" %
                             (read_end, read.qlen))
        
        return snp_idx_list, read_idx_list
            

                
                
                
                


        
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
                
    

    
def filter_reads(files, max_seqs=MAX_SEQS_DEFAULT, max_snps=5):
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = SNPTable()

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

            if len(snp_idx) > 1:
                write_read(read, snp_tab, snp_idx, read_idx)
                for r in unique_reads:
                    sys.stderr.write("REMAP:%s\n" % r)
            
            # TODO: write read to file for remapping
        else:
            # no SNPs overlap read, write to keep file
            files.keep_bam.write(read)

    sys.stderr.write("read SNP ref matches: %d\n" % n_ref_match)
    sys.stderr.write("read SNP alt matches: %d\n" % n_alt_match)
    sys.stderr.write("read SNP mismatches: %d\n" % n_other)

    mismatch_pct = 100.0 * float(n_other) / (n_other + n_ref_match +
                                             n_alt_match)

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
    
