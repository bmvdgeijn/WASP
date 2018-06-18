

import argparse
import gzip
import numpy as np
import numpy.random
import sys
import string

import tables



class Haplotypes(object):
    def __init__(self, pos_array, ref_allele_array, alt_allele_array,
                 hap1_array, hap2_array):
        self.pos = pos_array
        self.ref_allele = ref_allele_array
        self.alt_allele = alt_allele_array
        self.hap1 = hap1_array
        self.hap2 = hap2_array

        

class ReadCoord(object):
    def __init__(self, chrom_name, left_start, left_end, right_start, right_end):

        if left_start > left_end:
            raise ValueError("left_start must be <= left_end")
        if left_start > right_start:
            raise ValueError("left_start must be <= right_start")
        if right_start > right_end:
            raise ValueError("right_start must be <= right_end")

        self.chrom_name = chrom_name
        self.left_start = left_start
        self.left_end = left_end
        self.right_start = right_start
        self.right_end = right_end
        
        

dna_comp = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global dna_comp

    if dna_comp is None:
        dna_comp = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(dna_comp)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]

    

def parse_options():
    parser=argparse.ArgumentParser(
        description="Simulates paired-end reads that can be used to "
        "test the mapping pipeline")

    parser.add_argument("--seq", required=True,
                        help="Path to HDF5 file containing "
                        "genome sequence. (Can be created "
                        "using fasta2h5 program)",
                        metavar="SEQ_H5_FILE")
    
    parser.add_argument("--n_reads", action='store',
                        required=True, type=int,
                        help="number of reads to simulate")

    parser.add_argument("--out_fastq1", action='store',
                        required=True,
                        help="output file to write read1 to")

    parser.add_argument("--out_fastq2", action='store',
                        required=True,
                        help="output file to write read2 to")

    parser.add_argument("--read_len", default=36, 
                        type=int)

    parser.add_argument("--insert_size_mean", default=100.0,
                        help="mean insert size (assumed to be normally distributed)",
                        type=float)
    
    parser.add_argument("--insert_size_sd", default=50.0,
                        help="standard devaition of insert size",
                        type=float)

    parser.add_argument("--chrom", default="chr22",
                        help="for now just simulate from one chrom"
                        "in future may simulate uniformly across entire genome...")

    parser.add_argument("--hap_file", required=True,
                        help="path to file containing haplotypes and alleles "
                        "The file should contain 5 columns:\n"
                        "  position RefAllele AltAllele hap1 hap2\n"
                        "  example: 16050984 C G 1 0")

    
    return parser.parse_args()





def read_haps(hap_file):
    if hap_file.endswith(".gz"):
        f = gzip.open(hap_file)
    else:
        f = open(hap_file)

    pos_list = []
    ref_allele_list = []
    alt_allele_list = []
    hap1_list = []
    hap2_list = []
    
    for line in f:
        words = line.rstrip().split()
        ref_allele = words[1].upper()
        alt_allele = words[2].upper()

        # ignore indels for now
        if ref_allele not in ("A", "C", "T", "G"):
            continue
        if alt_allele not in ("A", "C", "T", "G"):
            continue

        pos_list.append(int(words[0]))
        
        # append ascii code instead of char, as this is how sequence
        # is represented in HDF5 files
        ref_allele_list.append(ord(ref_allele))
        alt_allele_list.append(ord(alt_allele))

        hap1_list.append(int(words[3]))
        hap2_list.append(int(words[4]))

    pos_array = np.array(pos_list, dtype=np.int32)
    ref_allele_array = np.array(ref_allele_list, dtype=np.uint8)
    alt_allele_array = np.array(alt_allele_list, dtype=np.uint8)
    hap1_array = np.array(hap1_list, dtype=np.uint8)
    hap2_array = np.array(hap2_list, dtype=np.uint8)

    haplotypes = Haplotypes(pos_array, ref_allele_array, alt_allele_array,
                            hap1_array, hap2_array)
    
    f.close()

    return haplotypes


    

def write_reads(file1, file2, read_coord, left_seq, right_seq):
    """Outputs PE reads to two separate files in fastq format"""

    if len(right_seq) != len(left_seq):
        raise ValueError("length of right and left seqs does not match")
    
    read_len = len(left_seq)
        
    qual_str = "h" * read_len

    tile = np.random.randint(32767)
    
    # use chromosome number as lane number
    # random identifier as tile
    # start, end of fragment as x, y pixels
    base_id = "PE%d:%s:%d:%d:%d#0" % (read_len, read_coord.chrom_name, 
                                      tile, read_coord.left_start, 
                                      read_coord.right_end)

    # id1 = base_id + "/1"
    # id2 = base_id + "/2"
    id1 = id2 = base_id

    # Make left read "read1" 50 % of time
    if np.random.randint(2) == 0:
        # make left read read1
        file1.write("@%s\n%s\n+\n%s\n" % (id1, left_seq, qual_str))
        file2.write("@%s\n%s\n+\n%s\n" % (id2, right_seq, qual_str))
    else:
        file1.write("@%s\n%s\n+\n%s\n" % (id1, right_seq, qual_str))
        file2.write("@%s\n%s\n+\n%s\n" % (id2, left_seq, qual_str))
    




def gen_read_coords(options, haps, het_only=True):
    """generate coordinates for a read pair that overlaps a SNP"""
    # to make this more efficient, could observer random
    # vars for many reads at once, rather than single read at a time

    if het_only:
        # select a heterozygous site to overlap
        is_het = haps.hap1 != haps.hap2
        i = numpy.random.randint(np.sum(is_het))
        snp_pos = haps.pos[is_het][i]

        
        sys.stderr.write("selected SNP %d %s/%s\n" % 
                         (snp_pos, chr(haps.ref_allele[i]), chr(haps.alt_allele[i])))

        
    else:
        # select any SNP to overlap
        i = numpy.random.randint(haps.pos.shape[0])
        snp_pos = haps.pos[i]


    # at what read position should het site be?
    snp_read_pos = numpy.random.randint(options.read_len)

    # what should insert size be?
    insert_size = numpy.random.normal(options.insert_size_mean,
                                      options.insert_size_sd)
    insert_size = int(np.rint(insert_size))

    # insert size cannot be smaller than read size...
    insert_size = max(options.read_len, insert_size)
    
    # does left or right read overlap SNP?
    if numpy.random.randint(2) == 0:
        # left read overlaps SNP
        left_start = snp_pos - snp_read_pos + 1
        left_end = left_start + options.read_len - 1
        right_end = left_start + insert_size - 1
        right_start = right_end - options.read_len + 1
    else:
        # right read overlaps SNP
        right_start = snp_pos - snp_read_pos + 1
        right_end = right_start + options.read_len - 1
        left_start = right_end - insert_size + 1
        left_end = left_start + options.read_len - 1

    read_coord = ReadCoord(options.chrom, left_start, left_end, right_start, right_end)

    return read_coord



def gen_seqs(read_coord, hap1_seq, hap2_seq):
    """makes the sequence strings for each read, 
    choosing 1 haplotype at random to obtain sequence from"""
    # randomly select haplotype
    if np.random.randint(2) == 0:
        chrom_seq = hap1_seq
    else:
        chrom_seq = hap2_seq

    # sys.stderr.write("%s-%s %s-%s\n" % (str(read_coord.left_start), 
    #                                     str(read_coord.left_end),
    #                                     str(read_coord.right_start),
    #                                     str(read_coord.right_end)))
    
    s = read_coord.left_start - 1
    e = read_coord.left_end
    left_read_seq = chrom_seq[s:e]
    # sys.stderr.write("%s\n%s\n" % (hap1_seq[s:e], hap2_seq[s:e]))
    
    s = read_coord.right_start - 1
    e = read_coord.right_end
    right_read_seq = revcomp(chrom_seq[s:e])
    # sys.stderr.write("%s\n%s\n\n" % (hap1_seq[s:e], hap2_seq[s:e]))

    return right_read_seq, left_read_seq

    
    


def make_hap_seqs(haps, options):
    """Makes a chromosome sequence for each haplotype"""
    seq_h5 = tables.open_file(options.seq)
    
    node_name = "/%s" % options.chrom
    if node_name not in seq_h5:
        raise ValueError("chromosome %s is not in sequence h5 file" % options.chrom)
        
    seq_node = seq_h5.get_node("/%s" % options.chrom)

    seq_array1 = seq_node[:]
    seq_array2 = np.array(seq_node[:])

    is_alt = (haps.hap1 == 1)
    seq_array1[haps.pos[is_alt] - 1] = haps.alt_allele[is_alt]

    is_alt = (haps.hap2 == 1)    
    seq_array2[haps.pos[is_alt] - 1] = haps.alt_allele[is_alt]
    

    seq1 = "".join([chr(x) for x in seq_array1])
    seq2 = "".join([chr(x) for x in seq_array2])

    seq_h5.close()

    return seq1, seq2


    


def main():
    options = parse_options()

    sys.stderr.write("reading haplotype information\n")
    haplotypes = read_haps(options.hap_file)
    
    sys.stderr.write("making haplotype sequences\n")
    hap1_seq, hap2_seq = make_hap_seqs(haplotypes, options)

    if options.out_fastq1.endswith(".gz"):
        fastq1_file = gzip.open(options.out_fastq1, "w")
    else:
        fastq1_file = open(options.out_fastq1, "w")

    if options.out_fastq2.endswith(".gz"):
        fastq2_file = gzip.open(options.out_fastq2, "w")
    else:
        fastq2_file = open(options.out_fastq2, "w")
    
    
    sys.stderr.write("simulating reads\n")
    for i in range(options.n_reads):
        # generate coords for a read pair, where one read pair overlaps a SNP
        read_coord = gen_read_coords(options, haplotypes)

        # make read pair sequences 
        left_seq, right_seq = gen_seqs(read_coord, hap1_seq, hap2_seq)

        write_reads(fastq1_file, fastq2_file, read_coord, left_seq, right_seq)

    fastq1_file.close()
    fastq2_file.close()


    


if __name__ == "__main__":
    main()
