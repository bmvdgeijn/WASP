

import snptable
import gzip
import os

import numpy as np

import pysam




class Data():

    def __init__(self, data_dir="test_data",
                 snp_filename="test_data/snp_tab.txt.gz",
                 sam_filename="test_data/test.sam"):
        self.data_dir = data_dir
        self.snp_filename = snp_filename
        self.sam_filename = sam_filename

        self.snp_list = [(10, "A", "C"),
                         (20, "T", "G"),
                         (100, "A", "T")]
        
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)


    def setup(self):
        snp_file = gzip.open(self.snp_filename, "wt")

        for snp in self.snp_list:
            snp_file.write("%d %s %s\n" % (snp[0], snp[1], snp[2]))
        snp_file.close()

        
    
    def write_sam_header(self, sam_file):
        sam_file.write("@HD\tVN:1.0\tSO:unsorted\n")
        sam_file.write("@SQ\tSN:test_chrom\tLN:60\n")
        sam_file.write("@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.6"
                       "\tCL:\"/iblm/netapp/home/gmcvicker/anaconda2/bin/"
                       "bowtie2-align-s --wrapper basic-0 -x "
                       "test_data/test_genome -S test_data/test.sam -U "
                       "test_data/test_1.fq\"\n")


    def write_sam_read(self, sam_file, read_name="read1#0/1",
                       flag=0, chrom="test_chrom", pos=1,
                       mapq=30, cigar="30M", rnext="*",
                       pnext=0, tlen=0, seq="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                       qual="BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                       opt_fields=["AS:i:0", "XS:i:-5", "XN:i:0",
                                   "XM:i:0", "XO:i:0", "XG:i:0", "NM:i:0",
                                   "MD:Z:30", "YT:Z:UU"]):
        sam_file.write("\t".join([read_name, "%d" % flag, chrom, "%d" % pos,
                                 "%d" % mapq, cigar, rnext, "%d" % pnext,
                                 "%d" % tlen, seq, qual,
                                  "\t".join(opt_fields)]) + "\n")


class TestReadFile(object):
    
    def test_read_snps(self):
        data = Data()
        data.setup()
        
        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)

        # check snp_index set correctly
        assert len(snp_tab.snp_index) == 100
        assert snp_tab.snp_index[9] == 0
        assert snp_tab.snp_index[19] == 1
        assert snp_tab.snp_index[99] == 2
        # only 3 values of index should be non -1
        assert np.where(snp_tab.snp_index != -1)[0].shape[0] == 3

        # check snp_allele set correctly
        assert snp_tab.snp_allele1[0] == b"A"
        assert snp_tab.snp_allele2[0] == b"C"
        assert snp_tab.snp_allele1[1] == b"T"
        assert snp_tab.snp_allele2[1] == b"G"
        assert snp_tab.snp_allele1[2] == b"A"
        assert snp_tab.snp_allele2[2] == b"T"

        # check that snp_pos set correctly
        assert snp_tab.snp_pos[0] == 10
        assert snp_tab.snp_pos[1] == 20
        assert snp_tab.snp_pos[2] == 100


    def test_read_indels(self):
        data = Data()
        data.snp_list = [(10, "A", "-"), # 1bp deletion
                         (20, "A", "ATTG"), # 3bp insertion
                         (21, "A", "T"), # not an indel
                         (3, "AAA", "A")] # 2bp deletion
        
        data.setup()

        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)

        # check snp_index set correctly
        assert len(snp_tab.snp_index) == 21
        assert snp_tab.snp_index[9] == 0
        assert snp_tab.snp_index[19] == 1
        assert snp_tab.snp_index[2] == 3
        
        # only 4 values of index should be non -1
        assert np.where(snp_tab.snp_index != -1)[0].shape[0] == 4

        # check snp_allele set correctly
        assert snp_tab.snp_allele1[0] == b"A"
        assert snp_tab.snp_allele2[0] == b""
        assert snp_tab.snp_allele1[1] == b"A"
        assert snp_tab.snp_allele2[1] == b"ATTG"
        assert snp_tab.snp_allele1[3] == b"AAA"
        assert snp_tab.snp_allele2[3] == b"A"

        # check that snp_pos set correctly
        assert snp_tab.snp_pos[0] == 10
        assert snp_tab.snp_pos[1] == 20
        assert snp_tab.snp_pos[2] == 21
        assert snp_tab.snp_pos[3] == 3


class TestGetOverlappingSNPs:
        
    def test_get_overlapping_snps_simple(self):
        """Do a simple test of getting 2 overlapping SNPs
        with a read with 30 matches"""
        data = Data()
        data.setup()

        # write a single read with all matches to SAM
        sam_file = open(data.sam_filename, "w")
        data.write_sam_header(sam_file)
        data.write_sam_read(sam_file)
        sam_file.close()

        sam_file = pysam.Samfile(data.sam_filename)
        read = next(sam_file)

        # simple case where read has only one big match segment
        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

        # check that overlapping SNPs are found and in correct locations
        assert len(snp_idx) == 2
        assert snp_idx[0] == 0
        assert snp_idx[1] == 1
        
        assert snp_read_pos[0] == 10
        assert snp_read_pos[1] == 20

        assert len(indel_idx) == 0
        assert len(indel_read_pos) == 0

        

    def test_get_overlapping_snps_intron(self):
        """Test a read spanning an intron (N in CIGAR string)"""
        data = Data()
        data.setup()

        # write a single read with intron in CIGAR (N)
        sam_file = open(data.sam_filename, "w")
        data.write_sam_header(sam_file)
        data.write_sam_read(sam_file, cigar="10M85N20M")
        sam_file.close()
        
        sam_file = pysam.Samfile(data.sam_filename)
        read = next(sam_file)

        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

        # check that overlapping SNPs are found and in correct locations
        assert len(snp_idx) == 2
        assert snp_idx[0] == 0
        assert snp_idx[1] == 2
        
        assert snp_read_pos[0] == 10
        assert snp_read_pos[1] == 15

    
    
    def test_get_overlapping_snps_softclip(self):
        """Test that soft-clipped part of read is not used"""
        data = Data()
        data.setup()

        # write a single read with softclipping on left end
        sam_file = open(data.sam_filename, "w")
        data.write_sam_header(sam_file)
        data.write_sam_read(sam_file, cigar="10S20M")
        sam_file.close()
        
        sam_file = pysam.Samfile(data.sam_filename)
        read = next(sam_file)

        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

        # check that overlapping SNPs are found and in correct locations
        assert len(snp_idx) == 2
        assert snp_idx[0] == 0
        assert snp_idx[1] == 1
        assert snp_read_pos[0] == 20
        assert snp_read_pos[1] == 30


    def test_get_overlapping_indel(self):
        """Test that indels can be correctly obtained"""
        data = Data()
        data.snp_list = [(10, "A", "-")]
        data.setup()

        # write a single read with match
        sam_file = open(data.sam_filename, "w")
        data.write_sam_header(sam_file)
        data.write_sam_read(sam_file, cigar="30M")
        sam_file.close()
        
        sam_file = pysam.Samfile(data.sam_filename)
        read = next(sam_file)

        snp_tab = snptable.SNPTable()
        snp_tab.read_file(data.snp_filename)
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)

        # check that overlapping indel found in correct location
        assert len(snp_idx) == 0
        assert len(indel_idx) == 1
        assert indel_idx[0] == 0
        assert indel_read_pos[0] == 10
        
        
        

    

