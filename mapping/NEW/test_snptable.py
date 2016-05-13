

import snptable
import gzip
import os

import numpy as np

import pysam


def write_sam_header(sam_file):
    sam_file.write("@HD\tVN:1.0\tSO:unsorted\n")
    sam_file.write("@SQ\tSN:test_chrom\tLN:60\n")
    sam_file.write("@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.6\tCL:\"/iblm/netapp/home/gmcvicker/anaconda2/bin/bowtie2-align-s --wrapper basic-0 -x test_data/test_genome -S test_data/test.sam -U test_data/test_1.fq\"\n")


def write_sam_read(sam_file, read_name="read1#0/1",
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
        snp_tab = snptable.SNPTable()
        
        filename = "test_data/snp_tab.txt.gz"
        if not os.path.exists("test_data"):
            os.makedirs("test_data")
        snp_file = gzip.open(filename, "wb")
                
        snp_file.write("10 A C\n")
        snp_file.write("20 T G\n")
        snp_file.write("100 A T\n")
        snp_file.close()
        
        snp_tab.read_file(filename)

        # check snp_index set correctly
        assert len(snp_tab.snp_index) == 100
        assert snp_tab.snp_index[9] == 0
        assert snp_tab.snp_index[19] == 1
        assert snp_tab.snp_index[99] == 2
        # only 3 values of index should be non -1
        assert np.where(snp_tab.snp_index != -1)[0].shape[0] == 3


        # check snp_allele set correctly
        assert snp_tab.snp_allele1[0] == "A"
        assert snp_tab.snp_allele2[0] == "C"
        assert snp_tab.snp_allele1[1] == "T"
        assert snp_tab.snp_allele2[1] == "G"
        assert snp_tab.snp_allele1[2] == "A"
        assert snp_tab.snp_allele2[2] == "T"

        # check that snp_pos set correctly
        assert snp_tab.snp_pos[0] == 10
        assert snp_tab.snp_pos[1] == 20
        assert snp_tab.snp_pos[2] == 100

        
        # TODO: need to add test for indels


class TestGetOverlappingSNPs:
        
    def test_get_overlapping_snps(self):

        if not os.path.exists("test_data"):
            os.makedirs("test_data")

        snp_filename = "test_data/snp_tab.txt.gz"
        snp_file = gzip.open(snp_filename, "wb")
        snp_file.write("10 A C\n")
        snp_file.write("20 T G\n")
        snp_file.write("100 A T\n")
        snp_file.close()

        # write a single read with all matches to SAM
        sam_filename = "test_data/test.sam"
        sam_file = open(sam_filename, "w")
        write_sam_header(sam_file)
        write_sam_read(sam_file)
        sam_file.close()

        sam_file = pysam.Samfile(sam_filename)
        read = sam_file.next()

        # simple case where read has only one big match segment
        snp_tab = snptable.SNPTable()
        snp_tab.read_file(snp_filename)
        snp_idx_list, read_idx_list = snp_tab.get_overlapping_snps(read)

        # check that overlapping SNPs are found and in correct locations
        assert len(snp_idx_list) == 2
        assert snp_idx_list[0] == 0
        assert snp_idx_list[1] == 1
        
        assert read_idx_list[0] == 9
        assert read_idx_list[1] == 19

        
        #### TODO: test more complicated CIGAR strings:
        ####       indels, soft-clipping, etc.

        
        
        

        
        

    

