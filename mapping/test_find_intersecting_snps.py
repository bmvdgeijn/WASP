import glob
import gzip
import os
import subprocess

from find_intersecting_snps import *

def read_bam(bam):
    """
    Read a bam file into a list where each element of the list is a line from
    the bam file (with the newline stripped). The header is discarded.
    """
    res = subprocess.check_output('samtools view {}'.format(bam), shell=True)
    return res.strip().split('\n')

def cleanup():
    fns = (glob.glob('test_data/test*.keep.bam') +
           glob.glob('test_data/test*.remap.fq*.gz') + 
           glob.glob('test_data/test*.to.remap.bam') + 
           glob.glob('test_data/test*.to.remap.num.gz'))
    [os.remove(x) for x in fns]

class TestSNP:
    def test_init(self):
        """Test to see whether __init__ is working as expected."""
        snp = SNP('12670\tG\tC\n')
        assert snp.pos == 12670 - 1
        assert snp.alleles == ['G', 'C']
        assert snp.ptype == 'snp'
        assert snp.max_len == 1

    def test_add_allele(self):
        """Test to see whether we can add an allele."""
        snp = SNP('12670\tG\tC\n')
        snp.add_allele(['A'])
        assert snp.alleles == ['G', 'C', 'A']
    
    def test_add_allele_multiple(self):
        """Test to see whether we can add multiple alleles."""
        snp = SNP('12670\tG\tC\n')
        snp.add_allele(['A', 'T'])
        assert snp.alleles == ['G', 'C', 'A', 'T']

    # TODO: tests for adding insertions and deletions

class TestBamScanner:
    def test_init_single(self):
        is_paired_end = False
        max_window = 100000
        pref = 'test_data/test_single'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq.gz"]
        snp_dir = 'test_data/snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        assert bs.max_window == len(bs.snp_table)
        cleanup()

    def test_init_paired(self):
        is_paired_end = True
        max_window = 100000
        pref = 'test_data/test_paired'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq1.gz",
                       pref + ".remap.fq2.gz"]
        snp_dir = 'test_data/snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        assert bs.max_window == len(bs.snp_table)
        cleanup()

    def test_simple_single(self):
        is_paired_end = False
        max_window = 100000
        pref = 'test_data/test_single'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq.gz"]
        snp_dir = 'test_data/snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()

        # Verify fastq is correct. The second base of the first read should be
        # switched from a C to an A.  (14538   C   A)
        seq = ('CATCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        qual = ('BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFBFF<FFFFFFFFFBFFFFFFFFFFFFFFFFFFF')
        with gzip.open('test_data/test_single.remap.fq.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4
        assert lines[1] == seq
        assert lines[3] == qual

        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam('test_data/test_single.sort.bam')
        new_lines = read_bam('test_data/test_single.to.remap.bam')
        assert old_lines == new_lines
        
        cleanup()

    def test_simple_paired(self):
        is_paired_end = True
        max_window = 100000
        pref = 'test_data/test_paired'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq1.gz",
                       pref + ".remap.fq2.gz"]
        snp_dir = 'test_data/snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()
        # The second base should be switched from a C to an A.
        # 14538   C   A
        seq = ('CATCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        qual = ('BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFBFF<FFFFFFFFFBFFFFFFFFFFFFFFFFFFF')
        with gzip.open('test_data/test_paired.remap.fq1.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4
        assert lines[1] == seq
        assert lines[3] == qual

        # Shouldn't be any changes to the second read.
        seq = ('TCATGGAGCCCCCTACGATTCCCAGTCGTCCTCGTCCTCCTCTGCCTGTGGCTGCTGCGGTGG'
               'CGGCAGAGGAGGGATGGAGTCTGACACGCGGGCAAAG')
        qual = ('B//FF77BB<7/BB<7FBFFF<</FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB')
        with gzip.open('test_data/test_paired.remap.fq2.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4
        assert lines[1] == bs.reverse_complement(seq)
        assert lines[3] == qual
        
        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam('test_data/test_paired.sort.bam')
        new_lines = read_bam('test_data/test_paired.to.remap.bam')
        assert old_lines == new_lines
        
        cleanup()

    def test_two_snps_single(self):
        is_paired_end = False
        max_window = 100000
        pref = 'test_data/test_single'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq.gz"]
        snp_dir = 'test_data/two_snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()

        # Verify fastq is correct. The second base of the first read should be
        # switched from a C to an A and the third base of the first read should
        # be switched from T to G. (14538   C   A, 14539    T   G)
        with gzip.open('test_data/test_single.remap.fq.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12
        seq = ('CATCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[1] == seq
        seq = ('CCGCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[5] == seq
        seq = ('CAGCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[9] == seq
        qual = ('BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFBFF<FFFFFFFFFBFFFFFFFFFFFFFFFFFFF')
        for i in [3, 7, 11]:
            assert lines[i] == qual

        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam('test_data/test_single.sort.bam')
        new_lines = read_bam('test_data/test_single.to.remap.bam')
        assert old_lines == new_lines
        
        cleanup()

    def test_two_snps_paired(self):
        is_paired_end = True
        max_window = 100000
        pref = 'test_data/test_paired'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq1.gz",
                       pref + ".remap.fq2.gz"]
        snp_dir = 'test_data/two_snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()
        
        # Verify fastq is correct. The second base of the first read should be
        # switched from a C to an A and the third base of the first read should
        # be switched from T to G. (14538   C   A, 14539    T   G)
        with gzip.open('test_data/test_paired.remap.fq1.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12
        seq = ('CATCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[1] == seq
        seq = ('CCGCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[5] == seq
        seq = ('CAGCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCC'
               'TCCCAAGGAAGTAGGTCTGAGCAGCTTGTCCTGGCTGT')
        assert lines[9] == seq
        qual = ('BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFBFF<FFFFFFFFFBFFFFFFFFFFFFFFFFFFF')
        for i in [3, 7, 11]:
            assert lines[i] == qual

        # Shouldn't be any changes to the second read.
        seq = ('TCATGGAGCCCCCTACGATTCCCAGTCGTCCTCGTCCTCCTCTGCCTGTGGCTGCTGCGGTGG'
               'CGGCAGAGGAGGGATGGAGTCTGACACGCGGGCAAAG')
        qual = ('B//FF77BB<7/BB<7FBFFF<</FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
                'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB')
        with gzip.open('test_data/test_paired.remap.fq2.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12
        for i in [1, 5, 9]:
            assert lines[i] == bs.reverse_complement(seq)
            assert lines[i + 2] == qual
        
        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam('test_data/test_paired.sort.bam')
        new_lines = read_bam('test_data/test_paired.to.remap.bam')
        assert old_lines == new_lines
        
        cleanup()

    def test_2015_06_18_bug(self):
        """
        This was reported as a bug because one read pair that overlaps one SNP
        was resulting in multiple pairs of reads in the fastq files. However, it
        is not a bug because the reads overlap and both reads overlap the SNP.
        """
        is_paired_end = True
        max_window = 100000
        file_name = 'test_data/2015_06_18_bug.bam'
        pref = 'test_data/test_paired'
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq1.gz",
                       pref + ".remap.fq2.gz"]
        snp_dir = 'test_data/2015_06_18_bug_snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()
        
        # Verify fastq are correct. 
        with gzip.open('test_data/test_paired.remap.fq1.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12
        seq = ('CGAGCGCTCACTCAATATCACGAGAACAGCAAGGGGGAAGTCGGCCCCCANGAGCCAATNACC'
               'TCCCANNNGGTCCCTCCCACAACACTGGGAATTACAA')
        assert lines[1] == seq
        seq = ('CGAGCGCTCACTCAATATCACAAGAACAGCAAGGGGGAAGTCGGCCCCCANGAGCCAATNACC'
               'TCCCANNNGGTCCCTCCCACAACACTGGGAATTACAA')
        assert lines[5] == seq
        seq = ('CGAGCGCTCACTCAATATCACAAGAACAGCAAGGGGGAAGTCGGCCCCCANGAGCCAATNACC'
               'TCCCANNNGGTCCCTCCCACAACACTGGGAATTACAA')
        assert lines[9] == seq
        qual = ('</<<<FBFFFFFFFFFFFBBFFFFFFFF/BFFFFF/BF<BBF//FBFFFB#<<<BFF//#<<'
                'FBFBF/###<7<<FFFFFFFFFFFF<<B//B7F7BBFB')
        for i in [3, 7, 11]:
            assert lines[i] == qual

        with gzip.open('test_data/test_paired.remap.fq2.gz') as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12
        seq = ('ATCACAAGAACAGCAAGGGGGAAGTCGGCCCCCATGAGCCAATCACCTCCCACCAGGTCCCTC'
               'CCACAACACTGGGAATTACAATTTNACATNACATTTG')
        assert lines[1] == bs.reverse_complement(seq)
        seq = ('ATCACGAGAACAGCAAGGGGGAAGTCGGCCCCCATGAGCCAATCACCTCCCACCAGGTCCCTC'
               'CCACAACACTGGGAATTACAATTTNACATNACATTTG')
        assert lines[5] == bs.reverse_complement(seq)
        seq = ('ATCACAAGAACAGCAAGGGGGAAGTCGGCCCCCATGAGCCAATCACCTCCCACCAGGTCCCTC'
               'CCACAACACTGGGAATTACAATTTNACATNACATTTG')
        assert lines[9] == bs.reverse_complement(seq)
        qual = ('#B<BBFFFFFFFFFFFFFFFF<FFFFBFFFFFFFFF<F/FFFFF<FFFFFFFF<FFFFFFBF'
                'F<<FFFFFFFFFFFFFFFFFFFB<B#F<<<#FFBBBBB')
        for i in [3, 7, 11]:
            assert lines[i] == qual

        cleanup()

# TODO: Add in the read pairs that I think are causing bugs.
