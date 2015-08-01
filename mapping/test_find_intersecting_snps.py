import glob
import os

from find_intersecting_snps import *

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
        """Test to see whether we can multiple alleles."""
        snp = SNP('12670\tG\tC\n')
        snp.add_allele(['A', 'T'])
        assert snp.alleles == ['G', 'C', 'A', 'T']

    #TODO: tests for adding insertions and deletions

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

# TODO: Add in the read pairs that I think are causing bugs.
