import glob
import gzip
import os
import subprocess

import sys


sys.path.append("../")

from mapping.filter_remapped_reads import *
from mapping.find_intersecting_snps import *


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
           glob.glob('test_data/test*.to.remap.num.gz') + 
           glob.glob('test_data/test_*_filtered.bam'))
    [os.remove(x) for x in fns]

class TestRun:
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

        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_single.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)

        lines = read_bam(keep_bam)
        assert len(lines) == 1

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
        
        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_paired.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 2

        cleanup()

    def test_simple_single_unmapped(self):
        """Test to make sure that if the read pair is unmapped in the remapping
        stage, it is not written to the final output file."""
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
        
        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_single_unmapped.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert lines == ['']

        cleanup()

    def test_simple_paired_unmapped(self):
        """Test to make sure that if the read pair is unmapped in the remapping
        stage, it is not written to the final output file."""
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
        
        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_paired_unmapped.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert lines == ['']

        cleanup()

    def test_simple_single_reverse(self):
        """Test to make sure that if the read is mapped correctly on the reverse
        strand in the remapping stage, it is written to the final output
        file."""
        is_paired_end = False
        max_window = 100000
        pref = 'test_data/test_single_reverse'
        file_name = pref + ".sort.bam"
        keep_file_name = pref + ".keep.bam"
        remap_name = pref + ".to.remap.bam"
        remap_num_name = pref + ".to.remap.num.gz"
        fastq_names = [pref + ".remap.fq.gz"]
        snp_dir = 'test_data/snps'
        bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
                        remap_name, remap_num_name, fastq_names, snp_dir)
        bs.run()
        
        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_single_reverse.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 1

        cleanup()

    def test_simple_paired_reverse(self):
        """Test to make sure that if the first read in the read pair is mapped
        correctly on the reverse strand in the remapping stage, it is written to
        the final output file."""
        cleanup()
        is_paired_end = True
        max_window = 100000
        pref = 'test_data/test_paired_reverse'
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
        
        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_paired_reverse.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 2

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

        keep_bam = pref + '_filtered.bam'
        run(remap_name, 'test_data/test_two_snps_single.remapped.bam', keep_bam,
            remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 1

        cleanup()

    # def test_two_snps_paired(self):
    #     # TODO: The remapped bam file should have six entries, only has two.
    #     is_paired_end = True
    #     max_window = 100000
    #     pref = 'test_data/test_paired'
    #     file_name = pref + ".sort.bam"
    #     keep_file_name = pref + ".keep.bam"
    #     remap_name = pref + ".to.remap.bam"
    #     remap_num_name = pref + ".to.remap.num.gz"
    #     fastq_names = [pref + ".remap.fq1.gz",
    #                    pref + ".remap.fq2.gz"]
    #     snp_dir = 'test_data/two_snps'
    #     bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
    #                     remap_name, remap_num_name, fastq_names, snp_dir)
    #     bs.run()
    #     
    #     keep_bam = pref + '_filtered.bam'
    #     run(remap_name, 'test_data/test_two_snps_paired.remapped.bam', keep_bam,
    #         remap_num_name, is_paired_end)

    #     lines = read_bam(keep_bam)
    #     assert len(lines) == 2

    #     cleanup()

    # def test_issue_18(self):
    #     """
    #     This was reported as a bug because one read pair that overlaps one SNP
    #     was resulting in multiple pairs of reads in the fastq files. However, it
    #     is not a bug because the reads overlap and both reads overlap the SNP.
    #     """
    #     # TODO: The remapped bam file does not correspond to the input bam file.
    #     is_paired_end = True
    #     max_window = 100000
    #     file_name = 'test_data/issue_18.bam'
    #     pref = 'test_data/test_paired'
    #     keep_file_name = pref + ".keep.bam"
    #     remap_name = pref + ".to.remap.bam"
    #     remap_num_name = pref + ".to.remap.num.gz"
    #     fastq_names = [pref + ".remap.fq1.gz",
    #                    pref + ".remap.fq2.gz"]
    #     snp_dir = 'test_data/issue_18_snps'
    #     bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
    #                     remap_name, remap_num_name, fastq_names, snp_dir)
    #     bs.run()
    #     
    #     keep_bam = pref + '_filtered.bam'
    #     run(remap_name, 'test_data/test_issue_18.remapped.bam', keep_bam,
    #         remap_num_name, is_paired_end)

    #     lines = read_bam(keep_bam)
    #     assert len(lines) == 2

    #     cleanup()
    # 
    # def test_issue_23(self):
    #     """
    #     This was reported as a bug because WASP said the read pair mapped
    #     correctly yet wasn't written to the output file.
    #     """
    #     # TODO: I think the test above was copied down but the test for but the
    #     # test for issue 23 was not implemented. Maybe I accidentally copied
    #     # over it at some point.
    #     is_paired_end = True
    #     max_window = 100000
    #     file_name = 'test_data/issue_18.bam'
    #     pref = 'test_data/test_paired'
    #     keep_file_name = pref + ".keep.bam"
    #     remap_name = pref + ".to.remap.bam"
    #     remap_num_name = pref + ".to.remap.num.gz"
    #     fastq_names = [pref + ".remap.fq1.gz",
    #                    pref + ".remap.fq2.gz"]
    #     snp_dir = 'test_data/issue_18_snps'
    #     bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
    #                     remap_name, remap_num_name, fastq_names, snp_dir)
    #     bs.run()
    #     
    #     keep_bam = pref + '_filtered.bam'
    #     run(remap_name, 'test_data/test_issue_18.remapped.bam', keep_bam,
    #         remap_num_name, is_paired_end)

    #     lines = read_bam(keep_bam)
    #     assert len(lines) == 2

    #     cleanup()
    
    def test_bad_first_paired(self):
        """Test whether the correct read pairs are output if the first read pair
        is incorrectly remapped."""
        remap_name = 'test_data/bad_first_paired/bad_first_to_remap.bam'
        keep_bam = 'test_data/test_bad_first_filtered.bam'
        remapped_bam = 'test_data/bad_first_paired/bad_first_remapped.bam' 
        remap_num_name = 'test_data/bad_first_paired/bad_first.to.remap.num.gz'
        is_paired_end = True
        run(remap_name, remapped_bam, keep_bam, remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 6

        cleanup()

    def test_bad_middle_paired(self):
        """Test whether the correct read pairs are output if a read pair in the
        middle is incorrectly remapped."""
        remap_name = 'test_data/bad_middle_paired/bad_middle_to_remap.bam'
        keep_bam = 'test_data/test_bad_middle_filtered.bam'
        remapped_bam = 'test_data/bad_middle_paired/bad_middle_remapped.bam' 
        remap_num_name = 'test_data/bad_middle_paired/bad_middle.to.remap.num.gz'
        is_paired_end = True
        run(remap_name, remapped_bam, keep_bam, remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 6

    def test_bad_last_paired(self):
        """Test whether the correct read pairs are output if the last read pair
        is incorrectly remapped."""
        remap_name = 'test_data/bad_last_paired/bad_last_to_remap.bam'
        keep_bam = 'test_data/test_bad_last_filtered.bam'
        remapped_bam = 'test_data/bad_last_paired/bad_last_remapped.bam' 
        remap_num_name = 'test_data/bad_last_paired/bad_last.to.remap.num.gz'
        is_paired_end = True
        run(remap_name, remapped_bam, keep_bam, remap_num_name, is_paired_end)
        
        lines = read_bam(keep_bam)
        assert len(lines) == 6

        cleanup()

class TestCLI:
    def test_simple_single_cli(self):
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

        keep_bam = pref + '_filtered.bam'
        c = ('python ../filter_remapped_reads.py {} {} {} {}'.format(
            remap_name, 'test_data/test_single.remapped.bam', keep_bam,
            remap_num_name))
        subprocess.check_call(c, shell=True)

        lines = read_bam(keep_bam)
        assert len(lines) == 1

        cleanup()
