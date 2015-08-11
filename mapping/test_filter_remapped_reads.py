import glob
import gzip
import os
import subprocess

from find_intersecting_snps import *
from filter_remapped_reads import *

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
    # def test_simple_single(self):
    #     is_paired_end = False
    #     max_window = 100000
    #     pref = 'test_data/test_single'
    #     file_name = pref + ".sort.bam"
    #     keep_file_name = pref + ".keep.bam"
    #     remap_name = pref + ".to.remap.bam"
    #     remap_num_name = pref + ".to.remap.num.gz"
    #     fastq_names = [pref + ".remap.fq.gz"]
    #     snp_dir = 'test_data/snps'
    #     bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
    #                     remap_name, remap_num_name, fastq_names, snp_dir)
    #     bs.run()

    #     keep_bam = pref + '_filtered.bam'
    #     run(remap_name, 'test_data/test_single.remapped.bam', keep_bam,
    #         remap_num_name, is_paired_end)

    #     lines = read_bam(keep_bam)
    #     import pdb
    #     pdb.set_trace()

    #     cleanup()
    
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

    # def test_simple_paired_unmapped(self):
    #     """Test to make sure that if the read pair is unmapped in the remapping
    #     stage, it is not written to the final output file."""
    #     is_paired_end = True
    #     max_window = 100000
    #     pref = 'test_data/test_paired'
    #     file_name = pref + ".sort.bam"
    #     keep_file_name = pref + ".keep.bam"
    #     remap_name = pref + ".to.remap.bam"
    #     remap_num_name = pref + ".to.remap.num.gz"
    #     fastq_names = [pref + ".remap.fq1.gz",
    #                    pref + ".remap.fq2.gz"]
    #     snp_dir = 'test_data/snps'
    #     bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
    #                     remap_name, remap_num_name, fastq_names, snp_dir)
    #     bs.run()
    #     
    #     keep_bam = pref + '_filtered.bam'
    #     run(remap_name, 'test_data/test_paired_unmapped.remapped.bam', keep_bam,
    #         remap_num_name, is_paired_end)
    #     
    #     lines = read_bam(keep_bam)
    #     assert lines == ['']

    #     cleanup()

#     def test_two_snps_single(self):
#         is_paired_end = False
#         max_window = 100000
#         pref = 'test_data/test_single'
#         file_name = pref + ".sort.bam"
#         keep_file_name = pref + ".keep.bam"
#         remap_name = pref + ".to.remap.bam"
#         remap_num_name = pref + ".to.remap.num.gz"
#         fastq_names = [pref + ".remap.fq.gz"]
#         snp_dir = 'test_data/two_snps'
#         bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
#                         remap_name, remap_num_name, fastq_names, snp_dir)
#         bs.run()
# 
#         keep_bam = pref + '_filtered.bam'
#         run(remap_name, 'test_data/test_single_two_snps.remapped.bam', keep_bam,
#             remap_num_name, is_paired_end)
# 
#         cleanup()
# 
#     def test_two_snps_paired(self):
#         is_paired_end = True
#         max_window = 100000
#         pref = 'test_data/test_paired'
#         file_name = pref + ".sort.bam"
#         keep_file_name = pref + ".keep.bam"
#         remap_name = pref + ".to.remap.bam"
#         remap_num_name = pref + ".to.remap.num.gz"
#         fastq_names = [pref + ".remap.fq1.gz",
#                        pref + ".remap.fq2.gz"]
#         snp_dir = 'test_data/two_snps'
#         bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
#                         remap_name, remap_num_name, fastq_names, snp_dir)
#         bs.run()
#         
#         keep_bam = pref + '_filtered.bam'
#         run(remap_name, 'test_data/test_paired_two_snps.remapped.bam', keep_bam,
#             remap_num_name, is_paired_end)
# 
#         cleanup()
# 
#     def test_2015_06_18_bug(self):
#         """
#         This was reported as a bug because one read pair that overlaps one SNP
#         was resulting in multiple pairs of reads in the fastq files. However, it
#         is not a bug because the reads overlap and both reads overlap the SNP.
#         """
#         is_paired_end = True
#         max_window = 100000
#         file_name = 'test_data/2015_06_18_bug.bam'
#         pref = 'test_data/test_paired'
#         keep_file_name = pref + ".keep.bam"
#         remap_name = pref + ".to.remap.bam"
#         remap_num_name = pref + ".to.remap.num.gz"
#         fastq_names = [pref + ".remap.fq1.gz",
#                        pref + ".remap.fq2.gz"]
#         snp_dir = 'test_data/2015_06_18_bug_snps'
#         bs = BamScanner(is_paired_end, max_window, file_name, keep_file_name,
#                         remap_name, remap_num_name, fastq_names, snp_dir)
#         bs.run()
#         
#         keep_bam = pref + '_filtered.bam'
#         run(remap_name, 'test_data/test_2015_06_18_bug.remapped.bam', keep_bam,
#             remap_num_name, is_paired_end)
# 
#         cleanup()
