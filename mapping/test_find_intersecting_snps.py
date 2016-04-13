import glob
import gzip
import os
import subprocess
import sys


import find_intersecting_snps

def read_bam(bam):
    """
    Read a bam file into a list where each element of the list is a line from
    the bam file (with the newline stripped). The header is discarded.
    """
    res = subprocess.check_output('samtools view %s' % bam, shell=True)
    return res.strip().split('\n')



class Data(object):
    """This class creates data that can be used for the tests"""

    def __init__(self,
                 data_dir="test_data",
                 prefix="test_data/test",
                 read1_seqs =  ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"],
                 read2_seqs =  ["TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                 read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                 read2_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                 genome_seq =  ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"),
                 chrom_name = 'test_chrom',
                 snp_list = [['test_chrom', 1, "A", "C"]]):

        self.data_dir = data_dir
        self.prefix = prefix
        self.read1_seqs  = read1_seqs
        self.read1_quals = read1_quals
        self.read2_seqs  = read2_seqs
        self.read2_quals = read2_quals
        self.genome_seq = genome_seq
        self.snp_list = snp_list

        self.genome_prefix = self.prefix + "_genome"
        self.genome_filename = self.genome_prefix + ".fa"
        self.chrom_name = "test_chrom"

        self.fastq1_filename = self.prefix + "_1.fq"
        self.fastq2_filename = self.prefix + "_2.fq"

        self.sam_filename = self.prefix + ".sam"
        self.bam_sort_filename = self.prefix + ".sort.bam"
        self.bam_filename = self.prefix + ".bam"

        self.snp_dir = self.prefix + "_snps"

        # these are files that are written by find_intersecting_snps
        self.bam_keep_filename = self.prefix + ".keep.bam"
        self.bam_remap_filename = self.prefix + ".to.remap.bam"
        self.num_remap_filename = self.prefix + ".to.remap.num.gz"
        self.fastq_remap_filename = self.prefix + ".remap.fq.gz"

    
    def setup(self):
        """Create the test genome, test fastq and test SNP files"""
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        self.write_ref_genome()
        self.write_fastqs()
        self.write_snps()

        

    def cleanup(self):
        """remove files created by the tests"""
        filenames = [
            self.genome_filename,
            self.fastq1_filename,
            self.fastq2_filename,
            self.sam_filename,
            self.bam_filename,
            self.bam_sort_filename,
            self.bam_keep_filename,
            self.bam_remap_filename,
            self.num_remap_filename,
            self.fastq_remap_filename]

        index_filenames = glob.glob(self.genome_prefix + "*.bt2")
        filenames.extend(index_filenames)
        
        snp_filenames = glob.glob(self.snp_dir + "/*.snps.txt.gz")
        filenames.extend(snp_filenames)

        for fname in filenames:
            if os.path.exists(fname):
                os.remove(fname)

        if os.path.exists(self.snp_dir):
            os.rmdir(self.snp_dir)
    

    def write_ref_genome(self):
        f = open(self.genome_filename, "w")        
        f.write(">" + self.chrom_name + "\n" + self.genome_seq)
        f.close()
    
    
    def write_fastqs(self):        
        if self.read1_seqs:
            # write fastq1
            read_num = 1
            f = open(self.fastq1_filename, "w")
            for seq_str, qual_str in zip(self.read1_seqs, self.read1_quals):
                f.write("@read%d#0/1\n" % read_num)
                f.write(seq_str + "\n")
                f.write("+\n")
                f.write(qual_str + "\n")
                read_num += 1
            f.close()

        if self.read2_seqs:
            read_num = 1
            f = open(self.fastq2_filename, "w")
            for seq_str, qual_str in zip(self.read2_seqs, self.read2_quals):
                f.write("@read%d#0/2\n" % read_num)
                f.write(seq_str + "\n")
                f.write("+\n")
                f.write(qual_str + "\n")
                read_num += 1
            f.close()

            

    def index_genome_bowtie2(self):
        cmd = ['bowtie2-build', self.genome_filename, self.genome_prefix]
        # write stderr and stdout to /dev/null because bowtie2-build writes a
        # lot of output
        fnull = open(os.devnull, 'w')
        subprocess.check_call(cmd, stdout=fnull, stderr=subprocess.STDOUT)


    def map_single_bowtie2(self):
        cmd = ['bowtie2', '-x', self.genome_prefix, "-U", self.fastq1_filename,
               "-S", self.sam_filename]
        subprocess.check_call(cmd)


    def sam2bam(self):
        cmd = 'samtools view -S -b %s > %s' % \
              (self.sam_filename, self.bam_filename)
        subprocess.check_call(cmd, shell=True)


    def write_snps(self):
        files = {}

        if not os.path.exists(self.snp_dir):
            os.makedirs(self.snp_dir)
        
        for snp in self.snp_list:
            chrom, pos, allele1, allele2 = snp
            if chrom not in files:
                filename = self.snp_dir + "/" + chrom + ".snps.txt.gz"
                files[chrom] = gzip.open(filename, "wb")
            files[chrom].write("%d\t%s\t%s\n" % (pos, allele1, allele2))
        
        for f in files.values():
            f.close()





        
                      
class TestSingleEnd:
    """tests for single end read mapping"""
        
    
    def test_single_one_read_one_snp(self):
        """Simple test of whether 1 read overlapping 
        1 SNP works correctly"""
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
                        
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, max_window=100000,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        
        assert len(lines) == 4
        
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)
        
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]

        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # Verify that the .to.remap.num filename contains
        # the number of reads that need to be remapped
        #
        with gzip.open(test_data.num_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 1
        assert lines[0] == "1"
        
        test_data.cleanup()


    def test_single_one_read_one_indel(self):
        """Test whether 1 read overlapping indel works correctly"""

        #
        # Currently WASP discards reads that overlap
        # indels. We want to improve WASP to handle indels
        # correctly, but in the meantime this test just checks
        # that the read is discarded.
        #
        
        test_data = Data(snp_list=[['test_chrom', 2, "A", "CCC"]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, max_window=100000,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. Should be empty because only
        # read overlaps an indel
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        
        assert len(lines) == 0
        
        
        #
        # Verify to.remap bam is empty
        #
        new_lines = read_bam(test_data.bam_remap_filename)
        assert len(new_lines) == 1
        assert new_lines[0] == ""
        
        #
        # Verify that the keep file is empty.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # Verify that the .to.remap.num filename is empty
        #
        with gzip.open(test_data.num_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 0
        
        test_data.cleanup()
        
        
    def test_single_two_reads_one_snp(self):
        """Test whether 2 reads (one overlapping SNP, 
        one not overlapping SNP) works correctly"""

        read1_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                      "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"]
        read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                       "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"]

        test_data = Data(read1_seqs=read1_seqs,
                         read1_quals=read1_quals)

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, max_window=100000,
                                    is_paired_end=False, is_sorted=False)

        
        #
        # Verify new fastq is correct. It should only contain
        # a single read and, the first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        

        assert len(lines) == 4

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]
        
        # Verify first line in TO.REMAP bam is the same as the first
        # line in the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        
        new_lines = read_bam(test_data.bam_remap_filename)
        assert len(new_lines) == 1
        assert old_lines[0] == new_lines[0]
        
        # verify that first line in the KEEP bam is the same as the second
        # line in the input bam file
        new_lines = read_bam(test_data.bam_keep_filename)
        assert len(new_lines) == 1
        assert old_lines[1] == new_lines[0]

        #
        # Verify that the .to.remap.num filename contains
        # the number of reads that need to be remapped (1)
        #
        with gzip.open(test_data.num_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 1
        assert lines[0] == "1"
        
        
    def test_single_one_read_two_snps(self):
        """Test whether 1 read overlapping 2 SNPs works correctly"""
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
                        
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, max_window=100000,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]
        
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'G'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'G'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs


        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # Verify that the .to.remap.num filename contains
        # the number of reads that need to be remapped
        #
        with gzip.open(test_data.num_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 1
        assert lines[0] == "3"
        
        test_data.cleanup()


    def test_single_one_read_ten_snps(self):
        """Test whether 1 read overlapping 10 SNPs works correctly"""

        snp_list = [['test_chrom', x, "A", "C"] for x in range(1, 11)]
        

        test_data = Data(snp_list=snp_list)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
                        
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, max_window=100000,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. There should be 1023 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        
        assert len(lines) == 4*1023

        # get every 4th line, starting at line 1
        seqs = [lines[x] for x in range(1, 1024, 4)]

        # test a few combinations of alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'C'
        new_seq2 = "".join(l)

        # read with 3 non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'C'
        l[10] = 'C'
        new_seq3 = "".join(l)

        # read with 10 non-ref alleles
        l = list(test_data.read1_seqs[0])
        for i in range(10):
            l[i] = 'C'
        new_seq4 = "".join(l)
        
        assert len(seqs) == 1023
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs
        assert new_seq4 in seqs

        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # Verify that the .to.remap.num filename contains
        # the number of reads that need to be remapped
        #
        with gzip.open(test_data.num_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 1
        assert lines[0] == "1023"
        
        test_data.cleanup()
        


    def test_single_cli(self):
        """This test is to make sure the command line interface 
        for single-end read mapping works"""
        
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        command = ['python', 'find_intersecting_snps.py',
                   test_data.bam_filename, test_data.snp_dir]
        subprocess.check_call(command)

        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]        
        assert len(lines) == 4
        
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)
        
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]
        
        # Verify first line in to.remap bam is the same as the first
        # line in the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        
        new_lines = read_bam(test_data.bam_remap_filename)
        assert len(new_lines) == 1
        assert old_lines[0] == new_lines[0]
        
        test_data.cleanup()
