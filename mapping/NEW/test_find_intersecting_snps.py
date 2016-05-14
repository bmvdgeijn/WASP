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
        self.fastq_remap_filename = self.prefix + ".remap.fq.gz"
        self.fastq1_remap_filename = self.prefix + ".remap.fq1.gz"
        self.fastq2_remap_filename = self.prefix + ".remap.fq2.gz"


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
            self.fastq_remap_filename,
            self.fastq1_remap_filename,
            self.fastq2_remap_filename]

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


    def map_paired_bowtie2(self):
        cmd = ['bowtie2', '-x', self.genome_prefix, "-1", self.fastq1_filename,
               "-2", self.fastq2_filename, "-S", self.sam_filename]
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
                                    test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

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
                                    test_data.snp_dir,
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
                                    test_data.snp_dir, 
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


    def test_single_one_read_two_snps(self):
        """Test whether 1 read overlapping 2 SNPs works correctly"""
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]])

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir,
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
                                    test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False,
                                    max_seqs=10)

        #
        # Verify new fastq is correct. There should be no reads,
        # because reads with greater than 10 allelic combinations
        # are thrown out
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 0

        #
        # Verify to.remap bam is empty
        #
        lines = read_bam(test_data.bam_remap_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        #
        # re-run find intersecting SNPs but allow a max of 1024
        # allelic combinations (we expect 1023 new seqs with 10
        # bi-allelic SNPs)
        #
        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, 
                                    is_paired_end=False, is_sorted=False,
                                    max_seqs=1024)

        #
        # Verify new fastq is correct. There should be 1023 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4*1023

        # get every 4th line, which correspond to sequences starting at line 1
        seqs = [lines[x] for x in range(1, len(lines), 4)]

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
        l[9] = 'C'
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

        test_data.cleanup()



    def test_single_cli(self):
        """Make sure the command line interface
        for single-end read mapping"""

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







class TestPairedEnd:
    """tests for paired end read mapping"""

    def test_paired_two_reads_one_snp(self):
        """Simple test of whether 1 PE read with one end overlapping
        1 SNP works correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "AAAAAAATTTAAAA"]
        read2_seqs = ["AAGAAACAACACAA",
                      "AAGAAACAACACAA"]
        
        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read1[1]                      AAAAAAATTTAAAA
        # SNP                            ^
        genome_seq =  ("AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT")
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890
        
        snp_list = [['test_chrom', 18, "A", "C"]]
        
        test_data = Data(genome_seq=genome_seq,
                         read1_seqs=read1_seqs,
                         read2_seqs=read2_seqs,
                         read1_quals=read1_quals,
                         read2_quals=read2_quals,
                         snp_list=snp_list)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_paired_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 is correct.
        #
        with gzip.open(test_data.fastq1_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        l = list(test_data.read1_seqs[0])
        # last base of first read should be changed from A to C
        l[13] = 'C'
        new_seq = "".join(l)
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]

        # second base of second read should be changed from A to C
        l = list(test_data.read1_seqs[1])
        l[1] = "C"
        new_seq = "".join(l)
        assert(lines[5] == new_seq)
        assert(lines[7] == test_data.read1_quals[1])

        #
        # verify fastq2 is correct
        #
        with gzip.open(test_data.fastq2_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        # bases should be the same for the second half of
        # the reads, since no SNP overlap
        assert lines[1] == test_data.read2_seqs[0]
        assert lines[3] == test_data.read2_quals[0]
        assert lines[5] == test_data.read2_seqs[1]
        assert lines[7] == test_data.read2_quals[1]
        
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

        test_data.cleanup()


        
    def test_paired_two_interleaved_reads_one_snp(self):
        """Test whether PE reads still work correctly 
        when read pairs are interleaved"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "AAAAAAATTTAAAA"]
        read2_seqs = ["AAGAAACAACACAA",
                      "AAAAATAAAAAATA"]
        
        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read1[1]                      AAAAAAATTTAAAA
        # SNP                            ^
        genome_seq =  ("AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT")
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # POS           123456789012345678901234567890
        
        snp_list = [['test_chrom', 18, "A", "C"]]
        
        test_data = Data(genome_seq=genome_seq,
                         read1_seqs=read1_seqs,
                         read2_seqs=read2_seqs,
                         read1_quals=read1_quals,
                         read2_quals=read2_quals,
                         snp_list=snp_list)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_paired_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 is correct.
        #
        with gzip.open(test_data.fastq1_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        l = list(test_data.read1_seqs[0])
        # last base of first read should be changed from A to C
        l[13] = 'C'
        new_seq = "".join(l)
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]

        # second base of second read should be changed from A to C
        l = list(test_data.read1_seqs[1])
        l[1] = "C"
        new_seq = "".join(l)
        assert(lines[5] == new_seq)
        assert(lines[7] == test_data.read1_quals[1])

        #
        # verify fastq2 is correct
        #
        with gzip.open(test_data.fastq2_remap_filename) as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        # bases should be the same for the second half of
        # the reads, since no SNP overlap
        assert lines[1] == test_data.read2_seqs[0]
        assert lines[3] == test_data.read2_quals[0]
        assert lines[5] == test_data.read2_seqs[1]
        assert lines[7] == test_data.read2_quals[1]
        
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

        test_data.cleanup()


   
    def test_paired_two_reads_one_indel(self):
        """Test of whether 2 PE reads with one end of one read
        overlapping indel works correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "AAAAAAATTTAAAA"]
        read2_seqs = ["AAGAAACAACACAA",
                      "AAAAATAAAAAATA"]
        
        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        #                       10        20        30
        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read1[1]                      AAAAAAATTTAAAA
        # SNP                            ^
        genome_seq =  ("AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT")
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # INDEL                              ^         
        # POS           123456789012345678901234567890
        #                       40        50
        
        snp_list = [['test_chrom', 19, "A", "C"],
                    ['test_chrom', 52, "G", "GTTA"]]
        
        test_data = Data(genome_seq=genome_seq,
                         read1_seqs=read1_seqs,
                         read2_seqs=read2_seqs,
                         read1_quals=read1_quals,
                         read2_quals=read2_quals,
                         snp_list=snp_list)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_paired_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        # Currently reads overlapping indels are thrown out
        expect_reads = set([("ACAAAAATTTAAAA", "AAAAATAAAAAATA")])

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename) as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename) as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == len(expect_reads) * 4
                           
        for i in range(1, len(expect_reads), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0

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

        test_data.cleanup()

        

    def test_paired_two_reads_two_snps(self):
        """Test of whether 2 PE reads with both ends overlapping
        SNPs work correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "AAAAAAATTTAAAA"]
        read2_seqs = ["AAGAAACAACACAA",
                      "AAAAATAAAAAATA"]
        
        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        #                       10        20        30
        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read1[1]                      AAAAAAATTTAAAA
        # SNP                            ^
        genome_seq =  ("AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT")
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # SNP                   ^           ^         
        # POS           123456789012345678901234567890
        #                       40        50
        
        snp_list = [['test_chrom', 18, "A", "C"],
                    ['test_chrom', 39, "A", "G"],
                    ['test_chrom', 51, "G", "T"]]
        
        test_data = Data(genome_seq=genome_seq,
                         read1_seqs=read1_seqs,
                         read2_seqs=read2_seqs,
                         read1_quals=read1_quals,
                         read2_quals=read2_quals,
                         snp_list=snp_list)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_paired_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        expect_reads = set([("AACGAAAAGGAGAC", "AAGAAACAACACAA"),
                            ("AACGAAAAGGAGAC", "AAGAAACAAAACAA"),
                            ("AACGAAAAGGAGAA", "AAGAAACAAAACAA"),
                            ("ACAAAAATTTAAAA", "AAAAATAAAAAATA"),
                            ("ACAAAAATTTAAAA", "AAAAATACAAAATA"),
                            ("AAAAAAATTTAAAA", "AAAAATACAAAATA")])

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename) as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename) as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == len(expect_reads) * 4
        for i in range(1, len(expect_reads), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0

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

        test_data.cleanup()

        # TODO: test when only one half of read maps

