import glob
import gzip
import os
import os.path
import subprocess
import sys
import tables
import numpy as np
import pysam

import find_intersecting_snps

def read_bam(bam):
    """
    Read a bam file into a list where each element of the list is a line from
    the bam file (with the newline stripped). The header is discarded.
    """
    res = subprocess.check_output('samtools view %s' % bam, shell=True)
    return res.decode("utf-8").strip().split('\n')





class Data(object):
    """This class creates data that can be used for the tests"""

    def __init__(self,
                 data_dir="test_data",
                 prefix="test_data/test",
                 output_prefix=None,
                 read1_seqs =  ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"],
                 read2_seqs =  ["TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                 read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                 read2_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                 genome_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                 chrom_names = ['test_chrom'],
                 read1_names = None,
                 read2_names = None,
                 snp_list = [['test_chrom', 1, "A", "C"]],
                 hap_samples = ["samp1", "samp2", "samp3", "samp4"],
                 haplotypes = [[0, 1, 0, 1]],
                 haplotypes_phase = False):

        if output_prefix is None:
            self.output_prefix = prefix
        else:
            self.output_prefix = output_prefix
        
        self.data_dir = data_dir
        self.prefix = prefix
        self.read1_seqs  = read1_seqs
        self.read1_quals = read1_quals
        self.read2_seqs  = read2_seqs
        self.read2_quals = read2_quals
        self.genome_seqs = list(genome_seqs)
        self.snp_list = list(snp_list)

        self.read1_names = read1_names
        self.read2_names = read2_names
        
        self.genome_prefix = self.prefix + "_genome"
        self.genome_filename = self.genome_prefix + ".fa"
        self.chrom_names = list(chrom_names)

        self.hap_samples = hap_samples
        self.haplotypes = haplotypes
        if haplotypes_phase is False:
            # if the haplotypes_phase arg was not defined, explicitly assume
            # all haps are phased
            haplotypes_phase = []
            for snp_hap in self.haplotypes:
                haplotypes_phase.append([1] * int(len(snp_hap)/2))
        self.haplotypes_phase = haplotypes_phase

        self.fastq1_filename = self.prefix + "_1.fq"
        self.fastq2_filename = self.prefix + "_2.fq"

        self.sam_filename = self.prefix + ".sam"
        self.bam_sort_filename = self.prefix + ".sort.bam"
        self.bam_filename = self.prefix + ".bam"

        self.snp_dir = self.prefix + "_snps"

        # these are files that are written by find_intersecting_snps
        self.bam_keep_filename = self.output_prefix + ".keep.bam"
        self.bam_remap_filename = self.output_prefix + ".to.remap.bam"
        self.fastq_remap_filename = self.output_prefix + ".remap.fq.gz"
        self.fastq1_remap_filename = self.output_prefix + ".remap.fq1.gz"
        self.fastq2_remap_filename = self.output_prefix + ".remap.fq2.gz"
        
        self.snp_tab_filename = self.prefix + "_snp_tab.h5"
        self.snp_index_filename = self.prefix + "_snp_index.h5"
        self.haplotype_filename = self.prefix + "_haplotype.h5"
        


    def setup(self):
        """Create the test genome, test fastq and test SNP files"""
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        self.write_ref_genome()
        self.write_fastqs()
        self.write_snps()
        self.write_h5_files()



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
            self.fastq2_remap_filename,
            self.snp_index_filename,
            self.snp_tab_filename,
            self.haplotype_filename]

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
        for chrom_name, genome_seq in zip(self.chrom_names, self.genome_seqs):
            f.write(">" + chrom_name + "\n" + genome_seq + "\n")
        f.close()


    def write_fastqs(self):        
        if self.read1_seqs:
            # write fastq1
            
            if self.read1_names is None:
                names = ["read%d" % (x+1) for x in range(len(self.read1_seqs))]
            else:
                names = self.read1_names

            f = open(self.fastq1_filename, "w")
            i = 0
            for seq_str, qual_str in zip(self.read1_seqs, self.read1_quals):
                f.write("@%s\n" % names[i])
                f.write(seq_str + "\n")
                f.write("+%s\n" % names[i])
                f.write(qual_str + "\n")
                i += 1
                
            f.close()

        if self.read2_seqs:
            if self.read2_names is None:
                names = ["read%d" % (x+1) for x in range(len(self.read2_seqs))]
            else:
                names = self.read2_names

            f = open(self.fastq2_filename, "w")
            i = 0
            for seq_str, qual_str in zip(self.read2_seqs, self.read2_quals):
                f.write("@%s\n" % names[i])
                f.write(seq_str + "\n")
                f.write("+%s\n" % names[i])
                f.write(qual_str + "\n")
                i += 1
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


    def index_genome_bwa(self):
        cmd = ['bwa', 'index', self.genome_filename]
        subprocess.check_call(cmd)

    def map_single_bwa(self):
        cmd = "bwa aln -q 15 %s %s > %s.sai" % (self.genome_filename,
                                            self.fastq1_filename,
                                            self.output_prefix)
        subprocess.check_call(cmd, shell=True)

        cmd = "bwa samse %s %s.sai %s > %s" % (self.genome_filename,
                                               self.output_prefix,
                                               self.fastq1_filename,
                                               self.sam_filename)

        subprocess.check_call(cmd, shell=True)


    def map_paired_bwa(self):
        cmd = "bwa aln -q 15 %s %s > %s_1.sai" % (self.genome_filename,
                                                  self.fastq1_filename,
                                                  self.output_prefix)
        subprocess.check_call(cmd, shell=True)

        cmd = "bwa aln -q 15 %s %s > %s_2.sai" % (self.genome_filename,
                                                  self.fastq2_filename,
                                                  self.output_prefix)
        subprocess.check_call(cmd, shell=True)

        
        cmd = "bwa sampe %s %s_1.sai %s_2.sai %s %s > %s" % \
              (self.genome_filename,
               self.output_prefix, self.output_prefix,
               self.fastq1_filename, self.fastq2_filename,
               self.sam_filename)

        subprocess.check_call(cmd, shell=True)

        
        

        
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
                files[chrom] = gzip.open(filename, "wt")
            files[chrom].write("%d\t%s\t%s\n" % (pos, allele1, allele2))

        for f in list(files.values()):
            f.close()



    def write_hap_samples(self, h5f):
        """Write tables containing sample names to HDF5 file"""
        class SamplesTab(tables.IsDescription):
            name = tables.StringCol(64)

        for chrom_name in self.chrom_names:
            table = h5f.create_table(h5f.root, "samples_%s" % chrom_name,
                                    SamplesTab)

            for samp in self.hap_samples:
                row = table.row
                row['name'] = samp
                row.append()

            table.flush()
        
        
            
    def write_snp_tab_h5(self):
        snp_tab_h5 = tables.open_file(self.snp_tab_filename, "w")

        class SNPTab(tables.IsDescription):
            name = tables.StringCol(16)
            pos = tables.Int64Col()
            allele1 = tables.StringCol(100)
            allele2 = tables.StringCol(100)

        chrom_tables = {}
        snp_num = 0
        for snp in self.snp_list:
            if snp[0] in chrom_tables:
                table = chrom_tables[snp[0]]
            else:
                table = snp_tab_h5.create_table(snp_tab_h5.root, snp[0], SNPTab)
                chrom_tables[snp[0]] = table

            row = table.row
            snp_num += 1
            row['name'] = "snp%d" % snp_num
            row['pos'] = snp[1]
            row['allele1'] = snp[2]
            row['allele2'] = snp[3]
            row.append()
            table.flush()

        self.write_hap_samples(snp_tab_h5)
        
        snp_tab_h5.close()

        
    def get_chrom_lengths(self):
        chrom_lengths = {}
        
        for chrom_name, genome_seq in zip(self.chrom_names, self.genome_seqs):
            genome_seq = genome_seq.replace("\n", "").replace(" ", "")
            chrom_lengths[chrom_name] = len(genome_seq)

        return chrom_lengths
            
        

    def write_haplotype_h5(self):
        chrom_lengths = self.get_chrom_lengths()
        
        atom = tables.Int8Atom(dflt=0)
        zlib_filter = tables.Filters(complevel=1, complib="zlib")
        
        hap_h5 = tables.open_file(self.haplotype_filename, "w")    

        chrom_haps = {}
        chrom_haps_phase = {}
        snp_index = 0

        # group haplotypes by chromosome
        # include phase information if it exists
        if self.haplotypes_phase:
            for snp, hap, phase in zip(self.snp_list, self.haplotypes, self.haplotypes_phase):
                if snp[0] in chrom_haps:
                    chrom_haps[snp[0]].append(hap)
                    chrom_haps_phase[snp[0]].append(phase)
                else:
                    chrom_haps[snp[0]] = [hap]
                    chrom_haps_phase[snp[0]] = [phase]
        else:
            for snp, hap in zip(self.snp_list, self.haplotypes):
                if snp[0] in chrom_haps:
                    chrom_haps[snp[0]].append(hap)
                else:
                    chrom_haps[snp[0]] = [hap]
            chrom_haps_phase = None
        
        for chrom, haps in list(chrom_haps.items()):
            # add haplotypes
            hap_array = np.array(haps, dtype=np.int8)
            carray = hap_h5.create_carray(hap_h5.root,
                                         chrom, atom, hap_array.shape,
                                         filters=zlib_filter)
            carray[:] = haps

            # also add phase information if it exists
            if chrom_haps_phase:
                phase_shape = (hap_array.shape[0], int(hap_array.shape[1]/2))
                phase_carray = hap_h5.create_carray(
                    hap_h5.root, "phase_"+chrom, atom, phase_shape, filters=zlib_filter
                )
                phase_carray[:] = chrom_haps_phase[chrom]
            
        self.write_hap_samples(hap_h5)

        hap_h5.close()
        
            
                
    def write_snp_index_h5(self):
        atom = tables.Int16Atom(dflt=0)
    
        zlib_filter = tables.Filters(complevel=1, complib="zlib")
        
        snp_index_h5 = tables.open_file(self.snp_index_filename, "w")    

        snp_index = 0

        chrom_arrays = {}
        chrom_lengths = self.get_chrom_lengths()
        
        for snp in self.snp_list:
            if snp[0] in chrom_arrays:
                carray = chrom_arrays[snp[0]]
            else:
                # create CArray for this chromosome
                shape = [chrom_lengths[snp[0]]]
                carray = snp_index_h5.create_carray(snp_index_h5.root,
                                                   snp[0], atom, shape,
                                                   filters=zlib_filter)
                carray[:] = -1
                chrom_arrays[snp[0]] = carray

            pos = snp[1]
            carray[pos-1] = snp_index
            snp_index += 1
            
        self.write_hap_samples(snp_index_h5)

        snp_index_h5.close()


    def write_h5_files(self):
        self.write_snp_tab_h5()
        self.write_snp_index_h5()
        self.write_haplotype_h5()
        
            




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
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_dir=test_data.snp_dir)

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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


    def test_single_two_read_two_snp_two_chrom(self):
        """Test whether having two chromosomes works, with reads
        and SNPs on both works correctly"""
                        
        test_data = Data(read1_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                       "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                        "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT",
                                         "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n" +
                                         "CCCCCCCCCCGCCCCCCCCCCCCCCCCCCC"],
                         chrom_names = ['test_chrom1', 'test_chrom2'],
                         snp_list = [['test_chrom1', 1, "A", "C"],
                                     ['test_chrom2', 3, "G", "C"]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from C to an A, and the third base of second read
        # should be switched from C to G
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]

        l = list(test_data.read1_seqs[1])
        l[2] = 'C'
        new_seq = "".join(l)

        assert lines[5] == new_seq
        assert lines[7] == test_data.read1_quals[1]
        
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


    def test_chrom_with_no_snps(self):
        """Test whether having two chromosomes works, with no
        SNPs on second chromsome"""
                        
        test_data = Data(read1_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                       "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                        "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT",
                                         "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n" +
                                         "CCCCCCCCCCGCCCCCCCCCCCCCCCCCCC"],
                         chrom_names = ['test_chrom1', 'test_chrom2'],
                         snp_list = [['test_chrom1', 1, "A", "C"]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        # write an empty file with no SNPs for chr=2
        filename = test_data.snp_dir + "/test_chrom2.snps.txt.gz"
        f = gzip.open(filename, "wt")
        f.write("");
        f.close()

        
        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False,
                                    is_sorted=False)

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from C to an A, and the second read
        # should be unchanged
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)
        assert lines[1] == new_seq
        assert lines[3] == test_data.read1_quals[0]
        
        #
        # Verify to.remap bam has 1 read
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert len(new_lines) == 1 and new_lines[0].startswith("read1")

        #
        # Verify that the keep file has 1 read
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1 and lines[0].startswith("read2")

        test_data.cleanup()
        

    def test_single_gapD_read_two_snps(self):
        """Test whether read with D in alignment works correctly"""
        #
        # Currently WASP discards reads that overlap
        # indels. We want to improve WASP to handle indels
        # correctly, but in the meantime this test just checks
        # that the read is discarded.
        #

        # create read that will align with a 'D' in CIGAR string
        test_data = Data(read1_seqs =  ["ACTGACTGAAACTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["ACTGACTGAAAAACTGACTGACTGACTGAC\n" +
                                        "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                         chrom_names = ['test_chrom1'],
                         snp_list = [['test_chrom1', 1, "A", "C"],
                                     ['test_chrom1', 15, "T", "C"]])

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        
        # pos:     123456789012345678901234567890
        # genome:  ACTGACTGAAAAACTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTTTTTTTTT
        # snps:    ^             ^
        # read:    ACTGACTGAAA--CTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTTTTTTTTT
        # read_pos:12345678901--23456789012345678
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False)


        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[12] = 'C'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[12] = 'C'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs

        # Check the new reads are named correctly
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        
        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()



    
    def test_single_gapI_read_two_snps(self):
        """Test whether read with I in alignment works correctly"""
        #
        # Currently WASP discards reads that overlap
        # indels. We want to improve WASP to handle indels
        # correctly, but in the meantime this test just checks
        # that the read is discarded.
        #

        # create read that will align with a 'D' in CIGAR string
        test_data = Data(read1_seqs =  ["ACTGACTGAAAAAAACTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTT"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["ACTGACTGAAAAACTGACTGACTGACTGAC\n" +
                                        "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT"],
                         chrom_names = ['test_chrom1'],
                         snp_list = [['test_chrom1', 1, "A", "C"],
                                     ['test_chrom1', 15, "T", "C"]])

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        
        # pos:     1234567890123--45678901234567890
        # genome:  ACTGACTGAAAAA--CTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTTTTTTTTT
        # snps:    ^               ^
        # read:    ACTGACTGAAAAAAACTGACTGACTGACTGACTTTTTTTTTTATTTTTTTTTTTTTTTTTTT
        # read_pos:1234567890123456789012345678
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False)


        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[16] = 'C'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[16] = 'C'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs

        # Check the new reads are named correctly
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        
        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()



    def test_single_one_read_one_indel(self):
        """Test whether 1 read overlapping indel works correctly"""

        test_data = Data(snp_list=[['test_chrom', 2, "A", "CCC"]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. Should be empty because only
        # read overlaps an indel
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=False, is_sorted=False)


        #
        # Verify new fastq is correct. It should only contain
        # a single read and, the first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False)

        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        
        
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=False, is_sorted=False,
                                    max_seqs=10)

        #
        # Verify new fastq is correct. There should be no reads,
        # because reads with greater than 10 allelic combinations
        # are thrown out
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=False, is_sorted=False,
                                    max_snps=10,
                                    max_seqs=1024)

        #
        # Verify new fastq is correct. There should be 1023 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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


        
    def test_single_unaligned_reads(self):
        """Simple test of whether unmapped reads are handled
        properly"""
        test_data = Data(read1_seqs=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                     "ACTAGACATACATAACACATATACCCACCC"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_dir=test_data.snp_dir)

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
        assert len(old_lines) == 2
        assert len(new_lines) == 1

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped or are discarded. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()



        

class TestCLI:
    def test_single_cli(self):
        """Make sure the command line interface
        for single-end read mapping works"""

        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        command = ['python', 'find_intersecting_snps.py',
                   test_data.bam_filename, '--snp_dir', test_data.snp_dir]
        subprocess.check_call(command)

        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
        

    def test_new_output_dir_cli(self):
        """Make sure files can be written to directory of choice"""
        out_dir = "test_output_dir"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        test_data = Data(
            prefix="test_data/test",
            output_prefix=out_dir + "/test")
        test_data.setup()
        
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        command = ['python', 'find_intersecting_snps.py',
                   "--output_dir", out_dir,
                   test_data.bam_filename, "--snp_dir",
                   test_data.snp_dir]
        subprocess.check_call(command)

        # verify that all output files written to output dir
        out_files = ["test_output_dir/test.keep.bam",
                     "test_output_dir/test.to.remap.bam",
                     "test_output_dir/test.remap.fq.gz"]
        
        for out_file in out_files:
            sys.stderr.write("%s\n"%out_file)
            assert(os.path.exists(out_file))
            os.remove(out_file)

        os.remove("test_output_dir/test.sort.bam")
        os.rmdir(out_dir)



    def test_single_cli_haplotypes_samples(self):
        ##
        ## Test with all possible combinations of haplotypes
        ## present in data
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "G", "A"]],
                         read1_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         haplotypes=[[1, 0, 1, 0],
                                     [1, 0, 0, 1]],
                         hap_samples=["samp1", "samp2"])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        command = ['python', 'find_intersecting_snps.py',
                   test_data.bam_filename, 
                   '--haplotype', test_data.haplotype_filename,
                   '--snp_tab', test_data.snp_tab_filename,
                   '--snp_index', test_data.snp_index_filename,
                   '--samples', 'samp1,samp2']
        
        subprocess.check_call(command)

        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'A'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'A'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs


        # Check the new reads are named correctly
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        

        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()


        # repeat test, but use samples file instead of
        # sample names on command line
        hap_samples = ["samp1", "samp2"]
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "G", "A"]],
                         read1_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         haplotypes=[[1, 0, 1, 0],
                                     [1, 0, 0, 1]],
                         hap_samples=hap_samples)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        samp_filename = test_data.data_dir + "/samples.txt"
        f = open(samp_filename, 'w')
        for samp in hap_samples:
            f.write("%s testdata\n" % samp)

        f.close()
        command = ['python', 'find_intersecting_snps.py',
                   test_data.bam_filename, 
                   '--haplotype', test_data.haplotype_filename,
                   '--snp_tab', test_data.snp_tab_filename,
                   '--snp_index', test_data.snp_index_filename,
                   '--samples', samp_filename]
        
        subprocess.check_call(command)
        
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'A'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'A'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs

        os.unlink(samp_filename)
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890

        snp_list = [['test_chrom', 18, "A", "C"]]
        
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 is correct.
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
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
        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # POS           123456789012345678901234567890
        
        snp_list = [['test_chrom', 18, "A", "C"]]
        
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        expect_reads = {("AACGAAAAGGAGAC", "AAGAAACAACACAA"),
                            ("ACAAAAATTTAAAA", "AAAAATAAAAAATA")}

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines2 = [x.strip() for x in f.readlines()]

        # should be same number of lines in each file
        assert len(lines1) == len(lines2)

        # number of read records (4 lines each) should match expectation
        assert len(lines2) == len(expect_reads) * 4
        
        for i in range(1, len(lines1), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            sys.stderr.write("removing pair: %s\n" % repr(read_pair))
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0
        

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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # INDEL                              ^         
        # POS           123456789012345678901234567890
        #                       40        50
        
        snp_list = [['test_chrom', 18, "A", "C"],
                    ['test_chrom', 52, "G", "GTTA"]]
        
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        # Currently reads overlapping indels are thrown out
        expect_reads = {("ACAAAAATTTAAAA", "AAAAATAAAAAATA")}

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == len(expect_reads) * 4
                           
        for i in range(1, len(lines2), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0

   
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # SNP                   ^           ^         
        # POS           123456789012345678901234567890
        #                       40        50
        
        snp_list = [['test_chrom', 18, "A", "C"],
                    ['test_chrom', 39, "T", "G"],
                    ['test_chrom', 51, "G", "T"]]
        
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)

        expect_reads = {("AACGAAAAGGAGAC", "AAGAAACAACACAA"),
                            ("AACGAAAAGGAGAC", "AAGAAACAAAACAA"),
                            ("AACGAAAAGGAGAA", "AAGAAACAAAACAA"),
                            ("ACAAAAATTTAAAA", "AAAAATAAAAAATA"),
                            ("ACAAAAATTTAAAA", "AAAAATACAAAATA"),
                            ("AAAAAAATTTAAAA", "AAAAATACAAAATA")}

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == len(expect_reads) * 4
        for i in range(1, len(lines2), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0

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





## 
## TODO: test haplotypes when list of samples provided...
##

class TestHaplotypesSingleEnd:
    """tests for single end read mapping, using known haplotypes"""

    def test_haplotypes_single_one_read_one_snp(self):
        """Simple test of whether 1 read overlapping
        1 SNP works correctly"""
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)
                                    

        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
        # repeat test, but with haoplotypes only containing non-reference allele
        #
        test_data = Data(haplotypes=[[1, 1, 1, 1]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)
                                    
        #
        # Verify new fastq is correct. The first base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
        # repeat test, but with haplotypes only containing reference allele
        #
        test_data = Data(haplotypes=[[0, 0, 0, 0]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)
                                    
        #
        # Verify new fastq is correct. There should not be any reads
        # to remap since all haplotypes match reference
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 0

        #
        # Verify to.remap bam is empty
        #
        new_lines = read_bam(test_data.bam_remap_filename)
        assert len(new_lines) == 1
        assert new_lines[0] == ''

        #
        # Verify that the keep file contains
        # the read since does not need to be remapped
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_keep_filename)
        assert old_lines == new_lines


        # Test with SNP at one before last position in read
        # since there was a bug with this situation
        test_data = Data(snp_list = [['test_chrom', 29, "A", "C"]],
                         hap_samples = ["samp1", "samp2", "samp3", "samp4"],
                         haplotypes = [[0, 1, 0, 1]])

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)
                                    

        #
        # Verify new fastq is correct. The last base of the first read
        # should be switched from a C to an A.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4

        l = list(test_data.read1_seqs[0])
        l[28] = 'C'
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


        
    def test_haplotypes_single_one_read_two_snps(self):
        """Test whether 1 read overlapping 2 SNPs works correctly"""

        ##
        ## Test with all possible combinations of haplotypes
        ## present in data
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "G", "A"]],
                         read1_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         haplotypes=[[1, 0, 1, 0],
                                     [1, 0, 0, 1]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))


        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'A'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'A'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs


        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        
        
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

        ##
        ## Test with subset of possible combinations present
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]],
                         haplotypes=[[1, 1, 1, 1, 1, 1],
                                     [1, 0, 0, 0, 1, 0]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There should be 2 reads
        # with two haplotype configurations
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        seqs = [lines[1], lines[5]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'G'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq2 = "".join(l)

        assert len(seqs) == 2
        assert new_seq1 in seqs
        assert new_seq2 in seqs

        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.2"
        assert lines[4] == "@read1.1.2.2"
                
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



    
    def test_haplotypes_wrong_sample_name(self):
        """Simple test that error o
        1 SNP works correctly"""
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    samples=["sampX", "sampY"],
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        test_data.cleanup()
                                    

        
    def test_haplotypes_single_one_read_two_snps_samples(self):
        """Test whether 1 read overlapping 2 SNPs works correctly"""

        ##
        ## Test with all possible combinations of haplotypes
        ## present in data
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "G", "A"]],
                         read1_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         haplotypes=[[1, 0, 1, 0],
                                     [1, 0, 0, 1]],
                         hap_samples=["samp1", "samp2"])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    samples=["samp1", "samp2"])

        #
        # Verify new fastq is correct. There should be 3 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 12

        seqs = [lines[1], lines[5], lines[9]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))


        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'A'
        new_seq2 = "".join(l)

        # read with both non-ref alleles
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'A'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs


        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
        
        
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

        
        ##
        ## Test with all possible combinations of haplotypes
        ## present in data, but only running on one of two samples
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "G", "A"]],
                         read1_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs=["ACTGACTGACTGACTGACTGACTGACTGACTG"],
                         haplotypes=[[1, 0, 1, 0],
                                     [1, 0, 0, 1]],
                         hap_samples=["samp1", "samp2"])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        # run using only samp2 
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    samples=["samp2"])

        # Verify new fastq is correct. There should be 2 reads
        # representing both haplotypes from samp2
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        seqs = [lines[1], lines[5]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))

        # read representing 1st haplotype from samp2
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq1 = "".join(l)

        
        # read representing 2nd haplotype from samp2
        l = list(test_data.read1_seqs[0])
        l[3] = 'A'
        new_seq2 = "".join(l)

        assert len(seqs) == 2
        assert new_seq1 in seqs
        assert new_seq2 in seqs


        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.2"
        assert lines[4] == "@read1.1.2.2"
        
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
                

        # run using only samp1 
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    samples=["samp1"])

        # Verify new fastq is correct. There should be 1 read
        # representing non-reference haplotype from samp1
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4

        seqs = [lines[1]]
        sys.stderr.write("SEQS: %s\n" % repr(seqs))

        # read representing 1st haplotype from samp1
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'A'
        new_seq1 = "".join(l)

        assert new_seq1 in seqs

        # Check the new read is named correctly
        assert lines[0] == "@read1.1.1.1"
        
        # Verify to.remap bam is the same as the input bam file.
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        
        test_data.cleanup()


        
    def test_haplotypes_one_read_ten_snps(self):
        """Test whether 1 read overlapping 10 SNPs works correctly"""

        snp_list = [['test_chrom', x, "A", "C"] for x in range(1, 11)]


        # generate all possible 1024 haplotype configurations for 10 SNPs:
        import itertools
        haplotypes = np.array([x for x in itertools.product([0, 1], repeat=10)])
        haplotypes = haplotypes.T

        test_data = Data(snp_list=snp_list,
                         haplotypes=haplotypes)

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False, is_sorted=False,
                                    max_seqs=10,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There should be no reads,
        # because reads with greater than 10 allelic combinations
        # are thrown out
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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
                                    is_paired_end=False, is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    max_snps=10,
                                    max_seqs=1024)

        #
        # Verify new fastq is correct. There should be 1023 reads
        # with all possible configurations of the two alleles, except
        # for the original configuration.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
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

        
        #
        # TEST WITH 10 SNPs but a limited number of possible haplotypes
        #

        haplotypes = np.array([[1, 1, 1, 0, 0, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0],
                               [1, 0, 1, 0, 1, 0]])

        test_data = Data(snp_list=snp_list,
                         haplotypes=haplotypes)

        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
                                       
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False, is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    max_snps=10,
                                    max_seqs=1024)

        #
        # Verify new fastq is correct. There should be 3 reads
        # (there are 4 unique haplotypes, and one of them is 
        # the original configuration)
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4*3

        # get every 4th line, which correspond to sequences starting at line 1
        seqs = [lines[x] for x in range(1, len(lines), 4)]

        # test whether haplotypes present in read
        l = list(test_data.read1_seqs[0])
        for i in range(10):
            l[0] = 'C'            
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq2 = "".join(l)

        l = list(test_data.read1_seqs[0])
        for i in range(1, 10):
            l[0] = 'C'            
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




    def test_haplotypes_bwa_softclipped(self):
        """Simple test whether reads that are softclipped
        by BWA are handled correctly."""

        # make a read that we expect to be softclipped at 5' end
        # also make a SNP that is at end of read

        snp_pos = [1, 61]
        test_data = Data(
            read1_seqs =["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTCTG"],
            read1_quals=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                         "BBBBBBBBBBBBBBBBBBBBBBBB!!!!!!!"],
            genome_seqs=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTCTG"],
            snp_list=[['test_chrom', snp_pos[0], "A", "G"],
                      ['test_chrom', snp_pos[1], "T", "G"]],
            haplotypes=[[0, 1, 0, 1],
                        [0, 1, 0, 1]])
        
        test_data.setup()
        test_data.index_genome_bwa()
        test_data.map_single_bwa()
        test_data.sam2bam()
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #  read1:  AAAAAAAAAAAAAAAAAAAAAAAAAACCTGATTTTTTTTTTATTTTTTTTTTTTTTTTTCTG
        #  qual1:  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB!!!!!!!
        #     POS: 12345678901234567890123456789012345678901234567890123456789012
        #    SNPs: G                                                           G
        #          0       10        20        30        40        50        60
        #  genome: AAAAAAAAAAAAAAAAAAAAAAAAAACCTGATTTTTTTTTTATTTTTTTTTTTTTTTTTCTG
        #  
        
        # Verify new fastq is correct. The first base of the read
        # should be switched from a A to G, but base 61 should be not
        # be changed because it is in soft-clipped region of
        # alignment
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4

        l = list(test_data.read1_seqs[0])
        l[snp_pos[0]-1] = 'G'
        new_seq = "".join(l)
        test_data.read1_quals[0]

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


        ##################################################################

        # need to test what happens when read is rev-complemented so that
        # soft-clipped region is at start of aligned region
        snp_pos = [1, 61]
        test_data = Data(
            read1_seqs =["CAGAAAAAAAAAAAAAAAAATAAAAAAAAAA"
                         "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"],
            read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" +
                         "BBBBBBBBBBBBBBBBBBBBBBBB!!!!!!!"],
            genome_seqs=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTCTG"],
            snp_list=[['test_chrom', snp_pos[0], "A", "G"],
                      ['test_chrom', snp_pos[1], "T", "G"]],
            haplotypes=[[0, 1, 0, 1],
                        [0, 1, 0, 1]])
        
        test_data.setup()
        test_data.index_genome_bwa()
        test_data.map_single_bwa()
        test_data.sam2bam()
        
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        read1_revcomp = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTATTTTTTTTTTTTTTTTTCTG"
        qual1_rev = "!!!!!!!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
        
        #   read1:  CAGAAAAAAAAAAAAAAAAATAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
        #  qual1:   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB!!!!!!!
        # rcread1:  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTATTTTTTTTTTTTTTTTTCTG
        # rqual1:   !!!!!!!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
        #    SNPs:  G                                                           G
        #     POS:  12345678901234567890123456789012345678901234567890123456789012
        #           0       10        20        30        40        50        60
        #  genome:  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTATTTTTTTTTTTTTTTTTCTG
        #

        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 4

        # expect only second SNP to remain, since first SNP should be in
        # softclipped part of read
        l = list(read1_revcomp)
        
        l[snp_pos[1]-1] = 'G'
        new_seq = "".join(l)

        assert lines[1] == new_seq
        assert lines[3] == qual1_rev

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


class TestHaplotypesPairedEnd:
	
    def test_haplotypes_paired_two_reads_two_snps(self):
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
        # SNP                            ^     ^
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]        ATAAAAAATAAAAA
        # SNP                   ^                    
        # POS           123456789012345678901234567890
        #                       40        50

        snp_list = [['test_chrom', 18, "A", "C"],
                    ['test_chrom', 24, "T", "G"],
                    ['test_chrom', 39, "T", "G"]]



        haplotypes = np.array([[1, 1, 1, 1, 0, 1, 0, 1],
                               [1, 1, 1, 1, 0, 1, 0, 1],
                               [1, 1, 1, 0, 0, 0, 0, 1]])


        test_data = Data(genome_seqs=genome_seq,
                         read1_seqs=read1_seqs,
                         read2_seqs=read2_seqs,
                         read1_quals=read1_quals,
                         read2_quals=read2_quals,
                         snp_list=snp_list,
                         haplotypes=haplotypes)
	
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_paired_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename,
                                    is_paired_end=True, is_sorted=False)

        expect_reads = {("AACGAAAAGGAGAC", "AAGAAACAACACAA"),
                            ("AAAAAAATTTAAAA", "AAAAATACAAAATA"),
                            ("ACAAAAAGTTAAAA", "AAAAATACAAAATA"),
                            ("ACAAAAAGTTAAAA", "AAAAATAAAAAATA")}

        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == len(expect_reads) * 4

        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == len(expect_reads) * 4
        for i in range(1, len(lines2), 4):
            read_pair = (lines1[i], lines2[i])
            assert read_pair in expect_reads
            expect_reads.remove(read_pair)

        assert len(expect_reads) == 0

        #
        # Verify that the keep file is empty since only
        # read needs to be remapped. Note that the
        # read_bam still gives back one empty line.
        #
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()


    def test_paired_unaligned_reads(self):
        """Simple test of whether unmapped reads are handled
        properly"""


        # test when only one half of read maps, other half is present as
        # unaligned read

        # create 4 read pairs
        # read1: both ends map
        # read2: only first 1/2 maps
        # read3: neither 1/2 maps
        # read4: only second 1/2 maps
        
        test_data = Data(genome_seqs=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
                                      "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"],
                         read1_seqs=["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                     "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                     "ACTAGACATACATAACACATATACCCACCC",
                                     "ACTAGACATACATAACACATATACCCACCC"],
                         read2_seqs=["AAGAAACAACACAAAAAAATAAAAAATAAA",
                                     "ACTAGACATACATAACACATATACCCACCC",
                                     "ACTAGACATACATAACACATATACCCACCC",
                                     "AAGAAACAACACAAAAAAATAAAAAATAAA"],
                         read1_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         read2_quals=["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                      "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"])

        
                         
                         
        test_data.setup()
        # test_data.index_genome_bowtie2()
        # test_data.map_paired_bowtie2()
        test_data.index_genome_bwa()
        test_data.map_paired_bwa()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=True,
                                    is_sorted=False,
                                    snp_dir=test_data.snp_dir)


        #
        # Verify fastq1 and fastq2 have appropriate read pairs
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines1 = [x.strip() for x in f.readlines()]
        assert len(lines1) == 4

        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines2 = [x.strip() for x in f.readlines()]
        assert len(lines2) == 4
        
        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq = "".join(l)

        assert lines1[1] == new_seq
        assert lines1[3] == test_data.read1_quals[0]

        assert lines2[1] == test_data.read2_seqs[0]
        assert lines2[3] == test_data.read2_quals[0]

        test_data.cleanup()


class TestFiltering:
    """tests that bad reads are filtered correctly"""

    def test_remap_secondary_alignments(self):
        """Test that secondary alignments don't appear in remap files"""
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_dir=test_data.snp_dir)
        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # now, let's alter the flag of the read and
        # rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            read.flag = 256
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the remap file shouldn't have any reads
        lines = read_bam(test_data.bam_remap_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_keep_secondary_alignments(self):
        """Test that secondary alignments don't appear in keep files"""
        test_data = Data(read1_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                       "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                        "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT",
                                         "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n" +
                                         "CCCCCCCCCCGCCCCCCCCCCCCCCCCCCC"],
                         chrom_names = ['test_chrom1', 'test_chrom2'],
                         snp_list = [['test_chrom1', 1, "A", "C"]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        # write an empty file with no SNPs for chr=2
        filename = test_data.snp_dir + "/test_chrom2.snps.txt.gz"
        f = gzip.open(filename, "wt")
        f.write("");
        f.close()

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # verify that the keep file has one read
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1 and lines[0].startswith("read2")

        # now, let's alter the flag of the read and
        # rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            read.flag = 256
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the keep file shouldn't have any reads
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_remap_secondary_alignments_pe(self):
        """Test that secondary alignments don't appear in remap files
        for paired end reads"""
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890
        # first, initiate an empty snp list so that reads are kept
        snp_list = [['test_chrom', 18, "A", "C"]]
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)
        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # now, let's alter the flag of one of the reads in each pair
        # and rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            # make first in read pair a secondary alignment
            if read.qname == 'read1' and read.flag == 99:
                read.flag = 355
            # make the second in the read pair a secondary alignment
            elif read.qname == 'read2' and read.flag == 147:
                read.flag = 403
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the remap file shouldn't have any reads
        lines = read_bam(test_data.bam_remap_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_keep_secondary_alignments_pe(self):
        """Test that secondary alignments don't appear in keep files
        for paired end reads"""
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890
        # first, initiate an empty snp list so that reads are kept
        snp_list = []
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 4

        # now, let's alter the flag of one of the reads in each pair
        # and rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            # make first in read pair a secondary alignment
            if read.qname == 'read1' and read.flag == 99:
                read.flag = 355
            # make the second in the read pair a secondary alignment
            elif read.qname == 'read2' and read.flag == 147:
                read.flag = 403
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the keep file shouldn't have any reads
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_remap_supplementary_alignments(self):
        """Test that supplementary alignments don't appear in remap files"""
        test_data = Data()
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_dir=test_data.snp_dir)
        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # now, let's alter the flag of the read and
        # rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            read.flag = 2048
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the remap file shouldn't have any reads
        lines = read_bam(test_data.bam_remap_filename)
        assert len(lines) == 1
        assert lines[0] == ''

    def test_keep_supplementary_alignments(self):
        """Test that supplementary alignments don't appear in keep files"""
        test_data = Data(read1_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                       "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"],
                         read1_quals = ["BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                                        "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"],
                         genome_seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" +
                                         "TTTTTTTTTTATTTTTTTTTTTTTTTTTTT",
                                         "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n" +
                                         "CCCCCCCCCCGCCCCCCCCCCCCCCCCCCC"],
                         chrom_names = ['test_chrom1', 'test_chrom2'],
                         snp_list = [['test_chrom1', 1, "A", "C"]])
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()
        # write an empty file with no SNPs for chr=2
        filename = test_data.snp_dir + "/test_chrom2.snps.txt.gz"
        f = gzip.open(filename, "wt")
        f.write("");
        f.close()

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # verify that the keep file has one read
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1 and lines[0].startswith("read2")

        # now, let's alter the flag of the read and rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            read.flag = 2048
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the keep file shouldn't have any reads
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_remap_supplementary_alignments_pe(self):
        """Test that supplementary alignments don't appear in remap file
        for paired end reads"""
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890
        # first, initiate an empty snp list so that reads are kept
        snp_list = [['test_chrom', 18, "A", "C"]]
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)
        #
        # Verify to.remap bam is the same as the input bam file.
        #
        old_lines = read_bam(test_data.bam_filename)
        new_lines = read_bam(test_data.bam_remap_filename)
        assert old_lines == new_lines

        # now, let's alter the flag of one of the reads in each pair
        # and rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            # make first in read pair a secondary alignment
            if read.qname == 'read1' and read.flag == 99:
                read.flag = 2147
            # make the second in the read pair a secondary alignment
            elif read.qname == 'read2' and read.flag == 147:
                read.flag = 2195
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)

        # now, the remap file shouldn't have any reads
        lines = read_bam(test_data.bam_remap_filename)
        assert len(lines) == 1
        assert lines[0] == ''

        test_data.cleanup()

    def test_keep_supplementary_alignments_pe(self):
        """Test that supplementary alignments don't appear in keep file
        for paired end reads"""
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
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read2[0]                      AACACAACAAAGAA
        # read2[1]                      AACACAACAAAGAA
        # POS           123456789012345678901234567890
        # first, initiate an empty snp list so that reads are kept
        snp_list = []
        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir, 
                                    is_paired_end=True, is_sorted=False)
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 4

        # now, let's alter the flag of the read and rerun find_intersecting_snps
        bam = pysam.Samfile(test_data.bam_filename, "rb")
        new_bam_file = test_data.bam_filename[:-4] + "2.bam"
        new_bam = pysam.AlignmentFile(new_bam_file, "wb", template=bam)
        for read in bam:
            # make first in read pair a secondary alignment
            if read.qname == 'read1' and read.flag == 99:
                read.flag = 2147
            # make the second in the read pair a secondary alignment
            elif read.qname == 'read2' and read.flag == 147:
                read.flag = 2195
            new_bam.write(read)
        bam.close()
        new_bam.close()
        # replace old bam file
        os.rename(new_bam_file, test_data.bam_filename)

        find_intersecting_snps.main(test_data.bam_filename,
                                    snp_dir=test_data.snp_dir, is_paired_end=False,
                                    is_sorted=False)


        # now, the keep file shouldn't have any reads
        lines = read_bam(test_data.bam_keep_filename)
        assert len(lines) == 1
        assert lines[0] == ''
        
        
class TestOverlappingPEReads:
    """tests that bad reads are filtered correctly"""

    def test_overlap_paired_two_reads_one_snp(self):
        """Simple test of whether 2 pairs of reads with both pairs
        overlapping a SNP works correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "TTTATTTTTTATTT"]
        read2_seqs = ["TTTTAAATTTTTTT",  # AAAAAAATTTAAAA
                      "ACAACACAAAAAAA"]  # TTTTTTTGTGTTGT

        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read2[0]                      AAAAAAATTTAAAA
        # SNP                            ^
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read1[1]      TTTATTTTTTATTT
        # read2[1]                 TTTTTTTGTGTTGT
        # SNP                       ^
        # POS           123456789012345678901234567890

        snp_list = [['test_chrom', 18, "A", "C"],
                    ['test_chrom', 43, "T", "G"]]

        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 is correct.
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert(len(lines) == 8)

        l = list(test_data.read1_seqs[0])
        # last base of first read should be changed from A to C
        l[13] = 'C'
        new_seq = "".join(l)
        assert(lines[1] == new_seq)
        assert(lines[3] == test_data.read1_quals[0])

        l = list(test_data.read1_seqs[1])
        # second to last base of second read should be changed from T to G
        l[12] = 'G'
        new_seq = "".join(l)
        assert(lines[5] == new_seq)
        assert(lines[7] == test_data.read1_quals[1])

        #
        # verify fastq2 is correct
        #
        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert(len(lines) == 8)

        l = list(test_data.read2_seqs[0])
        # second to last base of first read should be changed from T to G (since G is the complement of C)
        l[12] = 'G'
        new_seq = "".join(l)
        assert(lines[1] == new_seq)
        assert(lines[3] == test_data.read2_quals[0])

        l = list(test_data.read2_seqs[1])
        # second to last base of second read should be changed from A to C (since C is the complement of G)
        l[12] = 'C'
        new_seq = "".join(l)
        assert(lines[5] == new_seq)
        assert(lines[7] == test_data.read2_quals[1])

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

    def test_overlap_paired_two_reads_two_snps(self):
        """Test of whether 2 pairs of reads with both pairs
        overlapping 2 SNPs works correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA",
                      "TTTATTTTTTATTT"]
        read2_seqs = ["TTTTAAATTTTTTT",  # AAAAAAATTTAAAA
                      "ACAACACAAAAAAA"]  # TTTTTTTGTGTTGT

        read1_quals = ["B" * len(read1_seqs[0]),
                       "C" * len(read1_seqs[1])]
        read2_quals = ["D" * len(read2_seqs[0]),
                       "E" * len(read2_seqs[1])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read2[0]                      AAAAAAATTTAAAA
        # SNPs                          ^^
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]
        # read1[1]      TTTATTTTTTATTT
        # read2[1]                 TTTTTTTGTGTTGT
        # SNPs                     ^^
        # POS           123456789012345678901234567890

        snp_list = [['test_chrom', 17, "A", "T"],
                    ['test_chrom', 18, "A", "C"],
                    ['test_chrom', 42, "T", "A"],
                    ['test_chrom', 43, "T", "G"]]

        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 and fastq2 is correct.
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            seqs1 = [x.strip() for x in f.readlines()][1::4]
            # there are two pairs of reads and two snps (each with two alleles)
            # leading to 2^2 combinations of alleles. however, one of those
            # combinations is the original pair (hence why we subtract by one)
            assert(len(set(seqs1)) == 2*(2**2-1))
        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            seqs2 = [x.strip() for x in f.readlines()][1::4]
            # there are two pairs of reads and two snps (each with two alleles)
            # leading to 2^2 combinations of alleles. however, one of those
            # combinations is the original pair (hence why we subtract by one)
            assert(len(set(seqs2)) == 2*(2**2-1))
        pairs = list(zip(seqs1, seqs2))

        expect_reads1 = [
            # only last base of first read should be changed from A to C
            "AACGAAAAGGAGAC",
            # only second to last base of first read should be changed from A
            # to T
            "AACGAAAAGGAGTA",
            # last base of first read should be changed from A to C
            # AND second to last base of first read should be changed from A
            # to T
            "AACGAAAAGGAGTC",
            # second to last base of second read should be changed from T to G
            "TTTATTTTTTATGT",
            # only third to last base of second read should be changed from T
            # to A
            "TTTATTTTTTAATT",
            # second to last base of second read should be changed from T to G
            # AND third to last base of second read should be changed from T
            # to A
            "TTTATTTTTTAAGT"
        ]
        expect_reads2 = [
            # only second to last base of first read should be changed from T
            # to G (since it is the complement of C)
            "TTTTAAATTTTTGT",
            # only last base of first read should be changed from T to A (since
            # it is the complement of T)
            "TTTTAAATTTTTTA",
            # last base of first read should be changed from T to A (since it
            # is the complement of T)
            # AND second to last base of first read should be changed from T
            # to G (since it is the complement of C)
            "TTTTAAATTTTTGA",
            # second to last base of second read should be changed from A to C
            # (since it is the complement of G)
            "ACAACACAAAAACA",
            # last base of second read should be changed from A to T (since it
            # is the complement of C)
            "ACAACACAAAAAAT",
            # last base of second read should be changed from A to T (since it
            # is the complement of A)
            # AND second to last base of second read should be changed from A
            # to C (since it is the complement of G)
            "ACAACACAAAAACT"
        ]
        expect_pairs = list(zip(expect_reads1, expect_reads2))

        # are all of the pairs there?
        assert(len(expect_pairs) == len(pairs))
        for pair in expect_pairs:
            assert(pair in pairs)

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

    def test_overlap_paired_one_read_one_snp_one_snp(self):
        """Test of whether 1 pair of reads with both pairs
        overlapping the same SNP and one other SNP works correctly"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA"]
        read2_seqs = ["TTTTAAATTTTTTT"]  # AAAAAAATTTAAAA

        read1_quals = ["B" * len(read1_seqs[0])]
        read2_quals = ["D" * len(read2_seqs[0])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read2[0]                      AAAAAAATTTAAAA
        # SNPs                   ^       ^^
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]

        snp_list = [['test_chrom', 10, "A", "T"],
                    ['test_chrom', 18, "A", "C"],
                    ['test_chrom', 19, "A", "C"]]

        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 and fastq2 are correct.
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            seqs1 = [x.strip() for x in f.readlines()][1::4]
            # there is one pair of reads and two snps (each with two alleles),
            # leading to 2^2 combinations of alleles. one of those
            # combinations is the original pair (hence why we subtract by one)
            assert(len(set(seqs1)) == 2**2)
        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            seqs2 = [x.strip() for x in f.readlines()][1::4]
            # there is one pair of reads and two snps (each with two alleles),
            # leading to 2^2 combinations of alleles. one of those
            # combinations is the original pair (hence why we subtract by one)
            assert(len(set(seqs2)) == 2**2)
        pairs = list(zip(seqs1, seqs2))

        # which pairs of reads do we expect to see?
        # the reads below are grouped by whether they contain the shared snp
        expect_reads1 = [
            # reads with the reference allele on the shared snp
            set([
                # original read
                'AACGAAAAGGAGAA',
                # ninth to last base of read should be changed from A to T
                'AACGATAAGGAGAA'
            ]),
            # reads with the alternate allele on the shared snp
            set([
                # last base of read should be changed from A to C
                'AACGAAAAGGAGAC',
                # last base of read should be changed from A to C
                # AND ninth to last base of read should be changed from A to T
                'AACGATAAGGAGAC'
            ])
        ]
        expect_reads2 = [
            # reads with the reference allele on the shared snp
            set([
                # original read
                'TTTTAAATTTTTTT',
                # third to last base of read should be changed from T to G
                'TTTTAAATTTTGTT'
            ]),
            # reads with the alternate allele on the shared snp
            set([
                # second to last base of read should be changed from T to G
                'TTTTAAATTTTTGT',
                # second to last base of read should be changed from T to G
                # AND third to last base of read should be changed from T to G
                'TTTTAAATTTTGGT'
            ])
        ]
        from itertools import product
        # only compute combinations of reads from the same group
        expect_pairs = list(product(expect_reads1[0], expect_reads2[0])) + list(product(expect_reads1[1], expect_reads2[1]))
        # remove the original pair
        expect_pairs.remove((test_data.read1_seqs[0], test_data.read2_seqs[0]))

        # are all of the pairs there?
        assert(len(expect_pairs) == len(pairs))
        for pair in expect_pairs:
            assert(pair in pairs)

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

    def test_overlap_paired_one_read_one_discordant_snp(self):
        """Test of 1 pair of reads with both pairs
        overlapping the same SNP but with each pair having
        a different allele"""
        test_data = Data()

        read1_seqs = ["AACGAAAAGGAGAA"]
        read2_seqs = ["TTTTAAATTTTTGT"]  # ACAAAAATTTAAAA

        read1_quals = ["B" * len(read1_seqs[0])]
        read2_quals = ["D" * len(read2_seqs[0])]

        # POS           123456789012345678901234567890
        # read1[0]          AACGAAAAGGAGAA
        # read2[0]                      ACAAAAATTTAAAA
        # SNPs                           ^
        genome_seq =  ["AAAAAACGAAAAGGAGAAAAAAATTTAAAA\n"
                       "TTTATTTTTTATTTTTTTGTGTTGTTTCTT"]

        snp_list = [['test_chrom', 18, "A", "C"]]

        test_data = Data(genome_seqs=genome_seq,
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
                                    snp_dir=test_data.snp_dir,
                                    is_paired_end=True, is_sorted=False)

        #
        # Verify new fastq1 and fastq2 are correct.
        #
        with gzip.open(test_data.fastq1_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
            # the reads should have been tossed out
            assert(len(lines) == 0)
        with gzip.open(test_data.fastq2_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
            # the reads should have been tossed out
            assert(len(lines) == 0)

        #
        # Verify to.remap bam is empty since all reads were
        # tossed. Note that read_bam still gives back one empty line.
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

        test_data.cleanup()


class TestUnphasedHaplotypes:
    """tests that haplotypes are handled correctly when not all phased"""

    def test_haplotypes_phase_single_one_read_two_snps(self):
        """Test whether 1 read overlapping 2 SNPs works correctly"""

        ##
        ## Test with subset of possible combinations present but without phase
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]],
                         haplotypes=[[1, 1, 0, 1],
                                     [1, 0, 0, 1]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There should be 2 reads
        # with two haplotype configurations
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 8

        seqs = [lines[1], lines[5]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'G'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq2 = "".join(l)

        assert len(seqs) == 2
        assert new_seq1 in seqs
        assert new_seq2 in seqs

        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.2"
        assert lines[4] == "@read1.1.2.2"
                
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

        ##
        ## Test with subset of possible combinations present but with phase now
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]],
                         haplotypes=[[1, 1, 0, 1],
                                     [1, 0, 0, 1]],
                         haplotypes_phase=[[1, 0],
                                           [0, 1]])
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There would usually be 2 reads
        # with two haplotype configurations, but now there should be one more
        # read with a new haplotype configuration because of the phase info we
        # added.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 3*4

        seqs = [lines[1], lines[5], lines[9]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'G'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq2 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'G'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs

        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
                
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

        ##
        ## Test with subset of possible combinations present but without
        ## providing any phase information. The results should be the same as
        ## the previous test.
        ##
        test_data = Data(snp_list = [['test_chrom', 1, "A", "C"],
                                     ['test_chrom', 4, "A", "G"]],
                         haplotypes=[[1, 1, 0, 1],
                                     [1, 0, 0, 1]],
                         haplotypes_phase=None)
        
        test_data.setup()
        test_data.index_genome_bowtie2()
        test_data.map_single_bowtie2()
        test_data.sam2bam()

        find_intersecting_snps.main(test_data.bam_filename,
                                    is_paired_end=False,
                                    is_sorted=False,
                                    snp_tab_filename=test_data.snp_tab_filename,
                                    snp_index_filename=test_data.snp_index_filename,
                                    haplotype_filename=test_data.haplotype_filename)

        #
        # Verify new fastq is correct. There would usually be 2 reads
        # with two haplotype configurations, but now there should be one more
        # read with a new haplotype configuration because of all haplotypes
        # should be assumed unphased.
        #
        with gzip.open(test_data.fastq_remap_filename, "rt") as f:
            lines = [x.strip() for x in f.readlines()]
        assert len(lines) == 3*4

        seqs = [lines[1], lines[5], lines[9]]

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        l[3] = 'G'
        new_seq1 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[0] = 'C'
        new_seq2 = "".join(l)

        l = list(test_data.read1_seqs[0])
        l[3] = 'G'
        new_seq3 = "".join(l)

        assert len(seqs) == 3
        assert new_seq1 in seqs
        assert new_seq2 in seqs
        assert new_seq3 in seqs

        #
        # Check the new reads are named correctly
        #
        assert lines[0] == "@read1.1.1.3"
        assert lines[4] == "@read1.1.2.3"
        assert lines[8] == "@read1.1.3.3"
                
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

