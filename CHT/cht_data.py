import sys
import gzip
import os
import numpy as np

import util


class TestSNP:
    def __init__(self, name, geno_hap1, geno_hap2, AS_target_ref, AS_target_alt,
                 hetps, totals, counts):
        self.name = name
        self.geno_hap1 = geno_hap1
        self.geno_hap2 = geno_hap2
        self.AS_target_ref = AS_target_ref
        self.AS_target_alt = AS_target_alt
        self.hetps = hetps
        self.totals = totals
        self.counts = counts


    def is_het(self):
        """returns True if the test SNP is heterozygous"""
        return self.geno_hap1 != self.geno_hap2

    def is_homo_ref(self):
        """Returns True if test SNP is homozygous for reference allele"""
        return self.geno_hap1 == 0 and self.geno_hap2 == 0

    def is_homo_alt(self):
        """Returns True if test SNP is homozygous for non-reference allele"""
        return self.geno_hap1 == 1 and self.geno_hap2 == 1


dup_snp_warn = True


def parse_test_snp(snpinfo, shuffle=False):
    global dup_snp_warn
    snp_id = snpinfo[2]

    tot = 0 if snpinfo[16] == "NA" else float(snpinfo[16])

    if snpinfo[6] == "NA":
        geno_hap1 = 0
        geno_hap2 = 0
    else:
        geno_hap1 = int(snpinfo[6].strip().split("|")[0])
        geno_hap2 = int(snpinfo[6].strip().split("|")[1])

    count = 0 if snpinfo[15] == "NA" else int(snpinfo[15])

    if snpinfo[9].strip() == "NA" or geno_hap1 == geno_hap2:
        # SNP is homozygous, so there is no AS info
        return TestSNP(snp_id, geno_hap1, geno_hap2, [], [], [], tot, count)
    else:
        # positions of target SNPs
        snp_locs = np.array([int(y.strip()) for y in snpinfo[9].split(';')])

        # counts of reads that match reference overlapping linked 'target' SNPs
        snp_as_ref = np.array([int(y) for y in snpinfo[12].split(';')])

        # counts of reads that match alternate allele
        snp_as_alt = np.array([int(y) for y in snpinfo[13].split(';')])

        # heterozygote probabilities
        snp_hetps = np.array([np.float64(y.strip())
                              for y in snpinfo[10].split(';')])

        # linkage probabilities, not currently used
        snp_linkageps = np.array([np.float64(y.strip())
                                  for y in snpinfo[11].split(';')])


        # same SNP should not be provided multiple times, this
        # can create problems with combined test. Warn and filter
        # duplicate SNPs
        uniq_loc, uniq_idx = np.unique(snp_locs, return_index=True)

        if dup_snp_warn and uniq_loc.shape[0] != snp_locs.shape[0]:
            sys.stderr.write("WARNING: discarding SNPs that are repeated "
                                     "multiple times in same line\n")
            # only warn once
            dup_snp_warn = False

        snp_as_ref = snp_as_ref[uniq_idx]
        snp_as_alt = snp_as_alt[uniq_idx]
        snp_hetps = snp_hetps[uniq_idx]

        # linkage probabilities currently not used
        snp_linkageps = snp_linkageps[uniq_idx]

        if shuffle:
            # permute allele-specific read counts by flipping them randomly at
            # each SNP
            for y in range(len(snp_as_ref)):
                if random.randint(0, 1) == 1:
                    temp = snp_as_ref[y]
                    snp_as_ref[y] = snp_as_alt[y]
                    snp_as_alt[y] = temp

        return TestSNP(snp_id, geno_hap1, geno_hap2, snp_as_ref,
                       snp_as_alt, snp_hetps, tot, count)


    


def open_input_files(in_filename):
    if not os.path.exists(in_filename) or not os.path.isfile(in_filename):
        raise IOError("input file %s does not exist or is not a "
                      "regular file\n" % in_filename)

    # read file that contains list of input files
    in_file = open(in_filename)

    infiles = []
    for line in in_file:
        # open each input file and read first line
        filename = line.rstrip()
        sys.stderr.write(" " + filename + "\n")
        if (not filename) or (not os.path.exists(filename)) or \
           (not os.path.isfile(filename)):
            sys.stderr.write("input file '%s' does not exist or is not a "
                             "regular file\n" % in_file)
            exit(2)
        if util.is_gzipped(filename):
            f = gzip.open(filename, "rt")
        else:
            f = open(filename)

        # skip header
        f.readline()

        infiles.append(f)
    in_file.close()

    if len(infiles) == 0:
        sys.stderr.write("no input files specified in file '%s'\n" % in_filename)
        exit(2)

    return infiles





    
def read_count_matrices(input_filename, shuffle=False, skip=0,
                        min_counts=0, min_as_counts=0, sample=0):
    """Given an input file that contains paths to input files for all individuals, and returns 
    matrix of observed read counts, and matrix of expected read counts
    """
    infiles = open_input_files(input_filename)

    is_finished = False
    count_matrix = []
    expected_matrix = []
    line_num = 0
    skip_num = 0

    while not is_finished:
        is_comment = False
        line_num += 1
        count_line = []
        expected_line = []
        num_as = 0

        for i in range(len(infiles)):
            # read next row from this input file
            line = infiles[i].readline().strip()

            if line.startswith("#") or line.startswith("CHROM"):
                # skip comment lines and header line
                is_comment = True
            elif line:
                if is_finished:
                    raise IOError("All input files should have same number of lines. "
                                  "LINE %d is present in file %s, but not in all input files\n"
                                  % (line_num, infiles[i].name))
                if is_comment:
                    raise IOError("Comment and header lines should be consistent accross "
                                  "all input files. LINE %d is comment or header line in some input files "
                                  "but not in file %s" % (line_num, infiles[i].name))

                # parse test SNP and associated info from input file row
                new_snp = parse_test_snp(line.split(), shuffle=shuffle)
                if new_snp.is_het():
                    num_as += np.sum(new_snp.AS_target_ref) + \
                      np.sum(new_snp.AS_target_alt)

                count_line.append(new_snp.counts)
                expected_line.append(new_snp.totals)

            else:
                # out of lines from at least one file, assume we are finished
                is_finished = True

        if not is_finished and not is_comment:
            if skip_num < skip:
                # skip this row
                skip_num += 1
            else:
                if(sum(count_line) >= min_counts and num_as >= min_as_counts):
                    # this line exceeded minimum number of read counts and AS counts
                    count_matrix.append(count_line)
                    expected_matrix.append(expected_line)
                skip_num = 0

    count_matrix = np.array(count_matrix, dtype=int)
    expected_matrix = np.array(expected_matrix, dtype=np.float64)

    sys.stderr.write("count_matrix dimension: %s\n" % str(count_matrix.shape))
    sys.stderr.write("expect_matrix dimension: %s\n" % str(expected_matrix.shape))

    nrow = count_matrix.shape[0]
    if (sample > 0) and (sample < count_matrix.shape):
        # randomly sample subset of rows without replacement
        sys.stderr.write("randomly sampling %d target regions\n" % sample)
        samp_index = np.arange(nrow)
        np.random.shuffle(samp_index)
        samp_index = samp_index[:sample]
        count_matrix = count_matrix[samp_index,]
        expected_matrix = expected_matrix[samp_index,]

        sys.stderr.write("new count_matrix dimension: %s\n" % str(count_matrix.shape))
        sys.stderr.write("new expect_matrix dimension: %s\n" % str(expected_matrix.shape))

    return count_matrix, expected_matrix
