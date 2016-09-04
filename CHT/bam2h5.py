#!/bin/env python
#
# Copyright 2013-2014 Graham McVicker and Bryce van de Geijn
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
#
"""
This program reads BAM files and counts the number of reads that match
the alternate and reference allele at every SNP position in the provided
SNP HDF5 data files. The read counts are stored in specified HDF5 output
files.

Additionally counts of all reads are stored in another track (at the
left-most position of the reads).

This program does not perform filtering of reads based on mappability.
It is assumed that the inpute BAM files are filtered appropriately prior to
calling this script.

Reads that overlap known indels are not included in allele-specific
counts.


usage: bam2h5.py OPTIONS BAM_FILE1 [BAM_FILE2 ...]

BAM Files:
     Aligned reads are read from one or more BAM files. The provided
     BAM files must be sorted and indexed.

Input Options:
     --chrom CHROM_TXT_FILE [required]
	   Path to chromInfo.txt file (may be gzipped) with list of
	   chromosomes for the relevant genome assembly. Each line
	   in file should contain tab-separated chromosome name and
	   chromosome length (in basepairs). chromInfo.txt files can
	   be downloaded from the UCSC genome browser. For example,
	   a chromInfo.txt.gz file for hg19 can be downloaded from
	   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

     --snp_index SNP_INDEX_H5_FILE [required]
       Path to HDF5 file containing SNP index. The SNP index is
       used to convert the genomic position of a SNP to its
       corresponding row in the haplotype and snp_tab
       HDF5 files.

     --snp_tab SNP_TABLE_H5_FILE [required]
       Path to HDF5 file to read SNP information from. Each row of SNP
       table contains SNP name (rs_id), position, allele1, allele2.

     --haplotype HAPLOTYPE_H5_FILE [optional]
       Path to HDF5 file to read phased haplotypes from.
       If supplied, when read overlaps multiple SNPs counts are randomly
       assigned to ONE of the overlapping HETEROZYGOUS SNPs; if not supplied
       counts are randomly assigned to ONE of overlapping SNPs (regardless of
       their genotype).

     --samples SAMPLES_TXT_FILE [optional]
       Path to text file containing a list of individual identifiers. The
       ordering of individuals must be consistent with the haplotype
       file. The samples file is assumed to have one identifier per line
       in the first column (other columns are ignored).

     --individual INDIVIDUAL [optional]
       Identifier for individual, used to determine which
       SNPs are heterozygous. Must be provided
       if --haplotype argument is provided and must match one of the
       individuals in the file provided with --samples argument.

Output Options:
     --data_type uint8|uint16
       Data type of stored counts; uint8 takes up less disk
       space but has a maximum value of 255 (default=uint8).

     --ref_as_counts REF_AS_COUNT_H5_FILE [required]
       Path to HDF5 file to write counts of reads that match reference allele.
       Allele-specific counts are stored at the position of the SNP.

     --alt_as_counts ALT_AS_COUNT_H5_FILE [required]
       Path to HDF5 file to write counts of reads that match alternate allele.
       Allele-specific counts are stored at the position of the SNP.

     --other_as_counts OTHER_AS_COUNT_H5_FILE [required]
       Path to HDF5 file to write counts of reads that match neither reference
       nor alternate allele. Allele-specific counts are stored at the position
       of the SNP.

     --read_counts READ_COUNT_H5_FILE [required]
       Path to HDF5 file to write counts of all reads, regardless of whether
       they overlap a SNP. Read counts are stored at the left-most position
       of the mapped read.
"""

import sys
import os
import gzip

import tables
import argparse
import numpy as np

import pysam

import chromosome
import chromstat


# codes used by pysam for aligned read CIGAR strings
BAM_CMATCH     = 0 # M
BAM_CINS       = 1 # I
BAM_CDEL       = 2 # D
BAM_CREF_SKIP  = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD       = 6 # P
BAM_CEQUAL     = 7 # =
BAM_CDIFF      = 8 # X

BAM_CIGAR_DICT = {0 : "M",
                  1 : "I",
                  2 : "D",
                  3 : "N",
                  4 : "S",
                  5 : "H",
                  6 : "P",
                  7 : "=",
                  8 : "X"}


SNP_UNDEF = -1
MAX_UINT8_COUNT = 255
MAX_UINT16_COUNT = 65535



def create_carray(h5f, chrom, data_type):
    if data_type == "uint8":
        atom = tables.UInt8Atom(dflt=0)
    elif data_type == "uint16":
        atom = tables.UInt16Atom(dflt=0)
    else:
        raise NotImplementedError("unsupported datatype %s" % data_type)

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    # create CArray for this chromosome
    shape = [chrom.length]
    carray = h5f.createCArray(h5f.root, chrom.name,
                              atom, shape, filters=zlib_filter)

    return carray



def get_carray(h5f, chrom):
    return h5f.getNode("/%s" % chrom)




def is_indel(snp):
    if (len(snp['allele1']) != 1) or (len(snp['allele2'])) != 1:
        return True


def dump_read(f, read):
    cigar_str = " ".join(["%s:%d" % (BAM_CIGAR_DICT[c[0]], c[1])
                          for c in read.cigar])

    f.write("pos: %d\n"
            "aend: %d\n"
            "alen (len of aligned portion of read on genome): %d\n"
            "qstart: %d\n"
            "qend: %d\n"
            "qlen (len of aligned qry seq): %d\n"
            "rlen (read len): %d\n"
            "tlen (insert size): %d\n"
            "cigar: %s\n"
            "seq: %s\n"
            % (read.pos, read.aend, read.alen, read.qstart, read.qend,
               read.qlen, read.rlen, read.tlen, cigar_str, read.seq))




def get_sam_iter(samfile, chrom):
    try:
        sam_iter = samfile.fetch(reference=chrom.name,
                                 start=1, end=chrom.length)
    except ValueError as ve:
        sys.stderr.write("%s\n" % str(ve))
        # could not find chromosome, try stripping leading 'chr'
        # E.g. for drosophila, sometimes 'chr2L' is used but
        # othertimes just '2L' is used. Annoying!
        chrom_name = chrom.name.replace("chr", "")
        sys.stderr.write("WARNING: %s does not exist in BAM file, "
                         "trying %s instead\n" % (chrom.name, chrom_name))

        try:
            sam_iter = samfile.fetch(reference=chrom_name,
                                     start=1, end=chrom.length)
        except ValueError:
            # fetch can fail because chromosome is missing or because
            # BAM has not been indexed
            sys.stderr.write("WARNING: %s does not exist in BAM file, "
                             "or BAM file has not been sorted and indexed.\n"
                             "         Use 'samtools sort' and 'samtools index' to "
                             "index BAM files before running bam2h5.py.\n"
                             "         Skipping chromosome %s.\n" %
                             (chrom.name, chrom.name))
            sam_iter = iter([])

    return sam_iter




def choose_overlap_snp(read, snp_tab, snp_index_array, hap_tab, ind_idx):
    """Picks out a single SNP from those that the read overlaps.
    Returns a tuple containing 4 elements: [0] the index of the SNP in
    the SNP table, [1] the offset into the read sequence, [2] flag
    indicating whether the read was 'split' (i.e. was a spliced
    read), [3] flag indicating whether read overlaps known indel.
    If there are no overlapping SNPs or the read cannot be processed,
    (None, None, is_split, overlap_indel) is returned instead.
    """
    read_offsets = []
    snp_idx = []

    read_start_idx = 0
    genome_start_idx = read.pos

    n_match_segments = 0
    is_split = False
    overlap_indel = False

    for cig in read.cigar:
        op = cig[0]
        op_len = cig[1]

        if op == BAM_CMATCH:
            # this is a block of match/mismatch in read alignment
            read_end = read_start_idx + op_len
            genome_end = genome_start_idx + op_len

            # get offsets of any SNPs that this read overlaps
            idx = snp_index_array[genome_start_idx:genome_end]
            is_def = np.where(idx != SNP_UNDEF)[0]
            read_offsets.extend(read_start_idx + is_def)
            snp_idx.extend(idx[is_def])

            read_start_idx = read_end
            genome_start_idx = genome_end

            n_match_segments += 1
        elif op == BAM_CREF_SKIP:
            # spliced read, skip over this region of genome
            genome_start_idx += op_len
            is_split = True
        elif op == BAM_CSOFT_CLIP:
            # end of read is soft-clipped, which means it is
            # present in read, but not used in alignment
            read_start_idx += op_len
        elif op == BAM_CHARD_CLIP:
            # end of read is hard-clipped, so not present
            # in read and not used in alignment
            pass
        else:
            sys.stderr.write("skipping because contains CIGAR code %s "
                             " which is not currently implemented" %
                             BAM_CIGAR_DICT[op])

    # are any of the SNPs indels? If so, discard.
    for i in snp_idx:
        if is_indel(snp_tab[i]):
            overlap_indel = True
            return (None, None, is_split, overlap_indel)

    n_overlap_snps = len(read_offsets)
    if n_overlap_snps == 0:
        # no SNPs overlap this read
        return (None, None, is_split, overlap_indel)

    if hap_tab:
        # genotype info is provided by haplotype table
        # pull out subset of overlapping SNPs that are heterozygous
        # in this individual
        het_read_offsets = []
        het_snp_idx = []
        for (i, read_offset) in zip(snp_idx, read_offsets):
            haps = hap_tab[i, (ind_idx*2):(ind_idx*2 + 2)]

            if ind_idx*2 > hap_tab.shape[1]:
                raise ValueError("index of individual (%d) is >= number of "
                                 "individuals in haplotype_tab (%d). probably "
                                 "need to specify --population or use a different "
                                 "--samples_tab" % (ind_idx, hap_tab.shape[1]/2))

            if haps[0] != haps[1]:
                # this is a het
                het_read_offsets.append(read_offset)
                het_snp_idx.append(i)

        n_overlap_hets = len(het_read_offsets)

        if n_overlap_hets == 0:
            # none of the overlapping SNPs are hets
            return (None, None, is_split, overlap_indel)

        if n_overlap_hets == 1:
            # only one overlapping SNP is a het
            return (het_snp_idx[0], het_read_offsets[0], is_split, overlap_indel)

        # choose ONE overlapping HETEROZYGOUS SNP randomly to add counts to
        # we don't want to count same read multiple times
        r = np.random.randint(0, n_overlap_hets)
        return (het_snp_idx[r], het_read_offsets[r], is_split, overlap_indel)

    else:
        # We don't have haplotype tab, so we don't know which SNPs are
        # heterozygous in this individual. But we can still tell
        # whether read sequence matches reference or non-reference
        # allele. Choose ONE overlapping SNP randomly to add counts to
        if n_overlap_snps == 1:
            return (snp_idx[0], read_offsets[0], is_split, overlap_indel)
        else:
            r = np.random.randint(0, n_overlap_snps)
            return (snp_idx[r], read_offsets[r], is_split, overlap_indel)




def add_read_count(read, chrom, ref_array, alt_array, other_array,
                   read_count_array, snp_index_array, snp_tab, hap_tab,
                   warned_pos, max_count, ind_idx):

    # pysam positions start at 0
    start = read.pos+1
    end = read.aend

    if start < 1 or end > chrom.length:
        sys.stderr.write("WARNING: skipping read aligned past end of "
                         "chromosome. read: %d-%d, %s:1-%d\n" %
                         (start, end, chrom.name, chrom.length))
        return


    if read.qlen != read.rlen:
        sys.stderr.write("WARNING skipping read: handling of "
                         "partially mapped reads not implemented\n")
        return

    # look for SNPs that overlap mapped read position, and if there
    # are more than one, choose one at random
    snp_idx, read_offset, is_split, overlap_indel = \
      choose_overlap_snp(read, snp_tab, snp_index_array, hap_tab, ind_idx)

    if overlap_indel:
        return

    # store counts of reads at start position
    if read_count_array[start-1] < max_count:
        read_count_array[start-1] += 1
    else:
        if not start in warned_pos:
            sys.stderr.write("WARNING read count at position %d "
                             "exceeds max %d\n" % (start, max_count))
            warned_pos[start] = True


    if snp_idx is None:
        return

    snp = snp_tab[snp_idx]

    base = read.seq[read_offset]
    snp_pos = snp['pos']

    if base == snp['allele1']:
        # matches reference allele
        if ref_array[snp_pos-1] < max_count:
            ref_array[snp_pos-1] += 1
        elif not snp_pos in warned_pos:
            sys.stderr.write("WARNING ref allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, max_count))
            warned_pos[snp_pos] = True
    elif base == snp['allele2']:
        # matches alternate allele
        if alt_array[snp_pos-1] < max_count:
            alt_array[snp_pos-1] += 1
        elif not snp_pos in warned_pos:
            sys.stderr.write("WARNING alt allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, max_count))
            warned_pos[snp_pos] = True
    else:
        # matches neither
        if other_array[snp_pos-1] < max_count:
            other_array[snp_pos-1] += 1
        elif not snp_pos in warned_pos:
            sys.stderr.write("WARNING other allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, max_count))






def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--chrom",
                        help="Path to chromInfo.txt file (may be gzipped) "
                        "with list of chromosomes for the relevant genome "
                        "assembly. Each line in file should contain "
                        "tab-separated chromosome name and chromosome length "
                        "(in basepairs). chromInfo.txt files can be "
                        "downloaded from the UCSC genome browser. For "
                        "example, a chromInfo.txt.gz file for hg19 can "
                        "be downloaded from "
                        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/",
                        metavar="CHROM_TXT_FILE",
                        required=True)


    parser.add_argument("--snp_index",
                        help="Path to HDF5 file containing SNP index. The "
                        "SNP index is used to convert the genomic position "
                        "of a SNP to its corresponding row in the haplotype "
                        "and snp_tab HDF5 files.",
                        metavar="SNP_INDEX_H5_FILE",
                        required=True)

    parser.add_argument("--snp_tab",
                        help="Path to HDF5 file to read SNP information "
                        "from. Each row of SNP table contains SNP name "
                        "(rs_id), position, allele1, allele2.",
                        metavar="SNP_TABLE_H5_FILE",
                        required=True)

    parser.add_argument("--haplotype",
                        help=" Path to HDF5 file to read phased haplotypes "
                        "from. If supplied, when read overlaps multiple SNPs "
                        "counts are randomly assigned to ONE of the "
                        "overlapping HETEROZYGOUS SNPs; if not supplied "
                        "counts are randomly assigned to ONE of overlapping "
                        "SNPs (regardless of their genotype).",
                        metavar="HAPLOTYPE_H5_FILE",
                        default=None)

    parser.add_argument("--samples",
                        help="Path to text file containing a list of "
                        "individual identifiers. The ordering of individuals "
                        "must be consistent with the haplotype file. The "
                        "samples file is assumed to have one identifier per "
                        "line in the first column (other columns are "
                        "ignored).",
                        metavar="SAMPLES_TXT_FILE",
                        default=None)

    parser.add_argument("--individual",
                        help="Identifier for individual, used to determine "
                        "which SNPs are heterozygous. Must be provided "
                        "if --haplotype argument is provided and must "
                        "match one of the individuals in the file provided "
                        "with --samples argument.",
                        metavar="INDIVIDUAL",
                        default=None)

    parser.add_argument("--data_type",
                        help="Data type of counts stored in HDF5 files. "
                        "uint8 requires less disk space but has a "
                        "maximum value of 255."
                        "(default=uint8)", choices=("uint8", "uint16"),
                        default="uint8")

    parser.add_argument("--ref_as_counts",
                        help="Path to HDF5 file to write counts of reads "
                        "that match reference allele. Allele-specific counts "
                        "are stored at the position of the SNP."
                        "that match reference",
                        metavar="REF_AS_COUNT_H5_FILE",
                        required=True)

    parser.add_argument("--alt_as_counts",
                        help="Path to HDF5 file to write counts of reads "
                        "that match alternate allele. Allele-specific counts "
                        "are stored at the position of the SNP.",
                        metavar="ALT_AS_COUNT_H5_FILE",
                        required=True)

    parser.add_argument("--other_as_counts",
                        help="Path to HDF5 file to write counts of reads "
                        "that match neither reference nor alternate allele. "
                        "Allele-specific counts are stored at the position "
                        "of the SNP.",
                        metavar="OTHER_COUNT_H5_FILE",
                        required=True)

    parser.add_argument("--read_counts",
                       help="Path to HDF5 file to write counts of all "
                       "reads, regardless of whether they overlap a SNP. "
                       "Read counts are stored at the left-most position "
                       "of the mapped read.",
                       metavar="READ_COUNT_H5_FILE",
                       required=True)

    parser.add_argument("bam_filenames", action="store", nargs="+",
                        help="BAM file(s) to read mapped reads from. "
                        "BAMs must be sorted and indexed.")

    args = parser.parse_args()

    if args.haplotype and (args.individual is None or args.samples is None):
            parser.error("--indidivual and --samples arguments "
                         "must also be provided when --haplotype argument "
                         "is provided")


    return args




def lookup_individual_index(samples_file, ind_name, population=None):
    """Gets the index of individual that is used
    to lookup information in the genotype and haplotype tables"""
    f = open(samples_file)

    if population:
        p = population.lower()
    else:
        p = None

    idx = 0
    for line in f:
        if line.startswith("samples"):
            # header line
            continue

        words = line.rstrip().split()
        name = words[0].replace("NA", "")

        if len(words) > 1:
            pop = words[1].lower()
        else:
            pop = ""

        if len(words) > 2:
            group = words[2].lower()
        else:
            group = ""

        # if specified, only consider a single population or group
        if p and pop != p and group != p:
            continue

        if name == ind_name:
            f.close()
            return idx

        idx += 1


    raise ValueError("individual %s (with population=%s) "
                     "is not in samples file %s" %
                     (ind_name, population, samples_file))




def main():
    args = parse_args()

    snp_tab_h5 = tables.openFile(args.snp_tab, "r")
    snp_index_h5 = tables.openFile(args.snp_index, "r")

    if args.haplotype:
        hap_h5 = tables.openFile(args.haplotype, "r")
        ind_idx = lookup_individual_index(args.samples, args.individual)
    else:
        hap_h5 = None
        ind_idx = None

    ref_count_h5 = tables.openFile(args.ref_as_counts, "w")
    alt_count_h5 = tables.openFile(args.alt_as_counts, "w")
    other_count_h5 = tables.openFile(args.other_as_counts, "w")
    read_count_h5 = tables.openFile(args.read_counts, "w")

    output_h5 = [ref_count_h5, alt_count_h5, other_count_h5, read_count_h5]

    chrom_dict = {}

    # initialize every chromosome in output files
    chrom_list = chromosome.get_all_chromosomes(args.chrom)

    for chrom in chrom_list:
        for out_file in output_h5:
            create_carray(out_file, chrom, args.data_type)
        chrom_dict[chrom.name] = chrom

    count = 0
    dtype = None
    if args.data_type == "uint8":
        max_count = MAX_UINT8_COUNT
        dtype = np.uint8
    elif args.data_type == "uint16":
        max_count = MAX_UINT16_COUNT
        dtype = np.uint16
    else:
        raise NotImplementedError("unsupported datatype %s" % args.data_type)

    for chrom in chrom_list:
        sys.stderr.write("%s\n" % chrom.name)

        warned_pos = {}

        # fetch SNP info for this chromosome
        if chrom.name not in snp_tab_h5.root:
            # no SNPs for this chromosome
            continue

        sys.stderr.write("fetching SNPs\n")

        snp_tab = snp_tab_h5.getNode("/%s" % chrom.name)
        snp_index_array = snp_index_h5.getNode("/%s" % chrom.name)[:]
        if hap_h5:
            hap_tab = hap_h5.getNode("/%s" % chrom.name)
        else:
            hap_tab = None

        # initialize count arrays for this chromosome to 0
        ref_carray = get_carray(ref_count_h5, chrom)
        alt_carray = get_carray(alt_count_h5, chrom)
        other_carray = get_carray(other_count_h5, chrom)
        read_count_carray = get_carray(read_count_h5, chrom)

        ref_array = np.zeros(chrom.length, dtype)
        alt_array = np.zeros(chrom.length, dtype)
        other_array = np.zeros(chrom.length, dtype)
        read_count_array = np.zeros(chrom.length, dtype)

        # loop over all BAM files, pulling out reads
        # for this chromosome
        for bam_filename in args.bam_filenames:
            sys.stderr.write("reading from file %s\n" % bam_filename)

            samfile = pysam.Samfile(bam_filename, "rb")

            for read in get_sam_iter(samfile, chrom):
                count += 1
                if count == 10000:
                    sys.stderr.write(".")
                    count = 0

                add_read_count(read, chrom, ref_array, alt_array,
                               other_array, read_count_array,
                               snp_index_array, snp_tab, hap_tab,
                               warned_pos, max_count, ind_idx)

            # store results for this chromosome
            ref_carray[:] = ref_array
            alt_carray[:] = alt_array
            other_carray[:] = other_array
            read_count_carray[:] = read_count_array
            sys.stderr.write("\n")

            samfile.close()

    # set track statistics and close HDF5 files

    sys.stderr.write("setting statistics for each chromosome\n")
    for h5f in output_h5:
        chromstat.set_stats(h5f, chrom_list)
        h5f.close()

    snp_tab_h5.close()
    snp_index_h5.close()
    if hap_h5:
        hap_h5.close()


    sys.stderr.write("done\n")


main()



