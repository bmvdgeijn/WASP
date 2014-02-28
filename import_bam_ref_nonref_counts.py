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
usage: import_bam_ref_nonref_counts.py [-h] [--assembly ASSEMBLY]
                                       [--snp_index_track SNP_INDEX_TRACK]
                                       [--snp_track SNP_TRACK]
                                       [--hap_track HAP_TRACK]
                                       ref_as_count_track alt_as_count_track
                                       other_as_count_track read_count_track
                                       bam_filenames [bam_filenames ...]

positional arguments:
  ref_as_count_track    name of track to store counts of reads that match
                        reference
  alt_as_count_track    name of track to store counts of reads that match
                        alternate
  other_as_count_track  name of track to store counts of reads that match
                        neither reference or alternate
  read_count_track      name of track to store read counts in--positions of
                        LEFT end of read are used
  bam_filenames         bam file(s) to read mapped reads from

optional arguments:
  -h, --help            show help message and exit
  --assembly ASSEMBLY   genome assembly that reads were mapped to (e.g. hg19)
  --snp_index_track SNP_INDEX_TRACK
                        name of SNP index track
                        (default=1000genomes/snp_index)
  --snp_track SNP_TRACK
                        name of track containing table of SNPs
                        (default=1000genomes/)
  --data_type {uint8,uint16}
                        data type of stored counts; uint8 takes up less disk
                        space but has a maximum value of 255 (default=uint8)
  --hap_track HAP_TRACK
                        name of haplotype track; if supplied when read
                        overlaps multiple SNPs counts are randomly assigned to
                        ONE of the overlapping HETEROZYGOUS SNPs; if not
                        supplied counts are randomly assigned to ONE of
                        overlapping SNPs (regardless of genotype)

This program reads BAM files and counts the number of reads that match
the alternate and reference allele at every SNP position stored in the
/impute2/snps HDF5 track. The read counts are stored in specified
ref_track, alt_track and other_track HDF5 tracks. Additionally counts
of all reads are stored in another track (at the left-most position of
the reads).

This program does not perform filtering of reads based on mappability.
It is assumed that this filtering will be done prior to calling this
script.

Reads that overlap known indels are not included in allele-specific
counts.

We are currently working on an improved allele-specific mapping method
so that criteria 1 and 2 can be relaxed can be done at the level of
the BAM (before running this script).

This program requires the use of the genome library, which can be
downloaded from here:
https://github.com/gmcvicker/genome
"""

import sys
import os
import gzip

import tables
import argparse
import numpy as np

import pysam

import genome.db
import genome.coord
import genome.trackstat as trackstat


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


def create_carray(track, chrom, data_type):
    if data_type == "uint8":
        atom = tables.UInt8Atom(dflt=0)
    elif data_type == "uint16":
        atom = tables.UInt16Atom(dflt=0)
    else:
        raise NotImplementedError("unsupported datatype %s" % data_type)
        
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray



def get_carray(track, chrom):
    return track.h5f.getNode("/%s" % chrom)




def is_indel(snp):
    if (len(snp['allele1']) != 1) or (len(snp['allele2'])) != 1:
        return True


def dump_read(f, read):
    cigar_str = " ".join(["%s:%d" % (BAM_CIGAR_DICT[c[0]], c[1]) for c in read.cigar])

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
            sys.stderr.write("WARNING: %s does not exist in BAM file, "
                             "skipping chromosome\n" % chrom_name)
            sam_iter = iter([])

    return sam_iter




def choose_overlap_snp(read, snp_tab, snp_index_array, hap_tab):
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
        else:
            sys.stderr.write("skipping because contains CIGAR code %s "
                             " which is not currently implemented" % BAM_CIGAR_DICT[op])
            
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
        # pull out subset of overlapping SNPs that are heterozygous in this individual
        het_read_offsets = []
        for (i, read_offset) in snp_idx:
            haps = hap_tab[i, (ind_idx*2):(ind_idx*2 + 2)]
            if haps[0] != haps[1]:
                # this is a het
                het_read_offsets.append(read_offset)
                het_snp_idx.append(i)

        n_overlap_hets = len(het_read_offsets)

        if n_overlap_hets == 0:
            # none of the overlapping SNPs are hets
            return (None, None, is_split, overlap_indel)

        # choose ONE overlapping het SNP randomly to add counts to
        # we don't want to count same read multiple times
        r = np.random.randint(0, n_overlap_hets-1)
        return (het_snp_idx[r], het_read_offsets[r], is_split, overlap_indel)
    
    else:
        # We don't have haplotype tab, so we don't know which SNPs are
        # heterozygous in this individual. But we can still tell
        # whether read sequence matches reference or non-reference
        # allele. Choose ONE overlapping SNP randomly to add counts to
        if n_overlap_snps == 1:
            return (snp_idx[0], read_offsets[0], is_split, overlap_indel)
        else:
            r = np.random.randint(0, n_overlap_snps-1)
            return (snp_idx[r], read_offsets[r], is_split, overlap_indel)



    
def add_read_count(read, chrom, ref_array, alt_array, other_array,
                   read_count_array, snp_index_array, snp_tab, hap_tab, 
                   warned_pos, max_count):
    
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
      choose_overlap_snp(read, snp_tab, snp_index_array, hap_tab)

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

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg19)", default=None)

    parser.add_argument("--snp_index_track",
                        help="name of SNP index track "
                        "(default=1000genomes/snp_index)",
                        default="1000genomes/snp_index")

    parser.add_argument("--snp_track",
                        help="name of track containing table of SNPs"
                        " (default=1000genomes/snp_tab)",
                        default="1000genomes/snp_tab")

    parser.add_argument("--data_type",
                        help="data type of stored counts; uint8 takes "
                        "up less disk space but has a maximum value of 255 "
                        "(default=uint8)", choices=("uint8", "uint16"), 
                        default="uint8")
                        
    parser.add_argument("--hap_track",
                        help="name of haplotype track; if supplied when "
                        "read overlaps multiple SNPs counts are randomly "
                        "assigned to ONE of the overlapping HETEROZYGOUS "
                        "SNPs; if not supplied counts are randomly assigned "
                        "to ONE of overlapping SNPs (regardless of genotype)",
                        default=None)
        
    parser.add_argument("ref_as_count_track",
                        help="name of track to store counts of reads "
                        "that match reference")

    parser.add_argument("alt_as_count_track", 
                        help="name of track to store counts of reads "
                        "that match alternate")

    parser.add_argument("other_as_count_track", 
                        help="name of track to store counts of reads "
                        "that match neither reference or alternate")

    parser.add_argument("read_count_track",
                       help="name of track to store read counts in--"
                       "positions of LEFT end of read are used")
    
    parser.add_argument("bam_filenames", action="store", nargs="+",
                        help="BAM file(s) to read mapped reads from. "
                        "BAMs must be sorted and indexed.")
    

    args = parser.parse_args()

    return args
    


        
def main():
    args = parse_args()
    
    # create a database track
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    ref_count_track = gdb.create_track(args.ref_as_count_track)
    alt_count_track = gdb.create_track(args.alt_as_count_track)
    other_count_track = gdb.create_track(args.other_as_count_track)
    read_count_track = gdb.create_track(args.read_count_track)
    
    output_tracks = [ref_count_track, alt_count_track, 
                     other_count_track, read_count_track]

    snp_track = gdb.open_track(args.snp_track)
    snp_index_track = gdb.open_track(args.snp_index_track)
    if args.hap_track:
        hap_track = gdb.open_track(args.hap_track)
    else:
        hap_track = None

    chrom_dict = {}

    count = 0

    # initialize every chromosome in output tracks
    for chrom in gdb.get_all_chromosomes():
        for track in output_tracks:
            create_carray(track, chrom, args.data_type)
        chrom_dict[chrom.name] = chrom

    count = 0

    if args.data_type == "uint8":
        max_count = MAX_UINT8_COUNT
    elif args.data_type == "uint16":
        max_count = MAX_UINT16_COUNT
    else:
        raise NotImplementedError("unsupported datatype %s" % args.data_type)
    
    for chrom in gdb.get_chromosomes(get_x=False):
        sys.stderr.write("%s\n" % chrom.name)

        if chrom.name != "chr22":
            continue

        warned_pos = {}

        # initialize count arrays for this chromosome to 0
        ref_carray = get_carray(ref_count_track, chrom)
        alt_carray = get_carray(alt_count_track, chrom)
        other_carray = get_carray(other_count_track, chrom)
        read_count_carray = get_carray(read_count_track, chrom)
        
        ref_array = np.zeros(chrom.length, np.uint8)
        alt_array = np.zeros(chrom.length, np.uint8)
        other_array = np.zeros(chrom.length, np.uint8)
        read_count_array = np.zeros(chrom.length, np.uint8)

        # fetch SNP info for this chromosome
        sys.stderr.write("fetching SNPs\n")
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)
        snp_index_array = snp_index_track.get_nparray(chrom)
        if hap_track:
            hap_tab = hap_track.h5f.getNode("/%s" % chrom.name)
        else:
            hap_tab = None
        
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
                               warned_pos, max_count)

            # store results for this chromosome        
            ref_carray[:] = ref_array
            alt_carray[:] = alt_array
            other_carray[:] = other_array
            read_count_carray[:] = read_count_array
            sys.stderr.write("\n")

            samfile.close()

    # set track statistics and close HDF5 files
    sys.stderr.write("setting track statistics\n")
    for track in output_tracks:
        sys.stderr.write("%s\n" % track.name)
        trackstat.set_stats(gdb, track)
        track.close()        
    snp_track.close()
    snp_index_track.close()    
    if hap_track:
        hap_track.close()

    
main()
        
        
    
