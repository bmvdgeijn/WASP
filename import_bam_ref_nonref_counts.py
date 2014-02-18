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
usage: import_bam_ref_nonref_counts.py [-h] [--assembly ASSEMBLY] \
                                        ref_as_count_track \
                                        alt_as_count_track \
                                        other_as_count_track \
                                        read_count_track \
                                        bam_file1 [bam_file2 [...]]

positional arguments:
  ref_as_count_track    name of track to store counts of allele-specific reads that match
                        reference allele
  alt_as_count_track    name of track to store counts of allele-specific reads that match
                        alternate (non-reference) allele
  other                 name of track to store counts of allele-specific reads that match
                        neither allele
  read_count_track      name of track to store read counts in (stored at position of left end
                        of read)
  bam_file(s)           BAM files containing filtered and sorted mapped reads.
  input_file            bed-like file to read coordinates of test SNP and
                        target region from

optional arguments:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   genome assembly that reads were mapped to (e.g. hg18)


This program reads BAM files and counts the number of reads that
match the alternate and reference allele at every SNP position
stored in the /impute2/snps HDF5 track. The read counts are stored in
specified ref_track, alt_track and other_track HDF5 tracks. Additionally counts
of all reads are stored in another track (at the left-most position 
of the reads).

Currently this program filters out reads that:
  1) do not map uniquely according to Roger's 20bp mapping criteria
  2) reads that overlap known indels
  3) reads that overlap more than 1 known SNP

Filtering Criterium 1 uses a hardcoded track containing Roger's
mappability data. 

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

SNP_UNDEF = -1
SNP_TRACK_NAME = "impute2/snps"
SNP_INDEX_TRACK_NAME = "impute2/snp_index"
SNP_REF_MATCH_TRACK_NAME = "impute2/snp_match_ref"
DEL_TRACK_NAME = "impute2/deletions"

SEQ_TRACK_NAME = "seq"

MAP_TRACK_NAME = "mappability/roger_20bp_mapping_uniqueness"

MAX_VAL = 255


def create_carray(track, chrom):
    atom = tables.UInt8Atom(dflt=0)
    
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


    

def add_read_count(read, chrom, ref_array, alt_array, other_array,
                   read_count_array, snp_index_array, snp_tab, snp_ref_match,
                   del_array, map_array, filter_counts):

    # pysam positions start at 0
    start = read.pos+1
    end = read.aend

    if start < 1 or end > chrom.length:
        sys.stderr.write("WARNING: skipping read aligned past end of "
                         "chromosome. read: %d-%d, %s:1-%d\n" %
                         (start, end, chrom.name, chrom.length))
        return


    if read.qlen != read.rlen:
        sys.stderr.write("WARNING skipping read, because handling of "
                         "partially mapped reads not implemented "
                         "(maybe due to clipping?)\n")
        return
        
    map_code = map_array[start-1]    
    if map_code != 1:
        # skip, not uniquely mappable
        filter_counts['MAP'] += 1
        return

    if np.any(del_array[start-1:end]):
        # read overlaps bases that are deleted in non-reference
        filter_counts['DELETION'] += 1
        return
    
    idx = snp_index_array[start-1:end]
        
    # get offsets within read of any overlapping SNPs
    read_offsets = np.where(idx != SNP_UNDEF)[0]
    n_overlapping_snps = read_offsets.size

    if n_overlapping_snps > 1:
        # read overlaps multiple SNPs--discard
        #sys.stderr.write("discarding read that overlaps %d SNPs\n" %
        #                 n_overlapping_snps)
        filter_counts['MULTI_SNP'] += 1
        return
    
        
    if n_overlapping_snps == 1:
        # store the allele-specific counts at this SNP        
        offset = read_offsets[0]

        snp = snp_tab[idx[offset]]

        match_ref = snp_ref_match[idx[offset]]
        if not match_ref:
            # SNP was flagged because reference allele
            # does not match reference sequence
            filter_counts['SNP_MISMATCH_REF'] += 1
            return

        if is_indel(snp):
            # since we already checked for deletion, must be insertion
            filter_counts['INSERTION'] += 1
            return

        snp_pos = snp['pos']

        # check against reference sequence
        # there are a few SNPs that do not match the reference
        # and we'd like to throw those out
        # ref_base = chr(seq_vals[snp_pos-1])
        # if snp['allele1'] != ref_base:
        #     sys.stderr.write("discarding read that overlaps SNP that "
        #                      "does not match reference\n")
        #     return

        base = read.seq[offset]

        if base == snp['allele1']:
            # matches reference allele
            if ref_array[snp_pos-1] < MAX_VAL:
                ref_array[snp_pos-1] += 1
            else:
                sys.stderr.write("WARNING ref allele count at position %d "
                                 "exceeds max %d\n" % (snp_pos, MAX_VAL))
        elif base == snp['allele2']:
            # matches alternate allele
            if alt_array[snp_pos-1] < MAX_VAL:
                alt_array[snp_pos-1] += 1
            else:
                sys.stderr.write("WARNING alt allele count at position %d "
                                 "exceeds max %d\n" % (snp_pos, MAX_VAL))
        else:
            # matches neither
            if other_array[snp_pos-1] < MAX_VAL:
                other_array[snp_pos-1] += 1
            else:
                sys.stderr.write("WARNING other allele count at position %d "
                                 "exceeds max %d\n" % (snp_pos, MAX_VAL))
                

    # Store counts of reads that passed filtering. For this 
    # we are counting reads that overlap 0 or 1 SNPs
    if read_count_array[start-1] < MAX_VAL:
        read_count_array[start-1] += 1
    else:
        sys.stderr.write("WARNING read count at position %d "
                         "exceeds max %d\n" % (start-1, MAX_VAL))
        
    



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default=None)
    
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
                        help="bam file(s) to read data from")

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

    snp_track = gdb.open_track(SNP_TRACK_NAME)
    snp_index_track = gdb.open_track(SNP_INDEX_TRACK_NAME)
    ref_match_track = gdb.open_track(SNP_REF_MATCH_TRACK_NAME)
    del_track = gdb.open_track(DEL_TRACK_NAME)

    map_track = gdb.open_track(MAP_TRACK_NAME)

    chrom_dict = {}

    count = 0
    
    filter_counts = {'MAP' : 0,
                     'MULTI_SNP' : 0,
                     'SNP_MISMATCH_REF' : 0,
                     'DELETION' : 0,
                     'INSERTION' : 0}
    
    for chrom in gdb.get_chromosomes(get_x=False):
        sys.stderr.write("%s\n" % chrom.name)
        # initialize output tracks
        ref_carray = create_carray(ref_count_track, chrom)
        alt_carray = create_carray(alt_count_track, chrom)
        other_carray = create_carray(other_count_track, chrom)
        read_count_carray = create_carray(read_count_track, chrom)
        
        ref_array = np.zeros(chrom.length, np.uint8)
        alt_array = np.zeros(chrom.length, np.uint8)
        other_array = np.zeros(chrom.length, np.uint8)
        read_count_array = np.zeros(chrom.length, np.uint8)

        # fetch SNPs and indel info for this chromosome
        sys.stderr.write("fetching SNPs\n")
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)
        snp_index_array = snp_index_track.get_nparray(chrom)
        snp_ref_match = ref_match_track.h5f.getNode("/%s" % chrom.name)
        del_array = del_track.get_nparray(chrom)

        map_array = map_track.get_nparray(chrom)

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
                               snp_index_array, snp_tab,
                               snp_ref_match, del_array, 
                               map_array, filter_counts)


            # store results for this chromosome        
            ref_carray[:] = ref_array
            alt_carray[:] = alt_array
            other_carray[:] = other_array
            read_count_carray[:] = read_count_array
            sys.stderr.write("\n")

            samfile.close()

    sys.stderr.write("FILTERED reads:\n")
    sys.stderr.write("  non-unique mappability: %d\n" % 
                     filter_counts['MAP'])
    sys.stderr.write("  overlap multiple SNPs: %d\n" % 
                     filter_counts['MULTI_SNP'])
    sys.stderr.write("  overlap deletion: %d\n" % 
                     filter_counts['DELETION'])
    sys.stderr.write("  overlap insertion: %d\n" % 
                     filter_counts['INSERTION'])
    sys.stderr.write("  overlap SNP with incorrect ref allele: %d\n" %
                     filter_counts['SNP_MISMATCH_REF'])
        
    for track in output_tracks:
        track.close()
        
    snp_track.close()
    snp_index_track.close()
    ref_match_track.close()
    del_track.close()
    map_track.close()

    
main()
        
        
    
