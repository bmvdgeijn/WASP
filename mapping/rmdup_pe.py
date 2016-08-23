import random
import pysam
import os
import sys
import argparse



class ReadStats(object):

    def __init__(self):
        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0

        # paired reads map to different chromosomes
        self.discard_different_chromosome = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # reads where we expected to see other pair, but it was missing
        # possibly due to read-pairs with different names
        self.discard_missing_pair = 0

        # reads with only one paired mapped
        self.discard_single = 0

        # reads discarded because duplicated
        self.discard_dup = 0
        
        # number of read pairs kept
        self.keep_pair = 0
        

    def write(self, file_handle):
        sys.stderr.write("DISCARD reads:\n"
                         "  improper pair: %d\n"
                         "  different chromosome: %d\n"
                         "  secondary alignment: %d\n"
                         "  missing pairs (e.g. mismatched read names): %d\n"
                         "  not paired: %d\n"
                         "  duplicate pairs: %d\n"
                         "KEEP reads:\n"
                         "  pairs: %d\n"  %
                         (self.discard_improper_pair,
                          self.discard_different_chromosome,
                          self.discard_secondary,
                          self.discard_missing_pair,
                          self.discard_single,
                          self.discard_dup,
                          self.keep_pair))
        

                


def main(input_bam, output_bam):
    if input_bam.endswith(".sam") or input_bam.endswith("sam.gz"):
        infile = pysam.Samfile(input_bam, "r")
    else:
        # assume binary BAM file
        infile = pysam.Samfile(input_bam, "rb")

    if output_bam.endswith(".sam"):
        # output in text SAM format
        outfile = pysam.Samfile(output_bam, "w", template=infile)
    elif output_bam.endswith(".bam"):
        # output in binary compressed BAM format
        outfile = pysam.Samfile(output_bam, "wb", template=infile)
    else:
        raise ValueError("name of output file must end with .bam or .sam")

    filter_reads(infile, outfile)

    infile.close()
    outfile.close()



def filter_reads(infile, outfile):
    read_stats = ReadStats()
    
    cur_tid = None
    seen_chrom = set([])

    # name of reads to keep
    keep_cache = {}
    # name of reads to discard
    discard_cache = {}

    keep_cache_size = 0
    discard_cache_size = 0
    read_count = 0

    # current position on chromosome
    cur_pos = None
    # lists of reads at current position,
    # grouped by the mate pair position
    cur_by_mpos = {}
    
    for read in infile:
        read_count += 1
                            
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome
            cur_chrom = infile.getrname(read.tid)

            if keep_cache_size + discard_cache_size != 0:
                sys.stderr.write("WARNING: failed to find pairs for %d "
                                 "reads on this chromosome\n" %
                                 (keep_cache_size + discard_cache_size))
                read_stats.discard_missing_pair += keep_cache_size + \
                                                   discard_cache_size

            keep_cache = {}
            discard_cache = {}
            keep_cache_size = 0
            discard_cache_size = 0
            cur_pos = None
            cur_by_mpos = {}
            read_count = 0
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)
            sys.stderr.write("processing reads\n")

        if read.is_secondary:
            # this is a secondary alignment (i.e. read was aligned more than
            # once and this has align score that <= best score)
            read_stats.discard_secondary += 1
            continue

        if (not read.is_paired) or (read.next_reference_name is None):
            read_stats.discard_single +=1
            continue

        if (read.next_reference_name != cur_chrom) and \
           (read.next_reference_name != "="):
            # other side of pair mapped to different chromosome
            read.stats.discard_different_chromosome += 1
            continue

        if not read.is_proper_pair:
            read_stats.discard_improper_pair += 1
            continue

        if (cur_pos is not None) and (read.pos < cur_pos):
            raise ValueError("expected input BAM file to be sorted "
                             "but reads are out of order")
        
        if cur_pos is None or read.pos > cur_pos:
            # we have advanced to a new start position
            # decide which of reads at last position to keep or discard
            for mpos, read_list in cur_by_mpos.items():
                # only keep one read from list with same pos,mate_pos pair
                # shuffle order of reads in list and take first
                # as 'keep' read
                random.shuffle(read_list)
                keep_read = read_list.pop()
                if keep_read.qname in keep_cache:
                    raise ValueError("read %s is already "
                                     "in keep cache" % keep_read.qname)
                keep_cache[keep_read.qname] = keep_read
                keep_cache_size += 1
                
                # rest of reads get discarded
                for discard_read in read_list:
                    if discard_read in discard_cache:
                        raise ValueError("read %s is already in discard cache" %
                                         discard_read.qname)
                    discard_cache[discard_read.qname] = discard_read
                    discard_cache_size += 1

            # create new list of reads at current position
            cur_pos = read.pos
            cur_by_mpos = {}

        if read.qname in keep_cache:
            # we already saw prev side of pair, retrieve from cache
            read1 = keep_cache[read.qname]
            read2 = read
            del keep_cache[read.qname]
            keep_cache_size -= 1

            if read2.next_reference_start != read1.reference_start:
                sys.stderr.write("WARNING: read pair positions "
                                 "do not match for pair %s\n" % read.qname)

            read_stats.keep_pair += 1
            outfile.write(read1)
            outfile.write(read2)
            
        elif read.qname in discard_cache:
            # we already saw prev side of pair, but decided to discard
            # because read duplicated
            del discard_cache[read.qname]
            discard_cache_size -= 1
            read_stats.discard_dup += 1

        else:
            # we have not seen other side of this read yet
            # add read to list of those at current position
            # grouping by mate-pair position
            if read.mpos in cur_by_mpos:
                cur_by_mpos[read.mpos].append(read)
            else:
                cur_by_mpos[read.mpos] = [read]

    if (keep_cache_size + discard_cache_size) != 0:
        sys.stderr.write("WARNING: failed to find pairs for %d "
                         "keep reads and %d discard reads on this "
                         "chromosome\n" % (keep_cache_size, discard_cache_size))
        
        read_stats.discard_missing_pair += keep_cache_size + discard_cache_size

    read_stats.write(sys.stderr)
    
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help="input BAM or SAM file (must "
                        "be sorted!)")
    parser.add_argument("output_bam", help="output BAM or SAM file (not "
                        "sorted!)")
    
    options = parser.parse_args()
    
    main(options.input_bam, options.output_bam)
