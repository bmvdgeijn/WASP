def run(to_remap_bam, remap_bam, keep_bam, orig_num_file, is_paired_end):
    import gzip
    import sys

    import pysam

    to_remap_bam = pysam.Samfile(to_remap_bam, "rb")
    remap_bam = pysam.Samfile(remap_bam, "rb")
    keep_bam = pysam.Samfile(keep_bam, "wb", template=to_remap_bam)
    orig_num_file = gzip.open(orig_num_file)

    # correct_maps is a list of reads that mapped correctly. The read is
    # is represented by its "ID number" (which is the order that it was seen in
    # the original input bam file and also the first number in its new read
    # name, e.g. the 42 in @42:chr1:14536:1). For paired end, both reads will
    # have the same ID number. This means that the first read pair will have ID
    # number 1, the second will be 2, and so on. For a read or read pair to be
    # remapped correctly, the ID number X of that read/read pair should show up
    # N amount of times in correct_maps where N is the value in orig_num_file at
    # line X (one-based).
    correct_maps = []
    end_of_file = False
    
    # Get a list of reads that remapped correctly.
    remap_read = remap_bam.next()
    
    while not end_of_file:
        """
        The read names in the fastq files for remapping consist of four or
        five parts (depending on single (1:chr1:763045:1) or paired end data
        (1:chr1:763006:763045:3)) separated by colons. The first part is the
        remap number. The first read or read pair from the input file to be
        remapped is 1, the second read pair is 2, etc. The second part is the
        chromosome name. The third part is the alignment position of the
        original read pair. If paired end data, the third part is the
        alignment position of the left read (relative to the reference) and
        the fourth part is the position of the right read. This is
        irrespective of which read is R1 and R2. The left reads are written
        into the R1 fastq file and the right reads are written into the R2
        fastq file. The last part is the number of alternative sequences for
        the read pair (i.e. the number of new sequences with swapped alleles
        we have to map to see if they match with the alignment of the original
        sequence). Note that if a read pair overlaps multiple SNPs and has
        multiple alternate sequences, each of those alternate sequences will
        have the exact same read name.
        """

        chrm = remap_read.qname.strip().split(":")[1]
        
        if remap_read.is_reverse:
            if is_paired_end:
                pos = int(remap_read.qname.strip().split(":")[3])
            else:
                pos = int(remap_read.qname.strip().split(":")[2])
        else:
            pos = int(remap_read.qname.strip().split(":")[2])
        read_num = int(remap_read.qname.strip().split(":")[0])
        if (remap_read.tid != -1 and remap_read.pos == pos and 
            remap_bam.getrname(remap_read.tid) == chrm):
            dels = 0
            # Throw out the remapped read if it remapped with a deletion...for
            # now.
            for cig in remap_read.cigar:
                if not cig[0] in (0, 3, 4):
                    dels += 1
            if dels == 0:
                correct_maps.append(read_num)
        try:
            remap_read = remap_bam.next()
        except:
            end_of_file = True
    
    correct_maps.sort()

    # Pull out original aligned reads if all of the alternatives mapped
    # correctly.
    orig_read = to_remap_bam.next()
    # orig_num is the number of different reads generated from the original read
    # (pair) (depends on number of alleles it overlapped). For paired end data,
    # this is number is at least two (one for each read) and will always be a
    # multiple of two.
    orig_num = int(orig_num_file.readline().strip())
    # Line number of the remap_bam file (if single end data) or read pair
    # number if paired end data.
    line_num = 1
    
    # Index for walking through correct_maps.
    map_indx = 0
    # Number of correctly mapped reads for the current read (pair).
    correct = 0
    # Total number of correctly mapped read (pairs).
    total_correct = 0
    end_of_file = False
   
    # The idea here is that we will walk through remap_bam and check to see
    # whether all of the possible reads spawned from the original read (pair)
    # were mapped consistently with the original read (pair). Since correct_maps
    # is sorted, it should look something like [1, 1, 2, 2, 3, 3, 3, 3, ...].
    # If this example is for paired end data, this would mean that two reads
    # from the first pair mapped correctly, two reads from the second read pair
    # mapped correctly, four reads from the third read pair mapped correctly,
    # and so on. If the number of reads that mapped correctly for a given ID
    # number is equal to the corresponding orig_num, then everything is
    # consistent and we keep the read pair. For instance, in the example above,
    # correct_maps starts out with [1, 1]. This means that we have two correctly
    # mapped reads for the read pair with ID number 1 (the first read pair). If
    # the number on the first line of orig_num_file is 2, then we keep the read
    # pair because both reads mapped correctly. If the first number of
    # orig_num_file was 4, then we wouldn't keep the first read pair because
    # two of the reads didn't map correctly.
    while (not end_of_file and 
           (map_indx < len(correct_maps)) and 
           (line_num <= correct_maps[-1])):
        if line_num != correct_maps[map_indx]:
            # If we saw the correct number of remaps for the last read, we can
            # keep it.
            if correct == orig_num:
                total_correct += 1
                keep_bam.write(orig_read)
                # If the data is paired end, write out the paired read.
                if is_paired_end:
                    try:
                        orig_read = to_remap_bam.next()
                    except:
                        sys.stderr.write("File ended unexpectedly (no pair found).")
                        exit()
                    keep_bam.write(orig_read)
            else:
                try:
                    second_read = to_remap_bam.next()
                except:
                    end_of_file=True
                    break

            try:
                orig_read = to_remap_bam.next()
                orig_num = int(orig_num_file.readline().strip())
            except StopIteration:
                end_of_file = True
            line_num += 1
            correct = 0
        else:
            correct += 1
            map_indx += 1

    if correct == orig_num:
        total_correct += 1
        keep_bam.write(orig_read)
        # If the data is paired end, write out the paired read.
        if is_paired_end:
            try:
                orig_read = to_remap_bam.next()
            except:
                sys.stderr.write("File ended unexpectedly (no pair found).")
                exit()
            keep_bam.write(orig_read)
    sys.stderr.write('%d read (pair)s remapped to the correct position\n' %
                     total_correct)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', dest='is_paired_end', 
                        default=False, help=('Indicates that reads are '
                                             'paired-end (default is single).'))
    h = ('to.remap.bam file from find_intersecting_snps.py.')
    parser.add_argument("to_remap_bam", help=h)
    parser.add_argument("remap_bam", help='Remapped bam file.')
    parser.add_argument("keep_bam", help=('File to write correctly remapped '
                                          'reads to.'))
    h = 'to.remap.num.gz file from find_intersecting_snps.py.'
    parser.add_argument("orig_num_file", help=h)
    
    options = parser.parse_args()
    
    run(options.to_remap_bam, options.remap_bam, options.keep_bam,
        options.orig_num_file, options.is_paired_end)

if __name__ == '__main__':
    main()
