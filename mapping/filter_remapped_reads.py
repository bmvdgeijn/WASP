def run(orig_bam, remap_bam, keep_bam, orig_num_file):
    import gzip
    import sys

    import pysam

    orig_bam = pysam.Samfile(orig_bam, "rb")
    remap_bam = pysam.Samfile(remap_bam, "rb")
    keep_bam = pysam.Samfile(keep_bam, "wb", template=orig_bam)
    orig_num_file = gzip.open(orig_num_file)

    correct_maps = []
    end_of_file = False
    
    # Get a list of reads that remapped correctly.
    remap_read = remap_bam.next()
    
    while not end_of_file:    
        chrm = remap_read.qname.strip().split(":")[1]
        if remap_read.is_reverse:
            pos = int(remap_read.qname.strip().split(":")[3])
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
    
    sys.stderr.write(str(len(correct_maps)) + 
                     " reads remapped to the correct position\n")
    
    # Pull out original aligned reads if all of the alternatives mapped
    # correctly.
    orig_read = orig_bam.next()
    orig_num = int(orig_num_file.readline().strip())
    line_num = 1
    
    map_indx = 0
    correct = 0
    end_of_file = False
    
    while (not end_of_file and map_indx < len(correct_maps) and 
           line_num <= correct_maps[-1]):
        if line_num < correct_maps[map_indx]:
            if orig_num == correct:
                keep_bam.write(orig_read)
            if options.is_paired_end:
                try:
                    orig_read = orig_bam.next()
                except:
                    sys.stderr.write("File ended unexpectedly (no pair found)")
                    exit()
                if orig_num == correct:
                    keep_bam.write(orig_read)
            
            line_num += 1
            correct = 0
            try:
                orig_read = orig_bam.next()
                orig_num = int(orig_num_file.readline().strip())
            except:
                end_of_file = True
        elif line_num == correct_maps[map_indx]:
            correct += 1
            map_indx += 1
        else:
            sys.stderr.write("There was a problem with the index sorting\n")
            exit()

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', dest='is_paired_end', default=False)
    parser.add_argument("orig_bam")
    parser.add_argument("remap_bam")
    parser.add_argument("keep_bam")
    parser.add_argument("orig_num_file")
    
    options= parser.parse_args()
    
    run(orig_bam, remap_bam, keep_bam, orig_num_file)

if __name__ == '__main__':
    main()
