import sys
import os
import subprocess

import filter_remapped_reads
import util
import rmdup_pe

#
# rmdump_pe.py <input_bam> <output_bam>
#


def write_sam_header(f):
    f.write("@HD	VN:1.0	SO:coordinate\n")
    f.write("@SQ	SN:chr22	LN:51304566\n")
    f.write('@PG	ID:bowtie2	PN:bowtie2	VN:2.2.6	CL:"/iblm/netapp/home/gmcvicker/anaconda2/bin/bowtie2-align-s --wrapper basic-0 -x /iblm/netapp/data1/external/GRC37/combined/bowtie2_index/hg37 -1 /tmp/16686.inpipe1 -2 /tmp/16686.inpipe2\n')



def write_bam_pe(data_dir="test_data", bam_filename="test_data/rmdup_input.bam"):
    sam_lines = ["readpair1	163	chr22	100	12	101M	=	200	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                 # include a mate-unmapped read pair
                 "unmapped_readpair1	163	chr22	100	12	101M	=	100	0	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",                 
                 # duplicate of first read pair 
                 "dup_readpair1	163	chr22	100	12	101M	=	200	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                 "readpair2	163	chr22	150	12	101M	=	250	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",

                 # readpair 3 has same first read, but different second read as readpair2
                 "readpair3	163	chr22	150	12	101M	=	251	202	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",

                 # readpair4 has same positions as readpair2
                 "readpair4	163	chr22	150	12	101M	=	250	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                 
                 "dup_readpair1	83	chr22	200	12	101M	=	100	-201	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
                 
                 "readpair1	83	chr22	200	12	101M	=	100	-201	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
                
                 "readpair2	163	chr22	250	12	101M	=	150	-201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",

                 "readpair4	163	chr22	250	12	101M	=	150	-201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",

                 "readpair3	163	chr22	251	12	101M	=	150	-202	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",                 
                 # couple of read pairs that are completely overlapping
                 # (i.e. at same position)
                 "readpair5	163	chr22	500	12	101M	=	500	-101	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                 
                 "readpair5	163	chr22	500	12	101M	=	500	-101	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                "dup_readpair5	163	chr22	500	12	101M	=	500	-101	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
                 "dup_readpair5	163	chr22	500	12	101M	=	500	-101	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",

                 # include a couple of unmapped reads
                 "unmapped_read1	77	*	0	0	*	*	0	0	ACTAGACATACATAACACATATACCCACCC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                 "unmapped_read1	141	*	0	0	*	*	0	0	ACTAGACATACATAACACATATACCCACCC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                 
#                 "unmapped_readpair1	141	chr22	0	0	*	=	*	*	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP"
                

                 
                 
    ]

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        
    # write temporary file in SAM format, before converting to BAM
    sam_filename = data_dir + "/tmp.sam"
    f = open(sam_filename, "w")
    write_sam_header(f)
    for line in sam_lines:
        f.write(line + "\n")
    f.close()

    subprocess.check_call("samtools view -b %s > %s" % (sam_filename, bam_filename), shell=True)

    
def read_bam(bam):
    """
    Read a bam file into a list where each element of the list is a line from
    the bam file (with the newline stripped). The header is discarded.
    """
    res = subprocess.check_output('samtools view %s' % bam, shell=True)
    return res.decode("utf-8").strip().split('\n')


def test_rmdup_pe():
    test_dir = "test_data"
    rmdup_input_bam = "test_data/rmdup_input.bam"
    rmdup_output_bam = "test_data/rmdup_output.bam"

    # write test input data
    write_bam_pe(data_dir=test_dir, bam_filename=rmdup_input_bam)

    # remove duplicates
    rmdup_pe.main(rmdup_input_bam, rmdup_output_bam)
    
    # read in filtered reads
    lines = read_bam(rmdup_output_bam)
    
    # read lines from keep BAM file
    read_dict = {}
    for line in lines:
        words = line.split()
        read_name = words[0]
        if read_name in read_dict:
            read_dict[read_name].append(words)
        else:
            read_dict[read_name] = [words]
                    
    # verify that filtered reads look correct

    # expect one of readpair1 and dup_readpair (but not both)
    # expect one of readpair2 and readpair4 (but not both)
    # expect readpair3
    # expect one of readpair5 and dup_readpair5 (but not both)

    assert len(read_dict) == 4
    
    # expect one of readpair1 and dup_readpair (but not both)
    if "readpair1" in read_dict:
        assert "dup_readpair1" not in read_dict
        reads = read_dict["readpair1"]
    else:
        assert "dup_readpair1" in read_dict
        reads = read_dict["dup_readpair1"]
        
    assert len(reads) == 2

    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 100
    assert pos2 == 200

    # expect readpair2 OR readpair4 to be present
    if "readpair4" in read_dict:
        assert "readpair2" not in read_dict
        reads = read_dict["readpair4"]
    else:
        assert "readpair2" in read_dict
        reads = read_dict["readpair2"]
        
    assert len(reads) == 2
    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 150
    assert pos2 == 250

    # expect readpair3 to be present
    assert "readpair3" in read_dict
    reads =read_dict["readpair3"]
    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 150
    assert pos2 == 251


    # expect readpair5 OR dup_readpair5 to be present
    if "readpair5" in read_dict:
        assert "dup_readpair5" not in read_dict
        reads = read_dict["readpair5"]
    else:
        assert "dup_readpair5" in read_dict
        reads = read_dict["dup_readpair5"]

    

    
        

    
    
