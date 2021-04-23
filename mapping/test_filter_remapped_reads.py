import sys
import os
import subprocess

import filter_remapped_reads
import util
#
# filter_remapped_reads.py
#  INPUT FILES: 
#   to_remap_bam - input BAM file containing original set of reads
#                  that need to be remapped after having their alleles flipped
#
#   remap_bam - input BAM file containing remapped reads. Read names in this
#               file should be delimited with the '.' character and 
#               contain the following fields:
#                  <orig_name>.<coordinate>.<read_number>.<total_read_number>
#
#               For single-end reads <coordinate> is the left end of the read
#               (e.g. 16052611)
#               For paired-end reads the coordinate is the start of the 
#               the left read and start of the right read:
#               (e.g. 16052611-16052734)
#               
#
#
# OUTPUT FILES:
#   keep_bam - ouput BAM file containing reads that are retained
#              after filtering
#          
#

# DATA
to_remap_sam_lines = [
    "SRR1658224.34085432	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "SRR1658224.34085433	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085433	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "SRR1658224.34085434	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085434	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "SRR1658224.34975561	99	chr22	16071944	12	101M	=	16072163	320	ATTTATTTATTTATTTATTATTGGGACAGAGTCTCACTCTGTCCCCCAGACTGGAGTCCAGTGACATGATCTCAGCTCACTGCAACCTCTGCCTCGTGGGT	CCCFFFFFHHHHHJJJJJJJJJJJJIJJJJIEHIJJJJJJJIIJJJJJIJJJJJJJJJJIJHIJIJJJJIJJJJJHHHHHHFFFFFECEEEEDDDDDDBBD	AS:i:-5	XS:i:-22	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:89C11	YS:i:0	YT:Z:CP",
    "SRR1658224.34975561	147	chr22	16072163	12	101M	=	16071944	-320	GTCTCAAACTTCTGACCTCAGGTGATCCACCCACCTCGACCTCCCAAAGTGCTGGGATTACAGGCACTAGGTCCCTAAATTAGAGCCATATTCTTTAATGT	DDBCDEDCDCCDCC?DDDDDDDBACBDA<FFB:6HIIJIIJIIJJJJJJJJJJJJIJJIHJJJJJIJJJJJJJJJJJJJJJJJJJJJJHHHGGFFFFFCCC	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-5	YT:Z:CP",
    "SRR1658224.7462188	163	chr22	16235410	17	101M	=	16235625	316	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CC@FFFFFHHHHHJJJJJJJJJJJJJJJJIJBGIJJJJJJJJJJJJJIJIFIJJJJJJJJJHHHHGFFFFFFEEEEDEEDDDDDEED@CFFFEDDD?ABB?	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-5	YT:Z:CP",
    "SRR1658224.7462188	83	chr22	16235625	17	101M	=	16235410	-316	TTCAAAAGATGGTATATGCATTAATATTTTCATACAACTTCCAGCTTTTGTTTTTCTTCATTTAATTTTATTTATTTATTTATTTTTGAGATGGAGTCTCG	CBDDDDECEEDEFFFDFFFHHHHHHHJJIIJJIHIHFHGHJJJJJJJGJJJJJIJJJIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHHHHHFFFDFCCC	AS:i:-5	XS:i:-39	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:15G85	YS:i:0	YT:Z:CP",
    "SRR1658224.31153145	163	chr22	16235410	17	101M	=	16235625	316	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJJIJFHIJJJJJJJJJJJIJIJJFHIJJJJJJJJHHHHHFFFFFFEDEEEEEDDDDDEED@DEEEEDDDDDDB2	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-2	YT:Z:CP",
    "SRR1658224.31153145	83	chr22	16235625	17	101M	=	16235410	-316	TTCAAAAGATGGTATGTGCATTAATATTTTCATACAACTTCCAGTTTTTGTTTTTCTTCATTTAATTTTATTTATTTATTTATTTTTGAGATGGAGTCTCG	DDDDDDDDEEEEEEFFFFFFHHHHGHHJJIJJJIIJIJIHJHF@(JJJJJJJJJJJJIIIIJJJJJJJIJJJJJJJJJJJJJJJJJJJHHHHHFFFDFCCC	AS:i:-2	XS:i:-36	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:44C56	YS:i:0	YT:Z:CP",
    "SRR1658224.25014179	163	chr22	16236979	31	101M	=	16237137	259	ATGTTTTTTAAGATTTAATATTACTTTTTCCAACATCTTTTTATCCTCAAGTTTTTTATATTCCTGTTGTATTTTTTTATAGATAATAACTCCTGTTGAAT	CCCFFFFFHHHHFIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHGIJJJJJJJJIJJJJJJJHHHHHHHDCDDECDEEDDEDDDDDDDDDDCDC	AS:i:0	XS:i:-28	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:0	YT:Z:CP",
    "SRR1658224.25014179	83	chr22	16237137	31	101M	=	16236979	-259	TCATCGAACTACATTAATAAAATAATATAGCTTGATAATGAAGTAGGCTGAGAATAATCTCATACAAAACCAATAACAAATTTTGAAATACATTTACTTGC	CEFFFFFHHHHHHHHJJJJJJJJJIHJIJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIJJJIHJJJJJJIJJJJJJJJJJJJHHHHHFDDFFCCC	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:0	YT:Z:CP",
    "readpair1	163	chr22	100	12	101M	=	200	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair2	163	chr22	150	12	101M	=	250	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair1	83	chr22	200	12	101M	=	100	-201	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "readpair2	163	chr22	250	12	101M	=	150	-201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP"
]
remap_sam_lines = [
    # Read pair expected to map 2 times and maps to correct location 2 times
    "SRR1658224.34085432.16052611-16052734.1.2	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.1.2	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.2.2	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.2.2	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    # Read pair expected to map 2 times, but only maps 1 time
    "SRR1658224.34975561.16071944-16072163.2.2	99	chr22	16071944	12	101M	=	16072163	320	ATTTATTTATTTATTTATTATTGGGACAGAGTCTCACTCTGTCCCCCAGACTGGAGTCCAGTGACATGATCTCAGCTCACTGCAACCTCTGCCTCGTGGGT	CCCFFFFFHHHHHJJJJJJJJJJJJIJJJJIEHIJJJJJJJIIJJJJJIJJJJJJJJJJIJHIJIJJJJIJJJJJHHHHHHFFFFFECEEEEDDDDDDBBD	AS:i:-5	XS:i:-22	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:89C11	YS:i:0	YT:Z:CP",
    "SRR1658224.34975561.16071944-16072163.2.2	147	chr22	16072163	12	101M	=	16071944	-320	GTCTCAAACTTCTGACCTCAGGTGATCCACCCACCTCGACCTCCCAAAGTGCTGGGATTACAGGCACTAGGTCCCTAAATTAGAGCCATATTCTTTAATGT	DDBCDEDCDCCDCC?DDDDDDDBACBDA<FFB:6HIIJIIJIIJJJJJJJJJJJJIJJIHJJJJJIJJJJJJJJJJJJJJJJJJJJJJHHHGGFFFFFCCC	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-5	YT:Z:CP",
    # Read pair expected to map 2 times, but only 1/2 of 2nd pair maps back to same location
    "SRR1658224.7462188.16235410-16235625.1.2	163	chr22	16235410	17	101M	=	16235625	316	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CC@FFFFFHHHHHJJJJJJJJJJJJJJJJIJBGIJJJJJJJJJJJJJIJIFIJJJJJJJJJHHHHGFFFFFFEEEEDEEDDDDDEED@CFFFEDDD?ABB?	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-5	YT:Z:CP",
    "SRR1658224.7462188.16235410-16235625.1.2	83	chr22	16235625	17	101M	=	16235410	-316	TTCAAAAGATGGTATATGCATTAATATTTTCATACAACTTCCAGCTTTTGTTTTTCTTCATTTAATTTTATTTATTTATTTATTTTTGAGATGGAGTCTCG	CBDDDDECEEDEFFFDFFFHHHHHHHJJIIJJIHIHFHGHJJJJJJJGJJJJJIJJJIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHHHHHFFFDFCCC	AS:i:-5	XS:i:-39	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:15G85	YS:i:0	YT:Z:CP",
    "SRR1658224.7462188.16235410-16235625.2.2	163	chr22	16235410	17	101M	*	0	0	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CC@FFFFFHHHHHJJJJJJJJJJJJJJJJIJBGIJJJJJJJJJJJJJIJIFIJJJJJJJJJHHHHGFFFFFFEEEEDEEDDDDDEED@CFFFEDDD?ABB?	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-5	YT:Z:CP",
    # Read pair expected to map 2 times, but 1 pair maps to wrong location
    "SRR1658224.31153145.16235410-16235625.1.2	163	chr22	16235410	17	101M	=	16235625	316	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJJIJFHIJJJJJJJJJJJIJIJJFHIJJJJJJJJHHHHHFFFFFFEDEEEEEDDDDDEED@DEEEEDDDDDDB2	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-2	YT:Z:CP",
    "SRR1658224.31153145.16235410-16235625.1.2	83	chr22	16235625	17	101M	=	16235410	-316	TTCAAAAGATGGTATGTGCATTAATATTTTCATACAACTTCCAGTTTTTGTTTTTCTTCATTTAATTTTATTTATTTATTTATTTTTGAGATGGAGTCTCG	DDDDDDDDEEEEEEFFFFFFHHHHGHHJJIJJJIIJIJIHJHF@(JJJJJJJJJJJJIIIIJJJJJJJIJJJJJJJJJJJJJJJJJJJHHHHHFFFDFCCC	AS:i:-2	XS:i:-36	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:44C56	YS:i:0	YT:Z:CP",
    "SRR1658224.31153145.16235410-16235625.2.2	163	chr22	18235410	17	101M	=	16235625	316	AGATAATTGTCTTATTTTTTTAAAAAAAGAGTAACTTTATATTATGGAATTCATAATATTTGAGACTATAATGCATGACATAAATAGTATAAAGGAGAGAG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJJIJFHIJJJJJJJJJJJIJIJJFHIJJJJJJJJHHHHHFFFFFFEDEEEEEDDDDDEED@DEEEEDDDDDDB2	AS:i:0	XS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-2	YT:Z:CP",
    "SRR1658224.31153145.16235410-16235625.2.2	83	chr22	18235625	17	101M	=	16235410	-316	TTCAAAAGATGGTATGTGCATTAATATTTTCATACAACTTCCAGTTTTTGTTTTTCTTCATTTAATTTTATTTATTTATTTATTTTTGAGATGGAGTCTCG	DDDDDDDDEEEEEEFFFFFFHHHHGHHJJIJJJIIJIJIHJHF@(JJJJJJJJJJJJIIIIJJJJJJJIJJJJJJJJJJJJJJJJJJJHHHHHFFFDFCCC	AS:i:-2	XS:i:-36	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:44C56	YS:i:0	YT:Z:CP",
    # Read pair expected to map 2 times, but does not map at all
    # "SRR1658224.25014179"
    # Read pairs expected to map 1 times, with read-pairs interleaved
    "readpair1.100-200.1.2	163	chr22	100	12	101M	=	200	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair2.150-250.1.2	163	chr22	150	12	101M	=	250	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair1.100-200.1.2	83	chr22	200	12	101M	=	100	-201	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",        
    "readpair2.150-250.1.2	163	chr22	250	12	101M	=	150	-201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair1.100-200.2.2	163	chr22	100	12	101M	=	200	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair2.150-250.2.2	163	chr22	150	12	101M	=	250	201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "readpair1.100-200.2.2	83	chr22	200	12	101M	=	100	-201	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "readpair2.150-250.2.2	163	chr22	250	12	101M	=	150	-201	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    # Read pair is secondary
    "SRR1658224.34085433.16052611-16052734.1.1	419	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085433.16052611-16052734.1.1	339	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    # Read pair is supplementary
    "SRR1658224.34085434.16052611-16052734.1.1	2211	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085434.16052611-16052734.1.1	2131	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP"

]


to_remap_sam_lines_single= [
"single1	162	chr22	250	12	101M	*	0	0	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP"]
    
remap_sam_lines_single = [
 "single1.250.1.1	162	chr22	250	12	101M	*	0	0	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP"
]


#
# TODO: need to verify that interleaved read pairs handled appropriately
# TODO: need to test single end reads
#
#


def write_sam_header(f):
    f.write("@HD	VN:1.0	SO:coordinate\n")
    f.write("@SQ	SN:chr22	LN:51304566\n")
    f.write('@PG	ID:bowtie2	PN:bowtie2	VN:2.2.6	CL:"/iblm/netapp/home/gmcvicker/anaconda2/bin/bowtie2-align-s --wrapper basic-0 -x /iblm/netapp/data1/external/GRC37/combined/bowtie2_index/hg37 -1 /tmp/16686.inpipe1 -2 /tmp/16686.inpipe2\n')




def write_to_remap_bam(sam_lines, data_dir="test_data", bam_filename="test_data/test.to.remap.bam"):

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


    
def write_remap_bam(sam_lines, data_dir="test_data", bam_filename="test_data/test.remap.bam"):

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        
    # write temporary file in SAM format, before converting to BAM
    sam_filename = data_dir + "/tmp.sam"
    f = open(sam_filename, "wt")
    write_sam_header(f)
    for line in sam_lines:
        f.write(line + "\n")
    f.close()

    # write to temp bam file
    tmp_bam_filename = data_dir + "/tmp.bam"
    subprocess.check_call("samtools view -b %s > %s" % (sam_filename, tmp_bam_filename), shell=True)
    # sort the temp bam file
    util.sort_bam(tmp_bam_filename, data_dir + "/tmp")
    # remove temp bam
    os.remove(tmp_bam_filename)
    # rename sorted bam to output bam filename
    os.rename(data_dir + "/tmp.sort.bam", bam_filename)

    
def read_bam(bam):
    """
    Read a bam file into a list where each element of the list is a line from
    the bam file (with the newline stripped). The header is discarded.
    """
    lines = []
    res = subprocess.check_output('samtools view %s' % bam, shell=True)
    if res:
        lines += res.decode("utf-8").strip().split('\n')
    return lines



#
# Test single-end reads (Added 4/22/2021)
# TODO: current tests for single-end reads are very limited,
# and just check that remapped read is retained when it maps
# back to same position with same CIGAR. More should be added.
#
def test_filter_remapped_reads_single():
    test_dir = "test_data"
    to_remap_bam_filename = "test_data/test.to.remap.single.bam"
    remap_bam_filename = "test_data/test.remap.single.bam"
    keep_bam_filename = "test_data/keep.single.bam"

    # write test input data
    write_to_remap_bam(
        sam_lines=to_remap_sam_lines_single,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=remap_sam_lines_single,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )

    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )

    # read in filtered reads
    lines = read_bam(keep_bam_filename)

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
    read_name = "single1"

    sys.stderr.write("%s\n" % repr(read_dict))
    
    assert read_name in read_dict
    reads = read_dict[read_name]
    assert len(reads) == 1


    
#
# test paired-end reads
#
def test_filter_remapped_reads_pe():
    test_dir = "test_data"
    to_remap_bam_filename = "test_data/test.to.remap.bam"
    remap_bam_filename = "test_data/test.remap.bam"
    keep_bam_filename = "test_data/keep.bam"

    # write test input data
    write_to_remap_bam(
        sam_lines=to_remap_sam_lines,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=remap_sam_lines,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )

    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )

    # read in filtered reads
    lines = read_bam(keep_bam_filename)

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

    # we expect a read pair with this identifier:
    read_name = "SRR1658224.34085432"
    assert read_name in read_dict
    reads = read_dict[read_name]
    assert len(reads) == 2

    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 16052611
    assert pos2 == 16052734

    # expect these read pairs to be filtered out (not present)
    # only one version of read pair maps (expect 2)
    assert "SRR1658224.34975561" not in read_dict

    # 1/2 of second read pair missing
    assert "SRR1658224.7462188" not in read_dict

    # 1 pair maps to wrong location
    assert "SRR1658224.31153145" not in read_dict

    # neither pair maps
    assert "SRR1658224.25014179" not in read_dict

    # expect these (interleaved) read pairs to be kept
    read_name = "readpair1"
    assert read_name in read_dict
    reads = read_dict[read_name]
    assert len(reads) == 2
    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 100
    assert pos2 == 200

    sys.stderr.write("\n\nread_dict: %s\n\\n" % repr(read_dict))
    
    read_name = "readpair2"
    assert read_name in read_dict
    reads = read_dict[read_name]
    assert len(reads) == 2
    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 150
    assert pos2 == 250

    # secondary alignment
    assert "SRR1658224.34085433" not in read_dict

    # supplementary alignment
    assert "SRR1658224.34085434" not in read_dict


# CIGAR TESTING
# If reads have different CIGARs after the 2nd mapping,
# they should be discarded

# define data
to_remap_CIGAR_sam_lines = [
    "SRR1658224.34085432	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP"
]
remap_CIGAR_sam_lines = [
    # Read pair expected to map 2 times and maps to correct location 2 times
    "SRR1658224.34085432.16052611-16052734.1.2	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.1.2	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.2.2	163	chr22	16052611	12	101M	=	16052734	224	TGGAGACATAAAATGAGGCATATCTGACCTCCACTTCCAAAAACATCTGAGATAGGTCTCAGTTAATTAAGAAAGTTTGTTCTGCCTAGTTTAAGGACATG	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJIJJIHIJJJJEHIJJJHJJJJJJJJJJJJ=DHHHHHFFFFFFEEEEEEDDCDDDC	AS:i:-11	XS:i:-17	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:7G44C48	YS:i:0	YT:Z:CP",
    "SRR1658224.34085432.16052611-16052734.2.2	83	chr22	16052734	12	101M	=	16052611	-224	TCCTGACAGCATGTGCCCAAGGTGGTCAGGATACAGCTTGCTTCTATATATTTTAGGGAGAAAATACATCAGCCTGTAAACAAAAAATTAAATTCTAAGGT	DDDDDDDDDDDDDDEDEEEFFFFHHFHHIIFIIJJJJIJJJJJJJJJJIIJJJIIIJIJIJJJJIFIIIJJIJJJJJJJIIJJJJJJJHHHHHFFFFFCCC	AS:i:0	XS:i:-12	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:101	YS:i:-11	YT:Z:CP"
]

def test_filter_different_CIGAR():
    """test whether reads that map to the same location but have different
    CIGAR flags are still appropriately discarded"""
    test_dir = "test_data"
    to_remap_bam_filename = "test_data/test.to.remap.bam"
    remap_bam_filename = "test_data/test.remap.bam"
    keep_bam_filename = "test_data/keep.bam"

    # write test input data
    write_to_remap_bam(
        sam_lines=to_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )
    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )
    # read in filtered reads
    lines = read_bam(keep_bam_filename)
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
    # we expect a read pair with this identifier:
    read_name = "SRR1658224.34085432"
    assert read_name in read_dict
    reads = read_dict[read_name]
    assert len(reads) == 2
    pos1 = int(reads[0][3])
    pos2 = int(reads[1][3])
    assert pos1 == 16052611
    assert pos2 == 16052734

    # we know now that these pairs do indeed map to the same location
    # but what if the CIGAR changes after the second mapping?
    # then the reads should be discarded
    new_remap_CIGAR_sam_lines = [
        read.replace("\t101M\t", "\t101M1D\t")
        for read in remap_CIGAR_sam_lines
    ]
    # write test input data
    write_to_remap_bam(
        sam_lines=to_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=new_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )
    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )
    # read in filtered reads
    lines = read_bam(keep_bam_filename)
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
    # we expect a read pair with this identifier:
    assert "SRR1658224.34085432" not in read_dict

    # now what if the CIGAR changes in only one read but not its pair?
    # then both reads should be discarded
    new_remap_CIGAR_sam_lines = remap_CIGAR_sam_lines
    new_remap_CIGAR_sam_lines[1] = new_remap_CIGAR_sam_lines[1].replace("\t101M\t", "\t101M1D\t")
    # write test input data
    write_to_remap_bam(
        sam_lines=to_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=new_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )
    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )
    # read in filtered reads
    lines = read_bam(keep_bam_filename)
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
    # we expect a read pair with this identifier:
    assert "SRR1658224.34085432" not in read_dict

    # now what if the CIGARs are different between pairs
    # originally but become the same later?
    # then both reads should be discarded
    #
    # We are no longer handling this weird case correctly...
    # issue is that read1 and read2 are no longer consistently
    # labeled (read1 sometimes becomes read2 and vice-versa,)
    # making it harder to determine if CIGARs change.
    #
    # new_to_remap_CIGAR_sam_lines = to_remap_CIGAR_sam_lines
    # new_to_remap_CIGAR_sam_lines[1] = new_to_remap_CIGAR_sam_lines[1].replace("\t101M\t", "\t101M1D\t")
    # # change this back to what it was
    # new_remap_CIGAR_sam_lines = remap_CIGAR_sam_lines
    # write_to_remap_bam_pe(
    #     sam_lines=new_to_remap_CIGAR_sam_lines,
    #     data_dir=test_dir,
    #     bam_filename=to_remap_bam_filename
    # )
    # write_remap_bam_pe(
    #     sam_lines=new_remap_CIGAR_sam_lines,
    #     data_dir=test_dir,
    #     bam_filename=remap_bam_filename
    # )

    # sys.stderr.write("TO REMAP:\n%s\n\n" % "\n".join(new_to_remap_CIGAR_sam_lines))
    # sys.stderr.write("REMAPPED:\n%s\n\n" % "\n".join(new_remap_CIGAR_sam_lines))

    
    # # run filter remapped reads
    # filter_remapped_reads.main(
    #     to_remap_bam_filename,
    #     remap_bam_filename,
    #     keep_bam_filename
    # )
    # # read in filtered reads
    # lines = read_bam(keep_bam_filename)
    # # read lines from keep BAM file
    # read_dict = {}
    # for line in lines:
    #     words = line.split()
    #     read_name = words[0]
    #     if read_name in read_dict:
    #         read_dict[read_name].append(words)
    #     else:
    #         read_dict[read_name] = [words]
    # # verify that filtered reads look correct
    # # we expect a read pair with this identifier to be absent:
    # assert "SRR1658224.34085432" not in read_dict

    # now what if the CIGARs are different between pairs
    # originally and stay the same as they were after the second remapping
    # in both alternative reads?
    # then both of the pairs should be kept
    new_to_remap_CIGAR_sam_lines = to_remap_CIGAR_sam_lines
    new_to_remap_CIGAR_sam_lines[1] = new_to_remap_CIGAR_sam_lines[1].replace("\t101M\t", "\t101M1D\t")
    # change this, as well
    new_remap_CIGAR_sam_lines = remap_CIGAR_sam_lines
    new_remap_CIGAR_sam_lines[1] = new_remap_CIGAR_sam_lines[1].replace("\t101M\t", "\t101M1D\t")
    new_remap_CIGAR_sam_lines[3] = new_remap_CIGAR_sam_lines[3].replace("\t101M\t", "\t101M1D\t")
    # write test input data
    write_to_remap_bam(
        sam_lines=new_to_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=to_remap_bam_filename
    )
    write_remap_bam(
        sam_lines=new_remap_CIGAR_sam_lines,
        data_dir=test_dir,
        bam_filename=remap_bam_filename
    )
    # run filter remapped reads
    filter_remapped_reads.main(
        to_remap_bam_filename,
        remap_bam_filename,
        keep_bam_filename
    )
    # read in filtered reads
    lines = read_bam(keep_bam_filename)
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
    # we expect a read pair with this identifier:
    assert "SRR1658224.34085432" in read_dict
