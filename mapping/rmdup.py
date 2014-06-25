from random import choice
import pysam
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_bam', help="input BAM or SAM file (must be sorted!)")
parser.add_argument("output_bam", help="output BAM or SAM file")

options = parser.parse_args()

if options.input_bam.endswith(".sam") or options.input_bam.endswith("sam.gz"):
    infile = pysam.Samfile(options.input_bam, "r")
else:
    # assume binary BAM file
    infile = pysam.Samfile(options.input_bam, "rb")

if options.output_bam.endswith(".sam"):
    # output in text SAM format
    outfile = pysam.Samfile(options.output_bam, "w", template=infile)
elif options.output_bam.endswith(".bam"):
    # output in binary compressed BAM format
    outfile = pysam.Samfile(options.output_bam, "wb", template=infile)
else:
    raise ValueError("name of output file must end with .bam or .sam")


chr=""
pos=0
linelistplus=[]
linelistminus=[]
for line in infile:
    if line.rname!=chr or line.pos!=pos:
        if(len(linelistplus)>0):
            outfile.write(choice(linelistplus))
        if(len(linelistminus)>0):
            outfile.write(choice(linelistminus))
        chr=line.rname
        pos=line.pos
        linelistplus=[]
        linelistminus=[]
    if(line.flag==0):
        linelistplus.append(line)
    if(line.flag==16):
        linelistminus.append(line)
infile.close()
outfile.close()
