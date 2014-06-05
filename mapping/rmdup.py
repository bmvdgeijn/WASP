from random import choice
import pysam
import os
import sys
import optparse


parser = optparse.OptionParser()
parser.add_option("-i", action='store_true', dest='input_sam', default=False)
parser.add_option("-o", action='store_true', dest='output_sam', default=False)

options,args= parser.parse_args()

if options.input_sam:
    infile = pysam.Samfile(args[0],'r')
else:
    infile = pysam.Samfile(args[0],'rb')
    
if options.output_sam:
    outfile = pysam.Samfile(args[1],'w',template=infile)
else:
    outfile = pysam.Samfile(args[1],'wb',template=infile)

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
