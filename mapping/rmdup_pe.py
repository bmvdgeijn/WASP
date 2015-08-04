from random import choice
import pysam
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help="input BAM or SAM file (must be sorted!)")
    parser.add_argument("output_bam", help="output BAM or SAM file")
    
    max_window=10000
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

    readf=Read_filter(infile,outfile,max_window)
    infile.close()
    outfile.close()

class Read_filter:
    def __init__(self,infile,outfile,max_window):
        self.read_table=[[] for x in range(max_window)]
        self.num_reads=0
        self.cur_pos=0
        self.chr=""
        self.cur_read=infile.next()
        self.infile=infile
        self.outfile=outfile
        self.max_window=max_window
        self.finished=False
        while not self.finished:
            self.fill_table()
            self.empty_slot()
        self.empty_table()

    def fill_table(self):
        if self.cur_read.rname != self.chr:
            self.empty_table()
            self.chr=self.cur_read.rname

        if self.num_reads==0:
            self.cur_pos=self.cur_read.pos
        while not self.finished and self.cur_read.rname==self.chr and self.cur_read.pos<self.cur_pos+self.max_window:
            self.read_table[self.cur_read.pos % self.max_window].append(self.cur_read)
            self.num_reads+=1
            try:
                self.cur_read=self.infile.next()
            except:
                self.finished=True
        
    def empty_table(self):
        while self.num_reads>0:
            self.empty_slot()

    def empty_slot(self):
        ends=dict()
        
        for read in self.read_table[self.cur_pos % self.max_window]:
            mate_pos=read.mpos
            if mate_pos in ends:
                ends[mate_pos].append(read)
            else:
                ends[mate_pos]=[read]
        
        self.num_reads-=len(self.read_table[self.cur_pos % self.max_window])
        self.read_table[self.cur_pos % self.max_window] = []
        
        for end_key in ends:
            read_list=ends[end_key]
            if read_list[0].mpos < self.cur_pos:
                continue
            
            while len(read_list)>0:
                keep_indx=choice(range(len(read_list)))
                keep_read=read_list.pop(keep_indx)
                read_name=keep_read.qname
                mate_pos=keep_read.mpos
                found=False
                for mate_read in self.read_table[mate_pos % self.max_window]:
                    #sys.stderr.write("%s\t%s\n" % (mate_read.qname,keep_read.qname))
                    if mate_read.qname==keep_read.qname:
                        self.outfile.write(keep_read)
                        self.outfile.write(mate_read)
                        found=True
                        break
        if self.num_reads>0:
            while len(self.read_table[self.cur_pos % self.max_window])==0:
                self.cur_pos+=1



main()
