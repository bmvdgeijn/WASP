
import os
import argparse
import sys

import numpy
import tables

import chromosome
import util

    

def parse_args():
    parser = argparse.ArgumentParser(description="Creates an HDF5 file "
                                     "containing a table of chromosomes "
                                     "with information such as name and "
                                     "length. The data are parsed "
                                     "from a chromInfo.txt file, which "
                                     "can be downloaded from the UCSC "
                                     "genome browser")

    
    parser.add_argument("--chrom_info", required=True,
                        help="chromInfo.txt file from UCSC genome browser "
                        "database that contains list of chromosomes for "
                        "the genome assembly. For example chromInfo.txt.gz for hg19 can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/")

    
    parser.add_argument("--chrom_h5", required=True,
                        help="name of output HDF5 file")

    return parser.parse_args()


                

def convert_to_h5(input_chrom_filename, output_h5_filename):
    if os.path.exists(output_h5_filename):
        raise IOError("output file '%s' already exists\n" %
                      output_h5_filename)
                      
    chrom_tab = chromosome.ChromTab(output_h5_filename, "w")
    chrom_tab.load_from_file(input_chrom_filename)
    chrom_tab.close()



if __name__ == "__main__":
    args = parse_args()
    convert_to_h5(args.chrom_info, args.chrom_h5)


