
import os
import argparse
import sys

import numpy
import tables

import chromosome
import util
import vcf
    

def parse_args():
    parser = argparse.ArgumentParser(description="Read SNP information "
                                     "from a VCF file and store it "
                                     "in HDF5 files in the designated "
                                     "output directory")

    parser.add_argument("--chrom_h5", required=True,
                        help="path to HDF5 file containing chromosome info")
    
    parser.add_argument("--vcf", required=True, nargs="+",
                        help="VCF files containing SNPs information for each "
                        "chromosome. The chromosome name should be "
                        "part of the filename.")

    parser.add_argument("--snp_h5", required=True,
                        help="name of HDF5 file to write SNP data to")
                        
                        
    return parser.parse_args()



def convert_to_h5(vcf_files, chrom_h5_filename, snp_h5_filename):
    if os.path.exists(snp_h5_filename):
        raise IOError("output file '%s' already exists\n" %
                      snp_h5_filename)

    chrom_tab = chromosome.ChromTab(chrom_h5_filename)
    vcf.import_snps(chrom_tab, vcf_files, snp_h5_filename)
    

if __name__ == "__main__":
    args = parse_args()
    convert_to_h5(args.vcf, args.chrom_h5, args.snp_h5)



