import argparse
import sys

import chromosome
import util
import vcf
    

def parse_args():
    parser = argparse.ArgumentParser(description="extracts sample "
                                     "identifiers from a VCF file and writes "
                                     "them to stdout")
    
    parser.add_argument("vcf", help="VCF file containing sample identifiers")
                        
                        
    return parser.parse_args()

    

if __name__ == "__main__":
    args = parse_args()

    f = util.check_open(args.vcf)
    header, sample_ids = vcf.read_vcf_header(f)
    f.close()

    for samp in sample_ids:
        print samp
    
    
    


