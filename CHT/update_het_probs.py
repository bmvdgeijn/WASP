import gzip
import argparse
import math
import sys
from argparse import ArgumentParser

import tables

import util


def parse_options():
    
    parser = ArgumentParser(description="""This script adjusts
            heterozygote probabilities in CHT files to account for 
            possible genotyping errors. Total counts of reference 
            and alternative alleles are used to adjust the 
            probability. The read counts that are provided can be
            from the same experiment, combined across many different
            experiments, or (perhaps ideally) from DNA sequencing
            of the same individual.""")

    parser.add_argument("infile", action='store', 
                        help="Input file for Combined Haplotype "
                        "Test (CHT) that needs het probabilities "
                        "adjusted",
                        default=None)
    
    parser.add_argument("outfile", action='store',
                        help="Output CHT file with heterozygote "
                        "probabilities adjusted",
                        default=None)
    
    parser.add_argument("--ref_as_counts", 
                        action='store',
                        help="Path to HDF5 file containing counts "
                        "of reads that match reference allele",
                        metavar="REF_AS_COUNT_H5_FILE",
                        required=True)
    
    parser.add_argument("--alt_as_counts",
                        help="Path to HDF5 file containing counts "
                        "of reads that match alternate allele",
                        metavar="ALT_AS_COUNT_H5_FILE",
                        action='store', required=True)
    
    return parser.parse_args()



def main():
    error = 0.01
    args = parse_options()

    if util.is_gzipped(args.infile):
        infile = gzip.open(args.infile, "rt")
    else:
        infile = open(args.infile, "r")
        
    if args.outfile.endswith(".gz"):
        outfile = gzip.open(args.outfile,"w")
    else:
        outfile = open(args.outfile,"w")

    ref_count_h5 = tables.openFile(args.ref_as_counts)
    alt_count_h5 = tables.openFile(args.alt_as_counts)

    snp_line = infile.readline()
    if snp_line:
        outfile.write(snp_line)
    else:
        sys.stderr.write("The input file was empty.\n")
        exit(-1)

    snp_line = infile.readline()
    while snp_line:
        snpinfo = snp_line.strip().split()
        if snpinfo[9] == "NA":
            outfile.write(snp_line)
        else:
            new_hetps = process_one_snp(snpinfo, ref_count_h5, 
                                        alt_count_h5, error)
            outfile.write("\t".join(snpinfo[:10] + 
                                    [";".join(new_hetps)] + 
                                    snpinfo[11:]) + "\n")
        snp_line = infile.readline()

    ref_count_h5.close()
    alt_count_h5.close()
    

def process_one_snp(snpinfo, ref_count_h5, alt_count_h5, error):
    chrm = snpinfo[0]
    # positions of target SNPs
    snplocs = [int(y.strip()) for y in snpinfo[9].split(';')]
    
    # heterozygote probabilities of target SNPs
    hetps = [float(y.strip()) for y in snpinfo[10].split(';')]
    update_hetps = []

    ref_node = ref_count_h5.getNode("/%s" % chrm)
    alt_node = alt_count_h5.getNode("/%s" % chrm)
    
    for i in range(len(snplocs)):
        pos = snplocs[i]
        adr = ref_node[pos-1]
        ada = alt_node[pos-1]
        update_hetps.append(str(get_posterior_hetp(hetps[i], adr, 
                                                   ada, error)))
    return update_hetps


def get_posterior_hetp(hetp_prior, adr, ada, error):
    prior = min(0.99, hetp_prior)
    badlike = addlogs(math.log(error)*adr + 
                      math.log(1-error)*ada,
                      math.log(1-error)*adr +
                      math.log(error)*ada)
    goodlike = math.log(0.5)*adr + math.log(0.5)*ada
    if goodlike-badlike > 40:
        # avoid overflow (very close to 1.0)
        return 1.0
    else:
        return prior*math.exp(goodlike - badlike) / (prior*math.exp(goodlike - badlike) + (1.0 - prior))

    
def addlogs(loga, logb):
    return max(loga, logb) + math.log(1+math.exp(-abs(loga-logb)))


main()
