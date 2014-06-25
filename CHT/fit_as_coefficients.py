# Copyright 2013 Graham McVicker and Bryce van de Geijn
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os
import math
import time
import gzip
import argparse

from scipy.optimize import *
from scipy import cast
from scipy.special import gammaln
from scipy.special import betaln
import scipy.stats

import numpy as np
from random import shuffle
from random import randint

def open_input_files(in_filename):
    if not os.path.exists(in_filename) or not os.path.isfile(in_filename):
        sys.stderr.write("input file %s does not exist or is not a regular file\n" %
                         in_filename)
        exit(2)
    
    # read file that contains list of input files
    in_file = open(in_filename)

    infiles = []
    for line in in_file:
        # open each input file and read first line
        filename = line.rstrip()
        if not filename or not os.path.exists(filename) or not os.path.isfile(filename):
            sys.stderr.write("input file '%s' does not exist or is not a regular file\n" 
                             % line)
            exit(2)
        if filename.endswith(".gz"):
            f = gzip.open(filename)
        else:
            f = open(filename)
            
        # skip header
        f.readline()

        infiles.append(f)
    in_file.close()
    
    if len(infiles) == 0:
        sys.stderr.write("no input files specified in file '%s'\n" % options.infile_list)
        exit(2)

    return infiles
    
def main():
    options = parse_options()
    infiles = open_input_files(options.infile_list)
    outfile = open(options.out_file,"w")

    for cur_file in infiles:
        AS_ref=[]
        AS_alt=[]
        hetps=[]
        cur_line=cur_file.readline()
        while True:
            try:
                cur_line=cur_file.readline()
            except:
                break
            if not cur_line:
                break
            snpinfo=cur_line.strip().split()

            if snpinfo[12] != "NA":
                AS_ref = AS_ref + [int(y) for y in snpinfo[12].split(';')]
                AS_alt = AS_alt + [int(y) for y in snpinfo[13].split(';')]
                hetps = hetps + [float(y.strip()) for y in snpinfo[10].split(';')]
        dispersion=fmin(likelihood, 0.001, args=(AS_ref, AS_alt, hetps, options.read_error_rate))
        outfile.write(str(dispersion[0])+"\n")
        outfile.flush()

def likelihood(dispersion,AS_ref,AS_alt,hetps,error):
            cur_like=0
            for i in range(len(AS_ref)):
                cur_like=cur_like+AS_betabinom_loglike([math.log(0.5),math.log(0.5)], dispersion, AS_ref[i], AS_alt[i], hetps[i], error)
            return -cur_like
                
def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("-e", action='store', dest='read_error_rate',
                        help="estimate of error rate, used to update "
                        "heterozygous genotype probabilities "
                        "(currently this option disabled / not used)",
                        type=float, default=0.005)
    parser.add_argument("infile_list", action='store', default=None)
    parser.add_argument("out_file", action='store', default=None)

    return parser.parse_args()
                
def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))



#Given parameters, returns log likelihood.  Note that some parts have been cancelled out
def AS_betabinom_loglike(logps, sigma, AS1, AS2, hetp, error):
    a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
    b = math.exp(logps[1] + math.log(1/sigma**2 - 1))
    
    part1 = 0
    part1 += betaln(AS1 + a, AS2 + b)
    part1 -= betaln(a, b)
    
    if hetp==1:
        return part1        

    e1 = math.log(error) * AS1 + math.log(1 - error) * AS2
    e2 = math.log(error) * AS2 + math.log(1 - error) * AS1
    if hetp == 0:
        return addlogs(e1, e2)
    
    return addlogs(math.log(hetp)+part1, math.log(1-hetp) + addlogs(e1,e2))

def parse_test_snp(snpinfo, options):
    snp_id = snpinfo[2]
    if snpinfo[16] == "NA":
        # SNP is missing data
        tot = 0
    else:
        # rescale these to put totals in reasonable range
        # better approach might be to divide by minimum total
        # across individuals
        #if tot>10000:
        tot = float(snpinfo[16]) #/1000000

    if snpinfo[6] == "NA":
        geno_hap1 = 0
        geno_hap2 = 0
    else:
        geno_hap1 = int(snpinfo[6].strip().split("|")[0])
        geno_hap2 = int(snpinfo[6].strip().split("|")[1])
    
    if snpinfo[15] == "NA":
        count = 0
    else:
        count = int(snpinfo[15])

    if snpinfo[9].strip() == "NA" or geno_hap1 == geno_hap2:
        # SNP is homozygous, so there is no AS info
        return TestSNP(snp_id, geno_hap1, geno_hap2, [], [], [], tot, count)    
    else:
        # positions of target SNPs (not currently used)
        snplocs=[int(y.strip()) for y in snpinfo[9].split(';')]

        # counts of reads that match reference overlapping linked 'target' SNPs
        

        # linkage probabilities, not currently used
        linkageps = [float(y.strip()) for y in snpinfo[11].split(';')]

        return TestSNP(snp_id, geno_hap1, geno_hap2, AS_target_ref, 
                       AS_target_alt, hetps, tot, count)
main()
