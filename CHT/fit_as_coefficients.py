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
import gzip
import argparse

from scipy.optimize import *
from scipy.special import betaln
import scipy.stats

import numpy as np


def parse_options():
    parser = argparse.ArgumentParser(description="This script estimates the "
                                     "overdispersion parameter for the allele-specific "
                                     "(Beta-Binomial) half of the combined halotype test."
                                     " A single overdispersion parameter is estimated for "
                                     "each individual (across sites), under the assumption "
                                     "that all sites come from the null hypothesis (no "
                                     "genetic association)");
    
    
    parser.add_argument("-e", action='store', dest='read_error_rate',
                        help="sequence read error rate",
                        type=float, default=0.005)
    
    parser.add_argument("infile_list", 
                        help="Path to file containing list of CHT input "
                        "files (one for each individual)",
                        action='store', default=None)
    
    parser.add_argument("out_file", 
                        help="File to write overdispersion parameter estimates to",
                        action='store', default=None)
    
    return parser.parse_args()




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
        dispersion = fmin(likelihood, 0.001, 
                          args=(AS_ref, AS_alt, hetps, options.read_error_rate))
        outfile.write(str(dispersion[0])+"\n")
        outfile.flush()

        
def likelihood(dispersion,AS_ref,AS_alt,hetps,error):
    cur_like=0
    for i in range(len(AS_ref)):
        cur_like = cur_like + AS_betabinom_loglike([math.log(0.5), math.log(0.5)], 
                                                   dispersion, AS_ref[i], 
                                                   AS_alt[i], hetps[i], error)
    return -cur_like




def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))



#Given parameters, returns log likelihood.  Note that some parts have been cancelled out
def AS_betabinom_loglike(logps, sigma, AS1, AS2, hetp, error):
    if sigma >= 1.0 or sigma <= 0.0:
        return -99999999999.0
    else:
        a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
        b = math.exp(logps[1] + math.log(1/sigma**2 - 1))
    
    part1 = 0
    part1 += betaln(AS1 + a, AS2 + b)
    part1 -= betaln(a, b)
    
    if hetp == 1:
        return part1        

    e1 = math.log(error) * AS1 + math.log(1 - error) * AS2
    e2 = math.log(error) * AS2 + math.log(1 - error) * AS1
    if hetp == 0:
        return addlogs(e1, e2)
    
    return addlogs(math.log(hetp)+part1, math.log(1-hetp) + addlogs(e1,e2))

main()
