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


class Test_SNP:
    def __init__(self, geno_hap1, geno_hap2, AS_target_ref, AS_target_alt, 
                 hetps, totals, counts):
        self.geno_hap1 = geno_hap1
        self.geno_hap2 = geno_hap2
        self.AS_target_ref = AS_target_ref
        self.AS_target_alt = AS_target_alt
        self.hetps = hetps
        self.totals = totals
        self.counts = counts

        
    def is_het(self):
        """returns True if the test SNP is heterozygous"""        
        return self.geno_hap1 != self.geno_hap2

    def is_homo_ref(self):
        """Returns True if test SNP is homozygous for reference allele"""
        return self.geno_hap1 == 0 and self.geno_hap2 == 0

    def is_homo_alt(self):
        """Returns True if test SNP is homozygous for non-reference allele"""
        return self.geno_hap1 == 1 and self.geno_hap2 == 1



    
def main():
    options = parse_options()

    outfile = open(options.out_file, 'w')

    infiles = []
    snpinfo = []

    if not os.path.exists(options.infile_list) or \
      not os.path.isfile(options.infile_list):
        sys.stderr.write("input file %s does not exist or is not a regular file\n" %
                         options.infile_list)
        exit(2)
    
    # read file that contains list of input files
    infile_list = open(options.infile_list)
    for infile in infile_list:
        # open each input file and read first line
        infile = infile.rstrip()
        if not infile or not os.path.exists(infile) or not os.path.isfile(infile):
            sys.stderr.write("input file '%s' does not exist or is not a regular file\n" 
                             % infile)
            exit(2)
        infiles.append(open(infile))
        infiles[-1].readline() #skip the header
        snpinfo.append(infiles[-1].readline().strip().split())
    infile_list.close()

    if len(infiles) == 0:
        sys.stderr.write("no input files specified in file '%s'\n" % options.infile_list)
        exit(2)
    
    row_count = 0
    finished=False
    
    while not finished: 
        try:
            test_snps=[]
            # parse SNP info from each individual file
            for i in range(len(infiles)):
                test_snps.append(parse_test_snp(snpinfo[i], options))

            # how many allele-specific reads are there across all linked SNPs and
            # and individuals?
            totcounts = sum([np.sum(x.AS_target_ref) + np.sum(x.AS_target_alt) 
                             for x in test_snps])

            if totcounts <= options.min_counts:
                # skip, not enough allele-specific counts
                for i in range(len(infiles)):
                    line=infiles[i].readline().strip()
                    if line:
                        snpinfo[i] = line.split()
                    else:
                        # out of lines from at least one file, assume we are finished
                        finished = True
                continue

            row_count+=1
            if options.shuffle:
                # permute genotypes
                perm = range(len(test_snps))
                shuffle(perm)
                geno1temp = [test_snps[y].geno_hap1 for y in perm]
                geno2temp = [test_snps[y].geno_hap2 for y in perm]
                for i in range(len(test_snps)):
                    test_snps[i].geno_hap1 = geno1temp[i]
                    test_snps[i].geno_hap2 = geno2temp[i]
            t1=time.time()
            # maximize likelihood with alpha = beta (no difference between genotypes)
            best1par = fmin(loglikelihood,(20,10), args=(test_snps, options.is_bnb_only,
                                                         options.is_as_only,
                                                         options.bnb_sigma,
                                                         options.as_sigma, options.read_error_rate))
            
            print str(time.time()-t1)
            loglike1par = loglikelihood(best1par, test_snps, options.is_bnb_only,
                                        options.is_as_only, options.bnb_sigma,
                                        options.as_sigma, options.read_error_rate)
            t1=time.time()
            # maximize likelihood with alpha and beta as separate parameters
            best2par = fmin(loglikelihood, (best1par[0], best1par[0], best1par[1]),
                            args=(test_snps, options.is_bnb_only,options.is_as_only,
                                  options.bnb_sigma,options.as_sigma,options.read_error_rate))
            
            print str(time.time()-t1)
            
            sys.stdout.flush()
            
            loglike2par = loglikelihood(best2par, test_snps, options.is_bnb_only,
                                        options.is_as_only, options.bnb_sigma,
                                        options.as_sigma, options.read_error_rate)

            # compute likelihood ratio test statistic:
            chisq = 2*(loglike1par-loglike2par)

            # write result to output file
            outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], 
                                     str(chisq), str(best2par[0]),
                                     str(best2par[1]), str(best2par[2]), 
                                     str(totcounts)]) +'\n')
            outfile.flush()
        
        except:
            # an error occured, write to output file, but put 0s for all params and 
            outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], 
                                    "0", "0", "0", "0", "0"]) + '\n')

        # read next set of lines from input file
        for i in range(len(infiles)):
            line=infiles[i].readline().strip()
            if line:
                snpinfo[i] = line.split()
            else:
                # out of lines from at least one file, assume we are finished
                finished = True



def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("-a", action='store_true', 
                        dest='is_as_only', default=False,
                        help="only perform the allele-specific part (Beta Binomial) "
                        "part of the test")
    
    parser.add_argument("-d", action='store_true', 
                        dest='is_bnb_only', default=False,
                        help="only perform the association (Beta Negative Binomial) part "
                        "of the test")
    
    parser.add_argument("-o", action='store', 
                        dest='as_sigma', type=float, 
                        help="value for allele-specific (Beta Binomial) dispersion "
                        "parameter", default=0.00001)
    
    parser.add_argument("-b", action='store', dest='bnb_sigma', 
                        help="value for global Beta Negative Binomial dispersion parameter",
                        type=float, default=0.00001)

    parser.add_argument("-s", action='store_true', 
                        dest='shuffle', default=False,
                        help="permute genotypes")
    
    parser.add_argument("-e", action='store', dest='read_error_rate',
                        help="estimate of error rate, used to update "
                        "heterozygous genotype probabilities "
                        "(currently this option disabled / not used)",
                        type=float, default=0.005)
    
    parser.add_argument("-m", action='store', dest='min_counts', 
                        type=int, default=0,
                        help="only perform test when total number of allele-specific "
                        "read counts across individuals > MIN_COUNTS")

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


def BNB_loglike(k,mean,n,sigma):
    #Put variables in beta-NB form (n,a,b)
    logps = [math.log(n) - math.log(n + mean),
             math.log(mean) - math.log(n + mean)]
    
    a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
    b = math.exp(logps[1] + math.log(1/sigma**2 - 1))

    loglike = 0
    
    #Rising Pochhammer = gamma(k+n)/gamma(n)
    for j in range(k):
        loglike += math.log(j+n)

    #Add log(beta(a+n,b+k))
    loglike += betaln(a+n,b+k)
    
    #Subtract log(beta(a,b))
    loglike -= betaln(a,b)

    return loglike


def loglikelihood(x, test_snps, is_bnb_only, is_as_only, bnb_sigma, as_sigma,error): 
    if len(x) == 3:
        # model with separate alpha and beta params
        alpha = x[0]
        beta = x[1]
         # per-region dispersion param, r:
        r = x[2]
    else:
        # model with single param for alpha = beta
        alpha = x[0]
        beta = x[0]
        r = x[1]
    loglike = 0
    ratio = (alpha / (alpha + beta))

    #if input values are outside of reasonable range return a very high -loglike
    if alpha <= 0 or beta <= 0 or r <= 0:
        return 10000000000000000000000

    for i in range(len(test_snps)):
        if(test_snps[i].is_homo_ref()):
            m = 2*alpha*test_snps[i].totals
        elif(test_snps[i].is_homo_alt()):
            m = 2*beta*test_snps[i].totals
        else:
            m = (alpha+beta)*test_snps[i].totals
        if not is_bnb_only:
            for j in range(len(test_snps[i].AS_target_ref)):            
                if test_snps[i].hetps[j]>.9:
                    hetp = test_snps[i].hetps[j]
                    logps = [math.log(alpha) - math.log(alpha+beta),
                             math.log(beta) - math.log(alpha+beta)]
                    loglike += AS_betabinom_loglike(logps, as_sigma, 
                                                    test_snps[i].AS_target_ref[j],
                                                    test_snps[i].AS_target_alt[j],
                                                    hetp, error)
        if not is_as_only:
            l = BNB_loglike(test_snps[i].counts, m, r, bnb_sigma)
            loglike += l
    return -loglike





def parse_test_snp(snpinfo, options):
    tot = int(snpinfo[16])/100000
    geno_hap1 = int(snpinfo[6].strip().split("|")[0])
    geno_hap2 = int(snpinfo[6].strip().split("|")[1])
    count = int(snpinfo[15])

    if snpinfo[9].strip() == "NA":
        # SNP is homozygous, so there is no AS info
        return Test_SNP(geno_hap1, geno_hap2, [], [], [], tot, count)    
    else:
        # positions of target SNPs (not currently used)
        snplocs=[int(y.strip()) for y in snpinfo[9].split(';')]

        # counts of reads that match reference overlapping linked 'target' SNPs
        AS_target_ref = [int(y) for y in snpinfo[12].split(';')]

        # counts of reads that match alternate allele
        AS_target_alt = [int(y) for y in snpinfo[13].split(';')]

        # heterozygote probabilities
        hetps = [float(y.strip()) for y in snpinfo[10].split(';')]

        # linkage probabilities, not currently used
        linkageps = [float(y.strip()) for y in snpinfo[11].split(';')]

        if options.shuffle:
            # permute allele-specific read counts by flipping them randomly at
            # each SNP
            for y in range(len(AS_target_ref)):
                if randint(0,1) == 1:
                    temp=AS_target_ref[y]
                    AS_target_ref[y] = AS_target_alt[y]
                    AS_target_alt[y] = temp

        return Test_SNP(geno_hap1, geno_hap2, AS_target_ref, 
                       AS_target_alt, hetps, tot, count)
        
main()
