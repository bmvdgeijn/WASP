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

#import pdb

class TestSNP:
    def __init__(self, name, geno_hap1, geno_hap2, AS_target_ref, AS_target_alt, 
                 hetps, totals, counts):
        self.name = name
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
        sys.stderr.write(filename+"\n")
        if not filename or not os.path.exists(filename) or not os.path.isfile(filename):
            sys.stderr.write("input file '%s' does not exist or is not a regular file\n" 
                             % in_file)
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


def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("infile_list", action='store', default=None)
    parser.add_argument("outfile", action='store', default=None)
    return parser.parse_args()

def main():
    options=parse_options()
    infiles = open_input_files(options.infile_list)
    
    # add first row of each input file to snpinfo list
    snpinfo = []
    for f in infiles:
        f.readline()
        snpinfo.append(f.readline().strip().split())

    row_count = 0
    finished=False
    count_matrix=[]
    expected_matrix=[]
        
    while not finished: 
        count_line=[]
        expected_line=[]
        # parse test SNP and associated info from input file row
        for i in range(len(infiles)):
            new_snp=parse_test_snp(snpinfo[i], options)
            count_line.append(new_snp.counts)
            expected_line.append(new_snp.totals)
            try:
                line=infiles[i].readline().strip()
            except:
                finished=True
                break
            if line:
                snpinfo[i] = line.split()
            else:
                # out of lines from at least one file, assume we are finished
                finished = True
                continue
        if sum(count_line)>3000:
            count_matrix.append(count_line)
            expected_matrix.append(expected_line)
    
    count_matrix=np.array(count_matrix,dtype=int)
    expected_matrix=np.array(expected_matrix,dtype=np.float64)

    gw_fit_starts=[0.005,.01,.015,.02,.07,.09]
    old_ll=10000000000000000000000000000000
    best_start=-1
    if False:
        for i in range(len(gw_fit_starts)):
            gw_fits=[gw_fit_starts[i]]*count_matrix.shape[1]
            gene_fits=[1]*count_matrix.shape[0]
            gene_fits,fit_ll=get_gene_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits)
            gw_gits,fit_ll=get_gw_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits)
            if fit_ll<old_ll:
                old_ll=fit_ll
                best_start=i
    
    gw_fits=[gw_fit_starts[3]]*count_matrix.shape[1]
    gene_fits=[100]*count_matrix.shape[0]
 
    iteration=1
    while True:
        gene_fits,fit_ll=get_gene_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits,iteration)

        gw_fits,fit_ll=get_gw_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits)
        sys.stderr.write("%f\n"%fit_ll)
        outfile=open(options.outfile,"w")
        for i in gw_fits:
            outfile.write("%f\n" % i[0])
        outfile.close()
        iteration+=1

        if old_ll-fit_ll<.001:
            break
        old_ll=fit_ll

def get_gene_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits,iteration):
    fit_ll=0
    for gene_indx in range(count_matrix.shape[0]):
        ### fit new gene dispersion parameters
        if iteration<=2:
            starts=[gene_fits[gene_indx],1,1000]
        else:
            starts=[gene_fits[gene_indx]]
        best_fit=-1
        best_like=100000000
        for par_start in starts:
            new_fit= fmin(gene_like,par_start, 
                                   args=(count_matrix[gene_indx,:],
                                         expected_matrix[gene_indx,:],
                                         gw_fits),
                                   disp=False,maxiter=5000,maxfun=5000,ftol=0.0001,xtol=0.0001)
            new_like=gene_like(new_fit,count_matrix[gene_indx,:],expected_matrix[gene_indx,:],gw_fits)
            if new_like<best_like:
                best_like=new_like
                gene_fits[gene_indx]=new_fit
        if gene_indx % 20 == 1:
            like=gene_like(gene_fits[gene_indx],count_matrix[gene_indx,:],expected_matrix[gene_indx,:],gw_fits)
            sys.stderr.write("%f %f\n"%(gene_fits[gene_indx],like))
        fit_ll+=gene_like(gene_fits[gene_indx],count_matrix[gene_indx,:],expected_matrix[gene_indx,:],gw_fits)
    sys.stderr.write("\n")
    return gene_fits, fit_ll

def get_gw_overdisp(count_matrix,expected_matrix,gw_fits,gene_fits):
    fit_ll=0
    for indx in range(count_matrix.shape[1]):
        gw_fits[indx]= fmin(gw_like,gw_fits[indx], 
                              args=(count_matrix[:,indx],
                                    expected_matrix[:,indx],
                                    gene_fits),
                              disp=False,maxiter=50000,maxfun=50000,ftol=0.0001,xtol=0.0001)
        like=gw_like(gw_fits[indx],count_matrix[:,indx],expected_matrix[:,indx],gene_fits)
        fit_ll+=like
    return gw_fits,fit_ll


def gene_like(gene_fit,counts,expecteds,gw_fits):
    loglike=0
    if gene_fit<=0:
        return 100000
    for i in range(len(counts)):
        loglike+=BNB_loglike(counts[i],expecteds[i],gene_fit,gw_fits[i])
    return -loglike

def gw_like(gw_fit,counts,expecteds,gene_fits):
    loglike=0
    if gw_fit<=0.0001:
        return 100000
    for i in range(len(counts)):
        loglike+=BNB_loglike(counts[i],expecteds[i],gene_fits[i],gw_fit)
    return -loglike


def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))

def BNB_loglike(k,mean,n,sigma):
    #Put variables in beta-NB form (n,a,b)
    logps = [math.log(n) - math.log(n + mean),
             math.log(mean) - math.log(n + mean)]
    
    a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
    b = math.exp(logps[1] + math.log(1/sigma**2 - 1))

    loglike = 0
    
    #Rising Pochhammer = gamma(k+n)/gamma(n)
    if k>0:
        loglike=-betaln(n,k)-math.log(k)
    else:
        loglike=0
    
    #Add log(beta(a+n,b+k))
    loglike += betaln(a+n,b+k)
    
    #Subtract log(beta(a,b))
    loglike -= betaln(a,b)

    return loglike

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
        AS_target_ref = [int(y) for y in snpinfo[12].split(';')]

        # counts of reads that match alternate allele
        AS_target_alt = [int(y) for y in snpinfo[13].split(';')]

        # heterozygote probabilities
        hetps = [float(y.strip()) for y in snpinfo[10].split(';')]

        # linkage probabilities, not currently used
        linkageps = [float(y.strip()) for y in snpinfo[11].split(';')]

        return TestSNP(snp_id, geno_hap1, geno_hap2, AS_target_ref, 
                       AS_target_alt, hetps, tot, count)

main()
