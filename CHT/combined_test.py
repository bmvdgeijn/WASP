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

import pdb
#global log_table
#log_table=[float("-inf")]+[log(x) for x in range(1,1000000)]

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
    
    if options.pc_file:
        pc_matrix=load_covariates(options.pc_file)
        num_pcs=options.num_pcs
    else:
        pc_matrix=[]
        num_pcs=0

    if options.out_file.endswith(".gz"):
        outfile = gzip.open(options.out_file, "wb")
    else:
        outfile = open(options.out_file, 'w')

    infiles = open_input_files(options.infile_list)
    
    if (options.bnb_disp):
        disp_file=open(options.bnb_disp)
        line=disp_file.readline()
        bnb_sigmas=[]
        while line:
            bnb_sigmas.append(np.float64(line.strip()))
            line=disp_file.readline()
        disp_file.close()
    else:
        bnb_sigmas=[0.001]*len(infiles)

    if (options.as_disp):
        disp_file=open(options.as_disp)
        line=disp_file.readline()
        as_sigmas=[]
        while line:
            as_sigmas.append(np.float64(line.strip()))
            line=disp_file.readline()
        disp_file.close()
    else:
        as_sigmas=[0.001]*len(infiles)

    
    # add first row of each input file to snpinfo list
    snpinfo = []
    for f in infiles:
        snpinfo.append(f.readline().strip().split())

    row_count = 0
    finished=False
        
    while not finished: 
        try:
            test_snps=[]
            # parse test SNP and associated info from input file row
            for i in range(len(infiles)):
                test_snps.append(parse_test_snp(snpinfo[i], options))

            # how many allele-specific reads are there across all linked SNPs and
            # and individuals?
            totcounts = sum([np.sum(x.AS_target_ref) + np.sum(x.AS_target_alt) 
                             for x in test_snps])

            if totcounts < options.min_counts:
                
                if options.verbose:
                    sys.stderr.write("-----\nskipping SNP %s because "
                                     "total AS counts %d <= %d\n" % 
                                     (test_snps[0].name, totcounts, options.min_counts))

                # skip, not enough allele-specific counts
                for i in range(len(infiles)):
                    line=infiles[i].readline().strip()
                    if line:
                        snpinfo[i] = line.split()
                    else:
                        # out of lines from at least one file, assume we are finished
                        finished = True
                continue

            if options.verbose:
                sys.stderr.write("-----\ntesting SNP %s\n" % test_snps[0].name)

            
            row_count+=1
            old_genos=[test_snps[y].geno_hap1 + test_snps[y].geno_hap2 for y in range(len(test_snps))]
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

            starting_gene=[np.float64(x) for x in [0.1,0.001]] #np.float64(-4),np.float(2),np.float(10)]
            maxlike=10000000000
            for start in starting_gene:
                
                starts=[np.float64(0.5),np.float64(start)]
                # regress against the covariates and get residuals
                #fit_cov(test_snps,cov_table)

                # maximize likelihood with alpha = beta (no difference between genotypes)
                new_par = fmin(ll_one,starts, args=(test_snps, True, #options.is_bnb_only,
                                                options.is_as_only,
                                                bnb_sigmas,
                                                as_sigmas, 
                                                options.read_error_rate,
                                                [],
                                                pc_matrix),
                                disp=options.verbose,maxiter=50000,maxfun=50000)
                new_loglike = ll_one(new_par, test_snps, options.is_bnb_only,
                                 options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate,
                                 [], pc_matrix)
                if new_loglike<maxlike:
                    starting_par = new_par
                
            pc_coefs=[]
            for pc in range(num_pcs):
                new_coef=fmin(ll_pc,[np.float64(0)], args=(starting_par,test_snps, True, #options.is_bnb_only,
                                                options.is_as_only,
                                                bnb_sigmas,#options.bnb_sigma,
                                                as_sigmas, 
                                                options.read_error_rate,
                                                pc_coefs,
                                                pc_matrix),
                            disp=options.verbose,maxiter=500000,maxfun=500000)
                pc_coefs = np.concatenate([pc_coefs,new_coef])

            best1par = fmin(ll_one,starting_par, args=(test_snps, options.is_bnb_only,
                                                 options.is_as_only,
                                                 bnb_sigmas,#options.bnb_sigma,
                                                 as_sigmas, 
                                                 options.read_error_rate,
                                                 pc_coefs,
                                                 pc_matrix),
                            disp=options.verbose,maxiter=50000,maxfun=50000,ftol=1e-6, xtol=1e-6)

            if options.verbose:
                sys.stderr.write("null model optimization took %.3fs\n" % (time.time()-t1))
                
            loglike1par = ll_one(best1par, test_snps, options.is_bnb_only,
                                 options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate,
                                 pc_coefs, pc_matrix)

            #start=[best1par[0],best1par[0],best1par[1]] 
            start=[best1par[0],best1par[0],best1par[1]]
            t1=time.time()
            # maximize likelihood with alpha and beta as separate parameters
            best2par = fmin(ll_two, start, args=(test_snps,
                                                 options.is_bnb_only,
                                                 options.is_as_only,
                                                 bnb_sigmas,
                                                 as_sigmas,
                                                 options.read_error_rate,
                                                 pc_coefs,
                                                 pc_matrix),
                            disp=options.verbose,maxiter=50000,maxfun=50000,ftol=1e-6,xtol=1e-6)
            
            if options.verbose:
                sys.stderr.write("alternative model optimization took %.3fs\n" % (time.time()-t1))
            
            loglike2par = ll_two(best2par, test_snps, options.is_bnb_only,
                                 options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate,
                                 pc_coefs, pc_matrix)

            # compute likelihood ratio test statistic:
            chisq = 2*(loglike1par-loglike2par)
            

            loglike1par_up = ll_one([best1par[0],best2par[2]], test_snps, options.is_bnb_only,
                                    options.is_as_only, bnb_sigmas,
                                    as_sigmas, options.read_error_rate,
                                    pc_coefs, pc_matrix)

            if True: #2*abs(loglike2par-loglike1par) > 10:
                sys.stderr.write("%f %f %f %f %f %f %f %f %f\n" % (best1par[0],best2par[0],best2par[1],best1par[1],best2par[2],loglike1par,loglike1par_up,loglike2par,(loglike1par-loglike2par)*2))
                sys.stderr.write(", ".join([str(test_snps[i].counts) for i in range(len(test_snps))])+"\n")
                sys.stderr.write(", ".join([str(test_snps[i].totals) for i in range(len(test_snps))])+"\n")
                sys.stderr.write(", ".join([str(test_snps[i].geno_hap1 + test_snps[i].geno_hap2) for i in range(len(test_snps))])+"\n")
                sys.stderr.write(", ".join([str(old_genos[i]) for i in range(len(test_snps))])+"\n")
                

            all_counts=sum([test_snps[i].counts for i in range(len(test_snps))])
            # write result to output file
            outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], 
                                     str(chisq), str(best2par[0]),
                                     str(best2par[1]), str(best2par[2]), 
                                     str(totcounts),str(all_counts)]) + '\n')
            outfile.flush()

        except Exception as e:
            # an error occured, write to output file, but put 0s for all params and
            sys.stderr.write("An error occurred, writing line with 0s for SNP:\n%s\n" % str(e))
            
            outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], 
                                    "0", "0", "0", "0", "0"]) + '\n')
            #continue
            raise

        # read next set of lines from input file
        for i in range(len(infiles)):
            line = infiles[i].readline().strip()
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
    
    parser.add_argument("--pc-file", action='store', 
                        dest='pc_file', 
                        help="file containing PC covariates to include in the model"
                        ,default=None)

    parser.add_argument("-b", action='store', dest='bnb_disp', 
                        help="file containing depth (Beta Negative Binomial) dispersion parameters", default=None)


    parser.add_argument("-o", action='store', 
                        dest='as_disp', 
                        help="file containing allele-specific (Beta Binomial) dispersion "
                        "parameters", default=None)

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

    parser.add_argument("--num-pcs", action='store', dest='num_pcs', 
                        type=int, default=0,
                        help="designates the number of PCs to use as covariates")
    
    parser.add_argument("-v", action='store_true', dest='verbose', 
                        default=False, help="print extra information")


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

def betaln_asym(a,b):
    if b > a:
        a,b = b,a

    if a < 1e6:
        return betaln(a,b)

    l=gammaln(b)

    l -= b*math.log(a)
    l += b*(1-b)/(2*a)
    l += b*(1-b)*(1-2*b)/(12*a*a)
    l += -((b*(1-b))**2)/(12*a**3)

    return l

def BNB_loglike(k,mean,sigma,n):
    #Put variables in beta-NB form (n,a,b)
    #sys.stderr.write(str(sigma)+"\n")
    try:
        mean=max(mean,0.00001)
        logps = [math.log(n) - math.log(n + mean),
                 math.log(mean) - math.log(n + mean)]
    except:
        raise
        n_val=n
        pdb.set_trace()

    p=np.float64(n/(n+mean))

    if sigma < 0.00001: #> 18: #20:
        loglike=-betaln(n,k+1)-math.log(n+k)+n*logps[0]+k*logps[1]
        return loglike

    sigma=(1/sigma)**2 #+sigma*n
    sigma=sigma #+math.sqrt(sigma)/(p*(1-p))**2

    a = p*sigma+1
    b = (1-p)*sigma

    #Rising Pochhammer = gamma(k+n)/gamma(n)
    if k>0:
        loglike=-betaln_asym(n,k)-math.log(k)
    else:
        loglike=0

    #Add log(beta(a+n,b+k))
    loglike += betaln_asym(a+n,b+k)
    
    #Subtract log(beta(a,b))
    loglike -= betaln_asym(a,b)

    return loglike

#def ll_pc(x, test_snps, is_bnb_only, is_as_only, bnb_sigma, as_sigma,error, pcs):
#def fit_cov(test_snps, cov_table):
#    if len(cov_table)==0:
#        return
#    counts=[snp.counts for snp in test_snps]
#    lambdas=[snp.totals for snp in test_snps]
#    sys.stderr.write(str(counts)+"\n")
#    sys.stderr.write(str(lambdas)+"\n")
#    ys=[(counts[i]-lambdas[i])/lambdas[i] for i in range(len(counts))]
#    sys.stderr.write(str(ys)+"\n")
#    sys=np.log(ys)
#    fit=np.linalg.lstsq(ys,cov_table)
#    fitted_values=(ys-residuals+lambdas)*lambdas
#    
#    for i in range(len(fitted_values)):
#        test_snps[i].totals=fitted_values[i]

def ll_one(x, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix):
    alpha = x[0]
    beta = x[0]
    r = x[1]
    return loglikelihood(alpha,beta,r, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix)

def ll_pc(x, params, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,other_pc_coefs,pc_matrix):
    alpha = params[0]
    beta = params[0]
    r = params[1]
    pc_coefs=np.concatenate([other_pc_coefs,x])
    return loglikelihood(alpha,beta,r, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix)

def ll_two(x, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix):
    alpha = x[0]
    beta = x[1]
    r = x[2]
    #if len(x)>3:
    #    pc_fits=x[3:]
    #else:
    #    pc_fits=[]
    return loglikelihood(alpha,beta,r, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix)

def calc_pc_factor(pc_fits,pcs,i):
    if len(pc_fits)>0:
        return 1+sum(pc_fits*pcs[i,:len(pc_fits)])
    else:
        return 1

def loglikelihood(alpha,beta,r, test_snps, is_bnb_only, is_as_only, bnb_sigmas, as_sigmas,error,pc_coefs,pc_matrix): 
    loglike = 0

    #if input values are outside of reasonable range return a very high -loglike
    if alpha <= 0 or beta <= 0 or r<=0 or r>1:
        return 10000000
    #r=math.exp(r)/(1+math.exp(r)) #keep r between 0 and 1
    #r=(r+0.01)/(1.01)
    #r=math.exp(r)
    ratio = (alpha / (alpha + beta))

    for i in range(len(test_snps)):
        #if i in (17,36,49,51):
        #    continue
        if(test_snps[i].is_homo_ref()):
            m = 2*alpha*test_snps[i].totals *calc_pc_factor(pc_coefs,pc_matrix,i)
        elif(test_snps[i].is_homo_alt()):
            m = 2*beta*test_snps[i].totals *calc_pc_factor(pc_coefs,pc_matrix,i)
        else:
            m = (alpha+beta)*test_snps[i].totals *calc_pc_factor(pc_coefs,pc_matrix,i)
        if m<0:
            m=0.000001
        if not is_bnb_only:
            for j in range(len(test_snps[i].AS_target_ref)):            
                if test_snps[i].hetps[j]>.9:
                    hetp = min(.99,test_snps[i].hetps[j])
                    logps = [math.log(alpha) - math.log(alpha+beta),
                             math.log(beta) - math.log(alpha+beta)]
                    loglike += AS_betabinom_loglike(logps, as_sigmas[i], 
                                                    test_snps[i].AS_target_ref[j],
                                                    test_snps[i].AS_target_alt[j],
                                                    hetp, error)
        if not is_as_only:
            #sys.stderr.write(str(bnb_sigmas)+"\n")
            #sys.stderr.write(str(bnb_sigmas[i])+"\n")
            l = BNB_loglike(test_snps[i].counts, m, r, bnb_sigmas[i])
            loglike += l
    return -loglike


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

    #if snpinfo[9].strip() == "NA":
        # SNP is homozygous, so there is no AS info
    #    return TestSNP(snp_id, geno_hap1, geno_hap2, [], [], [], tot, count)    
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
        hetps = [np.float64(y.strip()) for y in snpinfo[10].split(';')]

        # linkage probabilities, not currently used
        linkageps = [np.float64(y.strip()) for y in snpinfo[11].split(';')]

        if options.shuffle:
            # permute allele-specific read counts by flipping them randomly at
            # each SNP
            for y in range(len(AS_target_ref)):
                if randint(0,1) == 1:
                    temp=AS_target_ref[y]
                    AS_target_ref[y] = AS_target_alt[y]
                    AS_target_alt[y] = temp

        return TestSNP(snp_id, geno_hap1, geno_hap2, AS_target_ref, 
                       AS_target_alt, hetps, tot, count)

def load_covariates(cov_file):
    infile=open(cov_file)
    cov_table=[]
    while True:
        line=infile.readline()
        if line:
            cov_table.append([np.float64(x) for x in line.strip().split()])
        else:
            break
    
    return np.array(cov_table,dtype=np.float64)
        
main()
