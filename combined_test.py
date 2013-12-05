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


import math
from scipy.optimize import *
from scipy import cast
from scipy.special import gammaln
from scipy.special import betaln
import sys
import numpy as np
from random import shuffle
from random import randint
import scipy.stats
import time
import gzip, argparser

class Test_SNP:
    def __init__(self, geno_hap1, geno_hap2, AS_target_ref, AS_target_alt, hetps, totals, counts):
        self.geno_hap1 = geno_hap1
        self.geno_hap2 = geno_hap2
        self.AS_target_ref = AS_target_ref
        self.AS_target_alt = AS_target_alt
        self.hetps = hetps
        self.totals = totals
        self.counts = counts

    def is_het(self):
        return self.geno_hap1 != self.geno_hap2

    def is_homo_ref(self):
        return self.geno_hap1 == 0 and self.geno_hap2 == 0

    def is_homo_alt(self):
        return self.geno_hap1 == 1 and self.geno_hap2 == 1
    
#helper function for addition of logs
def addlogs(loga,logb):
    return max(loga,logb)+math.log(1+math.exp(-abs(loga-logb)))

#Given parameters, returns log likelihood.  Note that some parts have been cancelled out
def AS_betabinom_loglike(logps,sigma,AS1,AS2,hetp,error):

    a=math.exp(logps[0]+math.log(1/sigma**2 - 1))
    b=math.exp(logps[1]+math.log(1/sigma**2 - 1))
    
    part1=0
    part1+=betaln(AS1+a,AS2+b)    
    part1-=betaln(a,b)
    
    if hetp==1:
        return part1        
    e1=math.log(error)*AS1+math.log(1-error)*AS2
    e2=math.log(error)*AS2+math.log(1-error)*AS1
    if hetp==0:
        return addlogs(e1,e2)
    return addlogs(math.log(hetp)+part1,math.log(1-hetp)+addlogs(e1,e2))


def BNB_loglike(k,mean,n,sigma):
    #Put variables in beta-NB form (n,a,b)

    logps=[math.log(n)-math.log(n+mean),math.log(mean)-math.log(n+mean)]
    a=math.exp(logps[0]+math.log(1/sigma**2 - 1))
    b=math.exp(logps[1]+math.log(1/sigma**2 - 1))

    loglike = 0
    
    #Rising Pochhammer = gamma(k+n)/gamma(n)
    for j in range(k):
        loglike+=math.log(j+n)

    #Add log(beta(a+n,b+k))
    loglike+=betaln(a+n,b+k)
    
    #Subtract log(beta(a,b))
    loglike-=betaln(a,b)

    return loglike
          
def loglikelihood(x,test_snps,is_nb_only,is_as_only,bnb_sigma,as_sigma,error): 
    if len(x)==3:
        alpha=x[0]
        beta=x[1]
        r=x[2]
    else:
        alpha=x[0]
        beta=x[0]
        r=x[1]
    loglike=0
    ratio=(alpha/(alpha+beta))

    #if input values are outside of reasonable range return a very high -loglike
    if alpha<=0 or beta<=0 or r<=0:
        return 10000000000000000000000

    for i in range(len(test_snps)):
        if(test_snps[i].is_homo_ref()):
            m=2*alpha*test_snps[i].totals
        elif(test_snps[i].is_homo_alt()):
            m=2*beta*test_snps[i].totals
        else:
            m=(alpha+beta)*test_snps[i].totals
        if not is_bnb_only:
            for j in range(len(test_snps[i].AS_target_ref)):            
                if test_snps[i].hetps[j]>.9:
                    hetp = test_snps[i].hetps[j]
                    logps=[math.log(alpha)-math.log(alpha+beta),math.log(beta)-math.log(alpha+beta)]
                    loglike+=AS_betabinom_loglike(logps,as_sigma,test_snps[i].AS_target_ref[j],test_snps[i].AS_target_alt[j],hetp,error)
        if not is_as_only:
            l=BNB_loglike(test_snps[i].counts,m,r,bnb_sigma)
            loglike+=l
    return -loglike

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-a", action='store_true', dest='is_as_only', default=False)
    parser.add_argument("-d", action='store_true', dest='is_bnb_only', default=False)
    
    parser.add_argument("-o", action='store', dest='as_sigma', type=float, default=0.00001)
    parser.add_argument("-b", action='store', dest='bnb_sigma', type=float, default=0.00001)

    parser.add_argument("-s", action='store_true', dest='shuffle', default=False)    
    parser.add_argument("-e", action='store', dest='read_error_rate', type=float, default=0.005)
    parser.add_argument("-m", action='store', dest='min_counts', type=int, default=0)

    parser.add_argument("infile_list", action='store')
    parser.add_argument("out_file", action='store')
    options=parser.parse_args()

    outfile=open(options.out_file,'w')

    infiles=[]
    snpinfo=[]

    infile_list=open(options.infile_list)
    for infile in infile_list:
        infiles.append(open(infile))
        snpinfo.append(infiles[-1].readline().strip().split())    
    infile_list.close()

    row_count=0
    finished=False
    while not finished:
        try:
            test_snps=[]
            for i in range(len(infiles)):
                tot=int(snpinfo[i][16])
                geno_hap1=int(snpinfo[i][6].strip().split("|")[0])
                geno_hap2=int(snpinfo[i][6].strip().split("|")[1])
                count=int(snpinfo[i][15])
               
                if snpinfo[i][9].strip()!="NA":
                    snplocs=[int(y.strip()) for y in snpinfo[i][9].split(';')]
                    AS_target_ref = [int(y) for y in snpinfo[i][12].split(';')]
                    AS_target_alt = [int(y) for y in snpinfo[i][13].split(';')]
                    totcounts+=np.sum(AS_target_ref)+np.sum(AS_target_alt)
                    
                    hetps=[float(y.strip()) for y in snpinfo[i][10].split(';')]
                    linkageps=[float(y.strip()) for y in snpinfo[i][11].split(';')]
                    
                    if options.shuffle:
                        for y in range(len(AS_target_ref)):
                            if randint(0,1) == 1:
                                temp=AS_target_ref[y]
                                AS_target_ref[y]=AS_target_alt[y]
                                AS_target_alt[y]=temp
                        
                    test_snps.append(Test_SNP(geno_hap1, geno_hap2, AS_target_ref, AS_target_alt, hetps, tot, count))
                else:
                    test_snps.append(Test_SNP(geno_hap1, geno_hap2, [], [], [], tot, count))

        if totcounts>options.min_counts:
            row_count+=1
            if shuff==1:
                perm=range(len(test_snps))
                shuffle(perm)
                geno1temp=[test_snps[perm[y]].geno_hap1 for y in range(len(perm))]
                geno2temp=[test_snps[perm[y]].geno_hap2 for y in range(len(perm))]
                for i in range(len(test_snps)):
                    test_snps[i].geno_hap1 = geno1temp[i]
                    test_snps[i].geno_hap2 = geno2temp[i]

            t1=time.time()
            best1par=fmin(loglikelihood,(20,10),args=(test_snps,options.is_bnb_only,options.is_as_only,options.bnb_sigma,options.as_sigma,error))
            print str(time.time()-t1)
            loglike1par= loglikelihood(best1par,test_snps,options.is_bnb_only,options.is_as_only,options.bnb_sigma,options.as_sigma,error)
            t1=time.time()
            best2par=fmin(loglikelihood,(best1par[0],best1par[0],best1par[1]),args=(test_snps,options.is_bnb_only,options.is_as_only,options.bnb_sigma,options.as_sigma,error))
            print str(time.time()-t1)
            sys.stdout.flush()
            loglike2par= loglikelihood(best2par,test_snps,options.is_bnb_only,options.is_as_only,options.bnb_sigma,options.as_sigma,error)
            outfile.write(snpinfo[0][0]+"\t"+snpinfo[0][1]+"\t"+str(2*(loglike1par-loglike2par))+'\t'+str(best2par[0])+'\t'+str(best2par[1])+"\t"+str(best2par[2])+"\t"+str(totcounts)+'\n')
        outfile.flush()
    except:
        outfile.write(snpinfo[0][0]+"\t"+snpinfo[0][1]+"\t"+str(0)+'\t'+str(0)+'\t'+str(0)+"\t"+str(0)+"\t"+str(0)+'\n')
    for i in range(len(infiles)):
        line=infiles[i].readline().strip()
        if line:
            snpinfo[i]=line.split()
        else:
            finished = True
