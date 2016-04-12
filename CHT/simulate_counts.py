
from random import choice, sample, random
from numpy.random import beta, negative_binomial, binomial
import numpy as np
import sys

import argparse


def parse_options():
    parser = argparse.ArgumentParser(description="simulate counts for combined haplotype test")

    parser.add_argument("--prefix", default=None, required=True,
                        help="prefix for output files")

    dflt_num_tests = 10000
    parser.add_argument("--num_tests", default=dflt_num_tests, type=int,
                        help="number of test SNPs / regions to simulate "
                        "(default=%d)" % dflt_num_tests)

    dflt_num_inds = 10
    parser.add_argument("--num_inds", default=dflt_num_inds, type=int,
                        help="number of individuals to simulate "
                        "(default=%d)" % dflt_num_inds)

    dflt_min_hets = 2
    parser.add_argument("--min_hets", default=dflt_min_hets, type=int,
                        help="minimum number of heterozygous individuals "
                        "per test SNP (default=2)" % min_hets)

    dflt_maf = 0.2
    parser.add_argument("--maf", default=dflt_maf, type=float,
                        help="minor allele frequency of test SNP "
                        "(default = %.2f)" % dflt_maf)

    dflt_mean_counts = 200
    parser.add_argument("--mean_counts", default=dflt_mean_counts, type=int,
                        help="mean number of read counts per region "
                        "(default=%d)" % dflt_mean_counts)

    dflt_as_counts = 20
    parser.add_argument("--as_counts", default=dflt_as_counts, type=int,
                        help="expected number of allele-specific read counts "
                        "per regions (default=%d)", dflt_as_counts)

    dflt_gene_disp = 0.01
    parser.add_argument("--gene_disp", default=dflt_gene_disp, type=float,
                        help="per-gene overdispersion parameter for "
                        "beta-negative binomial", dflt_gene_disp)

    # TODO: could take array of overdispersion parameters (one for each individual)
    dflt_ind_disp = 100.0
    parser.add_argument("--ind_disp", default=dflt_ind_disp, type=float,
                        help="per individual overdispersion parameter for "
                        "beta-negative binomial (default=%.2f)", dflt_ind_disp)

    dflt_as_disp = 0.2
    parser.add_argument("--as_disp", default=dflt_as_disp, type=float,
                        help="per individual allele-specific overdispersion "
                        "parameter for beta-binomial (default=%.2f)", dflt_as_disp)

    dflt_effect_size = 0.2
    parser.add_argument("--effect_size", default=dflt_effect_size, type=float,
                        help="effect size of true positives (default=%.2f)" %
                        dflt_effect_size)

    dflt_additivity = 1.0
    parser.add_argument("--additivity", default=dflt_additivity, type=float,
                        help="additivity of alleles (default=%.2f)" %
                        dflt_additivity)

    dflt_het_error_rate = 0.01
    parser.add_argument("--het_error_rate", default=dflt_het_error_rate,
                        type=float, help="rate of incorrect heterozygous genotype calls"
                        "(default=%.2f)" % dflt_het_error_rate)
    
    dflt_read_error_rate = 0.01
    parser.add_argument("--read_error_rate", default=dflt_read_error_rate,
                        type=float, help="rate of incorrect alleles in reads"
                        "(default=%.2f)" % dflt_read_error_rate)

    dflt_true_positives = 0.05
    parser.add_argument("--true_positives", default=dflt_true_positives,
                        type=float, help="fraction of test SNPs that are "
                        "true positives (default=%.2f)" % dflt_true_positives)

    dflt_sim_hom = False
    parser.add_argument("--sim_hom", action="store_true" dest="sim_hom",
                        type=bool, help="simulate homozygous test SNPs"
                        default=False)

    return parser.parse_args()
    


def main():
    options = parse_options()
    
    file_base=sys.argv[1]
    num_tests=int(sys.argv[2]) #10000
    num_inds=int(sys.argv[3]) #10
    min_hets=int(sys.argv[4]) #2
    MAF=float(sys.argv[5]) #0.2
    mean_counts=int(sys.argv[6]) #200
    as_counts=int(sys.argv[7]) #20
    gene_disp=float(sys.argv[8]) #0.01
    ind_disp=int(sys.argv[9]) #100
    AS_disp=float(sys.argv[10]) #0.2
    effect_size=float(sys.argv[11]) #0.2
    additivity=float(sys.argv[12]) #1
    het_error_rate=float(sys.argv[13]) #0.01
    read_error_rate=float(sys.argv[14]) #0.01
    true_positives=float(sys.argv[15]) #0.05
    simulate_hom_AS=True #False
    
    out_files=[]
    for i in range(1,num_inds+1):
        out_files.append(open("/mnt/gluster/home/bmvdgeijn/WASP_testing/SIMULATED_COUNTS/WASP/%s_%i.txt"%(file_base,i),"w"))
        out_files[i-1].write("Header\n")
    ASseq_Y=open("/mnt/gluster/home/bmvdgeijn/WASP_testing/SIMULATED_COUNTS/ASseq/%s_Y.txt" % (file_base),"w")
    ASseq_Y1=open("/mnt/gluster/home/bmvdgeijn/WASP_testing/SIMULATED_COUNTS/ASseq/%s_Y1.txt" % (file_base),"w")
    ASseq_Y2=open("/mnt/gluster/home/bmvdgeijn/WASP_testing/SIMULATED_COUNTS/ASseq/%s_Y2.txt" % (file_base),"w")
    ASseq_Z=open("/mnt/gluster/home/bmvdgeijn/WASP_testing/SIMULATED_COUNTS/ASseq/%s_Z.txt" % (file_base),"w")

    
    test=1
    while test <= num_tests:
        if random()>true_positives:
            effect=0
            alt_expr=1
            AS_frac=0.5
        elif random()<0.5:
            effect=0
            alt_expr=1+effect_size
            AS_frac=1/(2+additivity*effect_size)
        else:
            effect=1
            alt_expr=1/(1+effect_size)
            AS_frac=(1+additivity*effect_size)/(2+additivity*effect_size)

        snps=[]
        counts=[]
        num_hets=0
        for ind in range(num_inds):
            # Simulate the individual's haps=[0,0]
            minors=int(random()<MAF) + int(random()<MAF)
            if minors==0:
                haps=[0,0]
            elif minors==1:
                haps=[0,1]
                num_hets+=1
            else:
                haps=[1,1]

            # Expected number of reads based on genotypes
            ind_mean_counts=mean_counts*((2-sum(haps))+sum(haps)*alt_expr)
            sim_count=simulate_BNB(ind_mean_counts,gene_disp,ind_disp)

            if haps[0]!=haps[1]:
                if random()<het_error_rate:
                    if haps[0]==0:
                        ref,alt=simulate_BB(as_counts, read_error_rate, AS_disp)
                    else:
                        ref,alt=simulate_BB(as_counts, 1-read_error_rate, AS_disp)
                else:
                    ref,alt=simulate_BB(as_counts, AS_frac, AS_disp)
            else:
                if simulate_hom_AS:
                    ref,alt=simulate_BB(as_counts, 0.5, AS_disp)
                else:
                    ref,alt=0,0 #
            snps.append(TestSNP(effect,test,haps,sim_count,ref,alt,1-het_error_rate))
            counts.append(sim_count)
        mc=np.mean(counts)
        Y=[]
        Y1=[]
        Y2=[]
        Z=[]
        
        if num_hets>=min_hets:
            for snp_indx in range(len(snps)):
                snps[snp_indx].set_total_counts(mc)
                out_files[snp_indx].write(snps[snp_indx].print_snp())
                out_files[snp_indx].flush()
                Y.append(snps[snp_indx].count)
                Y1.append(snps[snp_indx].ref_count)
                Y2.append(snps[snp_indx].alt_count)
            
                if(snps[snp_indx].haps[0]==0 and snps[snp_indx].haps[1]==0):
                    Z.append(0)
                elif(snps[snp_indx].haps[0]==0 and snps[snp_indx].haps[1]==1):
                    Z.append(1)
                elif(snps[snp_indx].haps[0]==1 and snps[snp_indx].haps[1]==0):
                    Z.append(1)
                elif(snps[snp_indx].haps[0]==1 and snps[snp_indx].haps[1]==1):
                    Z.append(4)
                    
            ASseq_Y.write("\t".join(str(y) for y in Y)+"\n")
            ASseq_Y1.write("\t".join(str(y1) for y1 in Y1)+"\n")
            ASseq_Y2.write("\t".join(str(y2) for y2 in Y2)+"\n")    
            ASseq_Z.write("\t".join(str(z) for z in Z)+"\n")
            test+=1

class TestSNP:
    def __init__(self,effect,test_num,haps,count,as_ref,as_alt,hetp):
        self.chrm=effect
        self.pos=test_num
        self.ref_allele="A"
        self.alt_allele="T"
        self.count=count
        self.ref_count=as_ref
        self.alt_count=as_alt
        self.haps=haps
        if haps[0]!=haps[1]:
            self.hetp=hetp
        else:
            self.hetp=0
        self.genotype=sum(haps)
        self.tot_count=0

    def set_total_counts(self,tot_count):
        self.tot_count=tot_count

    def print_snp(self):
        return("%i %i %i %s %s %i %i|%i %i %i %i %f %i %i %i %i %i %i\n" % (self.chrm,self.pos,self.pos,self.ref_allele,self.alt_allele,self.genotype,self.haps[0],self.haps[1],self.pos,self.pos+1,self.pos,self.hetp,1,self.ref_count,self.alt_count,0,self.count,self.tot_count))

def simulate_BNB(mean,sigma,n):
    #sys.stderr.write("%f %f\n"%(n,mean))
    mean_p=np.float64(n)/(n+mean)
    sigma=(1/sigma)**2
    a=mean_p*(sigma)+1
    b=(1-mean_p)*sigma
    
    p=beta(a,b)
    #sys.stderr.write("%f %f\n"%(n,p))
    counts=negative_binomial(n,p)
    return counts

def simulate_BB(tot,mean_p,sigma):
    a = mean_p*(1/sigma**2 - 1)
    b = (1-mean_p)*(1/sigma**2 - 1)

    p=beta(a,b)
    counts=binomial(tot,p)
    #sys.stderr.write("%f %f %i\n"%(mean_p,p,counts))
    return counts,(tot-counts)

main()
