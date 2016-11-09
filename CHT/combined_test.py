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
from scipy.special import gammaln
from scipy.special import betaln
import scipy.stats

import numpy as np

import random

import util

# OPTIMIZER="BFGS"
OPTIMIZER="Nelder-Mead"


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
        if util.is_gzipped(filename):
            f = gzip.open(filename, "rt")
        else:
            f = open(filename, "r")

        # skip header
        f.readline()

        infiles.append(f)
    in_file.close()

    if len(infiles) == 0:
        sys.stderr.write("no input files specified in file '%s'\n" % options.infile_list)
        exit(2)

    return infiles



def write_header(outfile):
    outfile.write("\t".join(["TEST.SNP.CHROM", "TEST.SNP.POS",
                             "LOGLIKE.NULL", "LOGLIKE.ALT",
                             "CHISQ", "P.VALUE", "ALPHA", "BETA",
                             "PHI", "TOTAL.AS.READ.COUNT",
                             "REF.AS.READ.COUNT", "ALT.AS.READ.COUNT",
                             "TOTAL.READ.COUNT"]) + "\n")


def read_bnb_sigmas(options, infiles):
    """Read overdispersion parameters for beta-negative binomial.
    Expect one for each individual."""
    if (options.bnb_disp):
        disp_file = open(options.bnb_disp)
        line = disp_file.readline()
        bnb_sigmas = []
        while line:
            bnb_sigmas.append(np.float64(line.strip()))
            line = disp_file.readline()
        disp_file.close()

        if len(bnb_sigmas) != len(infiles):
            raise ValueError("expected %d values in bnb_disp file "
                             "(one for each input file) but got %d"
                             % (len(infiles), len(bnb_sigmas)))
    else:
        bnb_sigmas = [0.001]*len(infiles)

    return bnb_sigmas



def read_as_sigmas(options, infiles):
    """Read overdispersion parameters for allele-specific test
    (Beta-Binomial). Expect one for each individual."""

    if (options.as_disp):
        disp_file = open(options.as_disp)
        line = disp_file.readline()
        as_sigmas = []
        while line:
            val = np.float64(line.strip())
            if val < 0.0 or val > 1.0:
                raise ValueError("expected as_sigma values to be "
                                 " in range 0.0-1.0, but got %g" %
                                 val)
            as_sigmas.append(np.float64(line.strip()))
            line = disp_file.readline()

        disp_file.close()

        if len(as_sigmas) != len(infiles):
            raise ValueError("expected %d values in as_disp file "
                             "(one for each input file) but got "
                             "%d" % (len(infiles), len(as_sigmas)))

    else:
        as_sigmas = [0.001] * len(infiles)

    return as_sigmas



def write_results(outfile, snpinfo, loglike1par, loglike2par,
                  best2par, tot_as_counts, ref_as_counts, alt_as_counts,
                  all_counts):
    """Write result to output file. Tab-delimited columns are:
      1. chromosome,
      2. SNP position,
      3. Log likelihood 1 parameter model (Null)
      4. Log likelihood 2 parameter model (Alternative)
      3. Chi-squared statistic,
      4. P-value
      5. alpha parameter estimate (expression level
         of reference allele)
      6. beta parameter estimate (expression level of
         alternative allele)
      7. phi parameter estimate (beta-negative-binomial
         overdispersion
         parameter for this region)
      8. total number of allele-specific read counts for this
         region summed across individuals
      9. total number of reference haplotype allele-specific read counts
     10. total number of alt haplotype allele-specific read counts
     11. total number of mapped reads for this region,
         summed across individuals"""

    # compute likelihood ratio test statistic:
    chisq = 2 * (loglike1par - loglike2par)
    pval = (1-scipy.stats.chi2.cdf(chisq,1)),

    outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1],
                             "%.2f" % -loglike1par,
                             "%.2f" % -loglike2par,
                             "%.3f" % chisq,
                             "%g" % pval,
                             "%g" % best2par[0],
                             "%g" % best2par[1],
                             "%g" % best2par[2],
                             "%d" % tot_as_counts,
                             "%d" % ref_as_counts,
                             "%d" % alt_as_counts,
                             "%d" % all_counts]) + '\n')
    outfile.flush()


def write_empty_result(outfile, snpinfo):
    """Write all zeros in the even that the test failed"""
    outfile.write("\t".join([snpinfo[0][0], snpinfo[0][1], "0", "0",
                             "0", "NA", "0", "0", "0",
                             "0", "0", "0"]) + '\n')


def rescale_totals(test_snps):
    min_tot = min([s.totals for s in test_snps])

    if min_tot > 0:
        for s in test_snps:
            s.totals = s.totals / min_tot        
        
    


def main():
    options = parse_options()

    if options.pc_file:
        pc_matrix = load_covariates(options.pc_file)
        num_pcs = options.num_pcs
    else:
        pc_matrix = []
        num_pcs = 0

    if options.out_file.endswith(".gz"):
        outfile = gzip.open(options.out_file, "wb")
    else:
        outfile = open(options.out_file, 'w')

    if options.benchmark:
        if options.benchmark == "-":
            bench_file = sys.stderr
        else:
            bench_file = open(options.benchmark, "w")

        bench_file.write("TEST.TYPE TIME\n")

    write_header(outfile)

    # read list of input files (one for each individual)
    infiles = open_input_files(options.infile_list)

    # read dispersion parameters for each individual
    bnb_sigmas = read_bnb_sigmas(options, infiles)
    as_sigmas = read_as_sigmas(options, infiles)


    # add first row of each input file to snpinfo list
    snpinfo = []
    for f in infiles:
        snpinfo.append(f.readline().strip().split())

    row_count = 0
    finished=False

    options.dup_snp_warn = True
    
    while not finished:
        try:
            test_snps = []
            # parse test SNP and associated info from input file row
            for i in range(len(infiles)):
                test_snps.append(parse_test_snp(snpinfo[i], options))

            # rescale totals to put values into reasonable range for
            # alpha and beta parameter estimation
            rescale_totals(test_snps)
            

            # how many allele-specific reads are there across all
            # linked SNPs and and individuals?
            ref_as_counts = sum([np.sum(x.AS_target_ref) for x in test_snps])
            alt_as_counts = sum([np.sum(x.AS_target_alt) for x in test_snps])
            tot_as_counts = ref_as_counts + alt_as_counts
            
            all_counts = sum([test_snps[i].counts for i in range(len(test_snps))])

            if tot_as_counts < options.min_as_counts:
                if options.verbose:
                    sys.stderr.write("-----\nskipping SNP %s because "
                                     "total AS counts %d <= %d\n" %
                                     (test_snps[0].name, tot_as_counts, options.min_as_counts))

                # skip, not enough allele-specific counts
                for i in range(len(infiles)):
                    line = infiles[i].readline().strip()
                    if line:
                        snpinfo[i] = line.split()
                    else:
                        # out of lines from at least one file, assume we are finished
                        finished = True
                continue

            if options.verbose:
                sys.stderr.write("-----\ntesting SNP %s\n" % test_snps[0].name)

            row_count+=1
            old_genos = [test_snps[y].geno_hap1 + test_snps[y].geno_hap2
                         for y in range(len(test_snps))]

            if options.shuffle:
                # permute genotypes
                perm = range(len(test_snps))
                random.shuffle(perm)
                geno1temp = [test_snps[y].geno_hap1 for y in perm]
                geno2temp = [test_snps[y].geno_hap2 for y in perm]
                for i in range(len(test_snps)):
                    test_snps[i].geno_hap1 = geno1temp[i]
                    test_snps[i].geno_hap2 = geno2temp[i]


            if options.benchmark:
                # start timing test for NULL model
                bench_file.write("null %s %s %d %d " % (snpinfo[0][0], snpinfo[0][1],
                                                        tot_as_counts, all_counts))
                bench_file.flush()
            t1 = time.time()

            starting_gene = [np.float64(x) for x in [0.1, 0.001]]
            maxlike = 10000000000

            for start in starting_gene:
                starts = [np.float64(0.5), np.float64(start)]

                # regress against the covariates and get residuals
                #fit_cov(test_snps,cov_table)

                # maximize likelihood with alpha = beta (no difference between genotypes)
                res = minimize(ll_one, starts, args=(test_snps, True, #options.is_bnb_only,
                                                     options.is_as_only,
                                                     bnb_sigmas,
                                                     as_sigmas,
                                                     options.read_error_rate,
                                                     [],
                                                     pc_matrix),
                                options={"maxiter" : 50000, "disp" : options.verbose},
                                method=OPTIMIZER)

                new_par = res.x

                new_loglike = ll_one(new_par, test_snps, options.is_bnb_only,
                                     options.is_as_only, bnb_sigmas,
                                     as_sigmas, options.read_error_rate,
                                     [], pc_matrix)
                if new_loglike < maxlike:
                    starting_par = new_par

            pc_coefs=[]
            for pc in range(num_pcs):
                res = minimize(ll_pc, [np.float64(0)],
                               args=(starting_par, test_snps, True, #options.is_bnb_only,
                                     options.is_as_only, bnb_sigmas, as_sigmas,
                                     options.read_error_rate, pc_coefs, pc_matrix),
                               options={"maxiter" : 50000, "disp" : options.verbose},
                               method=OPTIMIZER)

                new_coef = res.x
                pc_coefs = np.concatenate([pc_coefs, new_coef])

            res = minimize(ll_one, starting_par,
                           args=(test_snps, options.is_bnb_only, options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate, pc_coefs, pc_matrix),
                           options={"maxiter" : 50000, "disp" : options.verbose},
                           method=OPTIMIZER)
            best1par = res.x

            time_taken = time.time() - t1
            if options.verbose:
                sys.stderr.write("null model optimization took %.3fs\n" % time_taken)
            if options.benchmark:
                bench_file.write("%.3f\n" % time_taken)
                bench_file.flush()


            loglike1par = ll_one(best1par, test_snps, options.is_bnb_only,
                                 options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate,
                                 pc_coefs, pc_matrix)


            start = [best1par[0], best1par[0], best1par[1]]

            if options.benchmark:
                # start timing test for ALT model
                bench_file.write("alt %s %s %d %d " % (snpinfo[0][0], snpinfo[0][1],
                                                       tot_as_counts, all_counts))
                bench_file.flush()

            t1 = time.time()

            # maximize likelihood with alpha and beta as separate parameters
            res = minimize(ll_two, start,
                           args=(test_snps, options.is_bnb_only, options.is_as_only,
                                 bnb_sigmas, as_sigmas, options.read_error_rate,
                                 pc_coefs, pc_matrix),
                           options={"maxiter" : 50000, "disp" : options.verbose},
                           method=OPTIMIZER)
            best2par = res.x

            time_taken = time.time() - t1
            if options.verbose:
                sys.stderr.write("alternative model optimization took %.3fs\n" % time_taken)
            if options.benchmark:
                bench_file.write("%.3f\n" % time_taken)
                bench_file.flush()

            loglike2par = ll_two(best2par, test_snps, options.is_bnb_only,
                                 options.is_as_only, bnb_sigmas,
                                 as_sigmas, options.read_error_rate,
                                 pc_coefs, pc_matrix)

            write_results(outfile, snpinfo, loglike1par, loglike2par, best2par,
                          tot_as_counts, ref_as_counts, alt_as_counts, all_counts)

        except Exception as e:
            write_empty_result(outfile, snpinfo)
            # an error occured, write to output file, but put 0s for all params and
            sys.stderr.write("An error occurred, writing line with 0s for SNP:\n%s\n" % str(e))
            # continue
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
    parser.add_argument("-a", "--as_only",
                        action='store_true',
                        dest='is_as_only', default=False,
                        help="only perform the allele-specific part (Beta Binomial) "
                        "part of the test")

    parser.add_argument("-d", "--bnb_only", action='store_true',
                        dest='is_bnb_only', default=False,
                        help="only perform the association (Beta Negative Binomial) part "
                        "of the test")

    parser.add_argument("--pc_file", action='store',
                        dest='pc_file',
                        help="file containing PC covariates to include in the model"
                        ,default=None)

    parser.add_argument("-b", "--bnb_disp", action='store', dest='bnb_disp',
                        help="file containing depth (Beta Negative Binomial)"
                        "dispersion parameters", default=None)

    parser.add_argument("-o", "--as_disp", action='store',
                        dest='as_disp',
                        help="file containing allele-specific (Beta Binomial) dispersion "
                        "parameters", default=None)

    parser.add_argument("-s", "--shuffle", action='store_true',
                        dest='shuffle', default=False,
                        help="permute genotypes")

    parser.add_argument("-e", "--read_error_rate", action='store', dest='read_error_rate',
                        help="estimate of error rate, used to update "
                        "heterozygous genotype probabilities "
                        "(currently this option disabled / not used)",
                        type=float, default=0.005)

    parser.add_argument("-m", "--min_as_counts", action='store', dest='min_as_counts',
                        type=int, default=0,
                        help="only perform test when total number of allele-specific "
                        "read counts across individuals > MIN_COUNTS")

    parser.add_argument("--num_pcs", action='store', dest='num_pcs',
                        type=int, default=0,
                        help="designates the number of PCs to use as covariates")

    parser.add_argument("-v", "--verbose", action='store_true', dest='verbose',
                        default=False, help="print extra information")

    parser.add_argument("--benchmark", dest="benchmark",
                        help="write information about time test is takes, number of optimization "
                        "functions, etc. to specified filename, or to stderr if '-' is specified")

    parser.add_argument("infile_list", action='store', default=None)
    parser.add_argument("out_file", action='store', default=None)

    return parser.parse_args()




def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))



def AS_betabinom_loglike(logps, sigma, AS1, AS2, hetp, error):
    """Given parameter, returns log likelihood of allele-specific
    part of test. Note that some parts of equation have been
    canceled out"""
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
        mean = max(mean,0.00001)
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



def ll_one(x, test_snps, is_bnb_only, is_as_only, bnb_sigmas,
           as_sigmas, error, pc_coefs, pc_matrix):
    alpha = x[0]
    beta = x[0]
    r = x[1]
    return loglikelihood(alpha, beta, r, test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, pc_coefs, pc_matrix)


def ll_pc(x, params, test_snps, is_bnb_only, is_as_only, bnb_sigmas,
          as_sigmas, error, other_pc_coefs, pc_matrix):
    alpha = params[0]
    beta = params[0]
    r = params[1]
    pc_coefs=np.concatenate([other_pc_coefs,x])
    return loglikelihood(alpha, beta, r, test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, pc_coefs,pc_matrix)


def ll_two(x, test_snps, is_bnb_only, is_as_only, bnb_sigmas,
           as_sigmas, error, pc_coefs, pc_matrix):
    alpha = x[0]
    beta = x[1]
    r = x[2]
    #if len(x)>3:
    #    pc_fits=x[3:]
    #else:
    #    pc_fits=[]
    return loglikelihood(alpha, beta, r, test_snps, is_bnb_only,
                         is_as_only, bnb_sigmas, as_sigmas, error, pc_coefs, pc_matrix)


def calc_pc_factor(pc_fits, pcs, i):
    if len(pc_fits) > 0:
        return 1 + sum(pc_fits * pcs[i,:len(pc_fits)])
    else:
        return 1


def loglikelihood(alpha, beta, r, test_snps, is_bnb_only,
                  is_as_only, bnb_sigmas, as_sigmas, error,
                  pc_coefs, pc_matrix):
    loglike = 0

    # if input values are outside of reasonable range return a
    # very high -loglike
    if alpha <= 0 or beta <= 0 or r <= 0 or r > 1:
        return 10000000

    ratio = (alpha / (alpha + beta))

    for i in range(len(test_snps)):
        if(test_snps[i].is_homo_ref()):
            m = 2*alpha*test_snps[i].totals * calc_pc_factor(pc_coefs, pc_matrix, i)
        elif(test_snps[i].is_homo_alt()):
            m = 2*beta*test_snps[i].totals * calc_pc_factor(pc_coefs, pc_matrix, i)
        else:
            m = (alpha+beta)*test_snps[i].totals * calc_pc_factor(pc_coefs, pc_matrix, i)
        if m<0:
            m = 0.000001
        if not is_bnb_only:
            for j in range(len(test_snps[i].AS_target_ref)):
                if test_snps[i].hetps[j]>.9:
                    hetp = min(0.99, test_snps[i].hetps[j])
                    logps = [math.log(alpha) - math.log(alpha+beta),
                             math.log(beta) - math.log(alpha+beta)]
                    loglike += AS_betabinom_loglike(logps, as_sigmas[i],
                                                    test_snps[i].AS_target_ref[j],
                                                    test_snps[i].AS_target_alt[j],
                                                    hetp, error)
        if not is_as_only:
            l = BNB_loglike(test_snps[i].counts, m, r, bnb_sigmas[i])
            loglike += l
    return -loglike


def parse_test_snp(snpinfo, options):    
    snp_id = snpinfo[2]
    if snpinfo[16] == "NA":
        # SNP is missing data
        tot = 0
    else:
        # these totals are later rescaled by dividing
        # by the minimum total across individuals to
        # put them into a reasonable range for
        # estimating alpha and beta
        tot = float(snpinfo[16])

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
        # positions of target SNPs
        snp_locs = np.array([int(y.strip()) for y in snpinfo[9].split(';')])

        # counts of reads that match reference overlapping linked 'target' SNPs
        snp_as_ref = np.array([int(y) for y in snpinfo[12].split(';')])

        # counts of reads that match alternate allele
        snp_as_alt = np.array([int(y) for y in snpinfo[13].split(';')])

        # heterozygote probabilities
        snp_hetps = np.array([np.float64(y.strip())
                          for y in snpinfo[10].split(';')])

        # linkage probabilities, not currently used
        snp_linkageps = np.array([np.float64(y.strip())
                                  for y in snpinfo[11].split(';')])


        # same SNP should not be provided multiple times, this
        # can create problems with combined test. Warn and filter
        # duplicate SNPs
        uniq_loc, uniq_idx = np.unique(snp_locs, return_index=True)

        if options.dup_snp_warn and uniq_loc.shape[0] != snp_locs.shape[0]:
            sys.stderr.write("WARNING: discarding SNPs that are repeated "
                                     "multiple times in same line\n")
            options.dup_snp_warn = False

        snp_as_ref = snp_as_ref[uniq_idx]
        snp_as_alt = snp_as_alt[uniq_idx]
        snp_hetps = snp_hetps[uniq_idx]
        snp_linkageps = snp_linkageps[uniq_idx]
                             
        if options.shuffle:
            # permute allele-specific read counts by flipping them randomly at
            # each SNP
            for y in range(len(snp_as_ref)):
                if random.randint(0, 1) == 1:
                    temp = snp_as_ref[y]
                    snp_as_ref[y] = snp_as_alt[y]
                    snp_as_alt[y] = temp

        return TestSNP(snp_id, geno_hap1, geno_hap2, snp_as_ref,
                       snp_as_alt, snp_hetps, tot, count)

def load_covariates(cov_file):
    infile=open(cov_file)
    cov_table=[]
    while True:
        line=infile.readline()
        if line:
            cov_table.append([np.float64(x) for x in line.strip().split()])
        else:
            break

    return np.array(cov_table, dtype=np.float64)

main()
