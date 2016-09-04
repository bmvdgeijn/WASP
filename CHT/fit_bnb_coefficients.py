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
from scipy.special import gammaln
from scipy.special import betaln
import scipy.stats
import time
import numpy as np

import cht_data


MIN_GENE_FIT = 0.0
MAX_GENE_FIT = 1e8

MIN_GW_FIT = 1.0
MAX_GW_FIT = 10e3


# MIN_GENE_FIT = 0.0
# MAX_GENE_FIT = 1.0

# MIN_GW_FIT = 1.0
# MAX_GW_FIT = 1e3


# bounds for parameters
MIN_MEAN_FIT = MIN_GENE_FIT
MAX_MEAN_FIT = MAX_GENE_FIT

GW_XTOL = 1e-2
GENE_XTOL = 1e-7
MEAN_XTOL = 1e-7




def parse_options():

    default_sample = 2000

    parser=argparse.ArgumentParser(description="This program estimates the "
                                   "genome-wide dispersion parameters for "
                                   "the read depth (beta-negative binomial, BNB) part "
                                   "of the combined haplotype test. These "
                                   "parameter estimates are not required "
                                   "for the allele-specific (beta-binomial) part "
                                   "of the test."
                                   "This program uses an iterative "
                                   "maximum likelihood approach to estimate "
                                   "parameters. To speed up parameter estimation "
                                   "a subset of test regions are randomly "
                                   "selected from the input file. By default %d"
                                   " regions are selected. At the end of each "
                                   "iteration, the current parameter estimates "
                                   "are written to stdout and to the "
                                   "specified output file. "
                                   "If you are tired of waiting for the "
                                   "program to converge, these intermediate "
                                   "estimates can be used before the program "
                                   "has finished running." % default_sample)




    parser.add_argument("infile_list", action='store', default=None)

    parser.add_argument("outfile", action='store', default=None)

    parser.add_argument("--min_as_counts", action='store', dest='min_as_counts',
                        type=int, default=0, metavar="MIN_COUNTS",
                        help="only use regions where total allele-specific "
                        "read counts across individuals > MIN_COUNTS")

    parser.add_argument("--min_counts", action='store', dest='min_counts',
                        type=int, default=0, metavar="MIN_COUNTS",
                        help="only use regions where total number of "
                        "read counts across individuals > MIN_COUNTS")

    parser.add_argument("--skip", action='store', dest='skip',
                        type=int, default=0,
                        help="skip n test region between each one used for fitting")

    parser.add_argument("--sample", type=int, default=default_sample,
                        help="Randomly sample a specified number of test "
                        "regions from the total set to speed up the "
                        "estimation procedure. This is set to %d "
                        "by default. To "
                        "disable sampling, set this to 0 (--sample 0), or use "
                        "the --no_sample option" % default_sample)

    parser.add_argument("--no_sample", action="store_true", dest="no_sample",
                        help="Disable random sampling of test regions",
                        default=False)

    parser.add_argument("--seed", help="Random seed which affects random "
                        "sampling of test regions. Useful for testing.",
                        type=int,
                default=-1)

    parser.add_argument("--fix_gene", type=float, default=-1.0,
                        help="fix the per-gene dispersion parameter estimates to "
                        "specified value, and only fit the genome-wide dispersion "
                        "parameter estimates")

    parser.add_argument("--fix_mean", type=float, default=-1.0,
                        help="fix the per-gene mean dispersion parameter estimates to "
                        "specified value, and only fit the genome-wide dispersion "
                        "parameter estimates")

    options = parser.parse_args()

    if options.no_sample or options.sample < 0:
        # do not perform sampling
        options.sample = 0

    return options





def main():
    options = parse_options()

    # set random seed:
    if options.seed >= 0:
        np.random.seed(seed=options.seed)

    # read input data
    sys.stderr.write("reading input data\n")
    count_matrix, expected_matrix = \
        cht_data.read_count_matrices(options.infile_list,
                                     shuffle=False,
                                     skip=options.skip,
                                     min_counts=options.min_counts,
                                     min_as_counts=options.min_as_counts,
                                     sample=options.sample)
    
    gene_fits = [np.float64(0.005)] * count_matrix.shape[0]
    mean_fits = [np.float64(1)] * count_matrix.shape[0]
    gw_fits = [np.float64(100)] * count_matrix.shape[1]

    old_ll = np.float64(1000000000000000000000000)

    fix_gene = options.fix_gene > 0.0
    fix_mean = options.fix_mean > 0.0

    if fix_gene:
        gene_fits = [np.float64(options.fix_gene)] * count_matrix.shape[0]
    if fix_mean:
        mean_fits = [np.float64(options.fix_mean)] * count_matrix.shape[0]

    # iteratively search for maximum likelihood parameter estimates
    iteration = 1
    while True:
        sys.stderr.write("\niteration %d\n" % iteration)

        # first fit over dispersion params for each region
        sys.stderr.write("fitting per-region overdispersion params\n")
        t1 = time.time()
        #gene_fits, mean_fits, fit_ll = get_gene_overdisp(count_matrix, expected_matrix,
        #                                                 gw_fits, gene_fits, mean_fits,
        #                                                 fix_gene=fix_gene,
        #                                                 fix_mean=fix_mean)

        gene_fit = gene_fits[0]
        gene_fit, mean_fits, fit_ll = get_single_param_gene_overdisp(count_matrix, expected_matrix,
                                                                     gw_fits, gene_fit, mean_fits)

        gene_fits[:] = [gene_fit] * len(gene_fits)


        time_taken = time.time() - t1
        sys.stderr.write("time: %.2fs\n" % time_taken)

        sys.stderr.write("min(gene_fits): %g\n" % np.min(gene_fits))
        sys.stderr.write("max(gene_fits): %g\n" % np.max(gene_fits))
        sys.stderr.write("mean(gene_fits): %g\n" % np.mean(gene_fits))

        sys.stderr.write("min(mean_fits): %g\n" % np.min(mean_fits))
        sys.stderr.write("max(mean_fits): %g\n" % np.max(mean_fits))
        sys.stderr.write("mean(mean_fits): %g\n" % np.mean(mean_fits))


        # then fit genome-wide overdispersion params for each individual
        sys.stderr.write("fitting genome-wide overdispersion params\n")
        t1 = time.time()
        gw_fits, fit_ll = get_gw_overdisp(count_matrix, expected_matrix,
                                          gw_fits, gene_fits, mean_fits)
        time_taken = time.time() - t1
        sys.stderr.write("time: %.2fs\n" % time_taken)
        sys.stderr.write("LL=-%f\n" % fit_ll)
        gw_str = " ".join("%.2f" % x for x in gw_fits)
        sys.stderr.write("current genome-wide overdispersion param estimates:\n"
                         "  %s\n" % gw_str)


        # write estimates to outfile each iteration
        # (this way outfile is written even if optimization terminated by user)
        outfile = open(options.outfile, "w")
        for i in gw_fits:
            outfile.write("%f\n" % i)
        outfile.close()

        iteration += 1

        if old_ll - fit_ll < 0.0001:
            # improvement is small, stop iterations
            break
        old_ll = fit_ll

    sys.stderr.write("done!\n")





def get_gene_overdisp(count_matrix, expected_matrix,
                      gw_fits, gene_fits, mean_fits, iteration=0,
                      fix_gene=False, fix_mean=False):
    fit_ll = 0
    gene_fit_improved = 0
    mean_fit_improved = 0

    for gene_indx in range(count_matrix.shape[0]):

        if not fix_mean:
            cur_like = mean_like(mean_fits[gene_indx],
                                 count_matrix[gene_indx,:],
                                 expected_matrix[gene_indx,:],
                                 gw_fits, gene_fits[gene_indx])

            xtol = min(mean_fits[gene_indx] * 1e-4, MEAN_XTOL)

            res = minimize_scalar(mean_like,
                                  bounds=(MIN_MEAN_FIT, MAX_MEAN_FIT),
                                  args=(count_matrix[gene_indx,:],
                                        expected_matrix[gene_indx,:],
                                        gw_fits, gene_fits[gene_indx]),
                                  options={"xatol" : xtol},
                                  method="Bounded")

            like_diff = cur_like - res.fun
            if like_diff >= 0.0:
                # update parameter
                mean_fits[gene_indx] = res.x
                mean_fit_improved += 1
            else:
                # likelihood got worse indicating failed to converge,
                # do not accept new param value
                pass

        if not fix_gene:
            cur_like = gene_like(gene_fits[gene_indx], count_matrix[gene_indx,:],
                                 expected_matrix[gene_indx,:], gw_fits,
                                 mean_fits[gene_indx])

            xtol = min(gene_fits[gene_indx] * 1e-4, GENE_XTOL)

            res = minimize_scalar(gene_like,
                                  bounds=(MIN_GENE_FIT, MAX_GENE_FIT),
                                  args=(count_matrix[gene_indx,:],
                                        expected_matrix[gene_indx,:],
                                        gw_fits, mean_fits[gene_indx]),
                                  options={"xatol" : xtol},
                                  method="Bounded")

            like_diff = cur_like - res.fun
            if like_diff >= 0.0:
                # update parameter
                gene_fits[gene_indx] = res.x
                gene_fit_improved += 1
            else:
                # likelihood got worse indicating failed to converge,
                # do not accept new param value
                pass


        fit_ll += gene_like(gene_fits[gene_indx],
                            count_matrix[gene_indx,:],
                            expected_matrix[gene_indx,:],
                            gw_fits, mean_fits[gene_indx])

    return gene_fits, mean_fits, fit_ll




def get_single_param_gene_overdisp(count_matrix, expected_matrix,
                                   gw_fits, gene_fit, mean_fits, iteration=0,
                                   fix_gene=False, fix_mean=False):
    fit_ll = 0
    gene_fit_improved = 0
    mean_fit_improved = 0

    # update mean fits
    for gene_indx in range(count_matrix.shape[0]):

        if not fix_mean:
            cur_like = mean_like(mean_fits[gene_indx],
                                 count_matrix[gene_indx,:],
                                 expected_matrix[gene_indx,:],
                                 gw_fits, gene_fit)

            xtol = min(mean_fits[gene_indx] * 1e-4, MEAN_XTOL)

            res = minimize_scalar(mean_like,
                                  bounds=(MIN_MEAN_FIT, MAX_MEAN_FIT),
                                  args=(count_matrix[gene_indx,:],
                                        expected_matrix[gene_indx,:],
                                        gw_fits, gene_fit),
                                  options={"xatol" : xtol},
                                  method="Bounded")

            like_diff = cur_like - res.fun
            if like_diff >= 0.0:
                # update parameter
                mean_fits[gene_indx] = res.x
                mean_fit_improved += 1
            else:
                # likelihood got worse indicating failed to converge,
                # do not accept new param value
                pass

    cur_like = single_param_gene_like(gene_fit, count_matrix,
                                      expected_matrix, gw_fits,
                                      mean_fits)

    xtol = min(gene_fit * 1e-4, GENE_XTOL)

    res = minimize_scalar(single_param_gene_like,
                          bounds=(MIN_GENE_FIT, MAX_GENE_FIT),
                          args=(count_matrix,
                                expected_matrix,
                                gw_fits, mean_fits),
                          options={"xatol" : xtol},
                          method="Bounded")

    like_diff = cur_like - res.fun

    if like_diff >= 0.0:
        # update parameter
        gene_fit = res.x
        fit_ll = res.fun
    else:
        # likelihood got worse indicating failed to converge,
        # do not accept new param value
        fit_ll = cur_like

    return gene_fit, mean_fits, fit_ll




def get_gw_overdisp(count_matrix, expected_matrix, gw_fits,
                    gene_fits, mean_fits):
    fit_ll = 0

    for indx in range(count_matrix.shape[1]):
        cur_like = gw_like(gw_fits[indx], count_matrix[:,indx],
                           expected_matrix[:,indx], gene_fits, mean_fits)

        init_gw_fit = gw_fits[indx]

        xtol = min(gw_fits[indx] * 1e-2, GW_XTOL)

        res = minimize_scalar(gw_like,
                              bounds=(MIN_GW_FIT, MAX_GW_FIT),
                              args=(count_matrix[:,indx],
                                    expected_matrix[:,indx],
                                    gene_fits, mean_fits),
                              options={'xatol' : xtol},
                              # tol=GW_LIKE_TOL,
                              method="Bounded")

        new_like = res.fun
        like_diff = cur_like - res.fun

        if like_diff >= 0.0:
            # likelhood improved
            gw_fits[indx] = res.x
            like = new_like
        else:
            # likelihood did not improve because failed to converge
            like = cur_like

        like = gw_like(gw_fits[indx], count_matrix[:,indx],
                       expected_matrix[:,indx], gene_fits, mean_fits)

        fit_ll += like
    return gw_fits, fit_ll



def mean_like(mean_fit, counts, expecteds, gw_fits, gene_fit):
    if mean_fit < MIN_MEAN_FIT or mean_fit > MAX_MEAN_FIT:
        return 1e8
    loglike=0
    for i in range(len(counts)):
        loglike += BNB_loglike(counts[i], expecteds[i]*mean_fit,
                               gw_fits[i], gene_fit)
    return -loglike



def gene_like(gene_fit, counts, expecteds, gw_fits, mean_fit):
    if gene_fit < MIN_GENE_FIT or gene_fit > MAX_GENE_FIT:
        return 1e8
    loglike=0
    for i in range(len(counts)):
        loglike += BNB_loglike(counts[i], expecteds[i]*mean_fit,
                               gw_fits[i], gene_fit)
    return -loglike


def single_param_gene_like(gene_fit, count_matrix, expected_matrix, gw_fits, mean_fits):
    ll = 0

    for gene_indx in range(count_matrix.shape[0]):
        ll += gene_like(gene_fit,
                        count_matrix[gene_indx,:],
                        expected_matrix[gene_indx,:],
                        gw_fits, mean_fits[gene_indx])

    return ll



def gw_like(gw_fit, counts, expecteds, gene_fits, mean_fits):
    if gw_fit > MAX_GW_FIT or gw_fit < MIN_GW_FIT:
        return 1e8
    loglike = 0
    for i in range(len(counts)):
        loglike += BNB_loglike(counts[i], expecteds[i]*mean_fits[i],
                               gw_fit, gene_fits[i])
    return -loglike





def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))



def lbeta_asymp(a,b):
    if b > a:
        a,b = b,a

    if a<1e6:
        return betaln(a,b)

    l = gammaln(b)

    l -= b*math.log(a)
    l += b*(1-b)/(2*a)
    l += b*(1-b)*(1-2*b)/(12*a*a)
    l += -((b*(1-b))**2)/(12*a**3)

    return l



def BNB_loglike(k,mean,n,sigma):
    #n=min(n,10000)
    #Put variables in beta-NB form (n,a,b)
    mean = max(mean, 0.00001)
    p = np.float64(n)/(n+mean)
    logps = [math.log(n) - math.log(n + mean),
             math.log(mean) - math.log(n + mean)]

    if sigma < 0.00001:
        loglike = -betaln(n,k+1)-math.log(n+k)+n*logps[0]+k*logps[1]
        return loglike

    sigma = (1/sigma)**2

    a = p * sigma + 1
    b = (1-p) * sigma

    loglike = 0

    #Rising Pochhammer = gamma(k+n)/gamma(n)
    #for j in range(k):
    #    loglike += math.log(j+n)
    if k>0:
        loglike = -lbeta_asymp(n,k) - math.log(k)
        #loglike=scipy.special.gammaln(k+n)-scipy.special.gammaln(n)
    else:
        loglike=0

    #Add log(beta(a+n,b+k))
    loglike += lbeta_asymp(a+n,b+k)

    #Subtract log(beta(a,b))
    loglike -= lbeta_asymp(a,b)

    return loglike



main()
