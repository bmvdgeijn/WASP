from random import choice, sample, random
from numpy.random import beta, negative_binomial, binomial
import numpy as np
import sys

import argparse


def write_options(f, options):
    f.write("# prefix = %s\n" % options.prefix)
    f.write("# num_tests = %d\n" % options.num_tests)
    f.write("# num_inds = %d\n" % options.num_inds)
    f.write("# min_hets = %d\n" % options.min_hets)
    f.write("# maf = %g\n" % options.maf)
    f.write("# mean_counts = %g\n" % options.mean_counts)
    f.write("# mean_counts_distr = %s\n" % options.mean_counts_distr)
    f.write("# as_counts = %g\n" % options.as_counts)
    f.write("# gene_disp = %g\n" % options.gene_disp)
    f.write("# gene_disp_distr = %s\n" % options.gene_disp_distr)
    f.write("# ind_disp = %s\n" % ",".join(["%g" % x for x in options.ind_disp]))
    f.write("# as_disp = %g\n" % options.as_disp)
    f.write("# effect_size = %g\n" % options.effect_size)
    f.write("# additivity = %g\n" % options.additivity)
    f.write("# het_error_rate = %g\n" % options.het_error_rate)
    f.write("# read_error_rate = %g\n" % options.read_error_rate)
    f.write("# true_positives = %g\n" % options.true_positives)
    f.write("# sim_hom_as = %s\n" % options.sim_hom_as)

    
def write_header(f):
    f.write("CHROM "
            "TEST.SNP.POS "
            "TEST.SNP.ID "
            "TEST.SNP.REF.ALLELE "
            "TEST.SNP.ALT.ALLELE "
            "TEST.SNP.GENOTYPE "
            "TEST.SNP.HAPLOTYPE "
            "REGION.START "
            "REGION.END "
            "REGION.SNP.POS "
            "REGION.SNP.HET.PROB "
            "REGION.SNP.LINKAGE.PROB "
            "REGION.SNP.REF.HAP.COUNT "
            "REGION.SNP.ALT.HAP.COUNT "
            "REGION.SNP.OTHER.HAP.COUNT "
            "REGION.READ.COUNT "
            "GENOMEWIDE.READ.COUNT\n")



def parse_options():
    parser = argparse.ArgumentParser(description="simulate counts for "
                                     "combined haplotype test")

    parser.add_argument("--prefix", default=None, required=True,
                        help="prefix for output files")

    dflt_num_tests = 10000
    parser.add_argument("--num_tests", default=dflt_num_tests, type=int,
                        help="number of regions to simulate "
                        "(default=%d)" % dflt_num_tests)

    dflt_num_inds = 10
    parser.add_argument("--num_inds", default=dflt_num_inds, type=int,
                        help="number of individuals to simulate "
                        "(default=%d)" % dflt_num_inds)

    dflt_min_hets = 2
    parser.add_argument("--min_hets", default=dflt_min_hets, type=int,
                        help="minimum number of heterozygous individuals "
                        "per test SNP (default=%d)" % dflt_min_hets)

    dflt_maf = 0.2
    parser.add_argument("--maf", default=dflt_maf, type=float,
                        help="minor allele frequency of test SNP "
                        "(default=%.2f)" % dflt_maf)

    dflt_mean_counts = 200.0
    parser.add_argument("--mean_counts", default=dflt_mean_counts, type=float,
                        help="mean number of read counts per region "
                        "(default=%.2f)" % dflt_mean_counts)

    dflt_mean_counts_distr = "POINT"
    parser.add_argument("--mean_counts_distr", default=dflt_mean_counts_distr,
                        help="distribution for mean number of "
                        "read counts per region (default=%s). "
                        " If EXPONENTIAL is specified, the value of "
                        "mean_counts is used as the scale parameter "
                        "(mean) of the distribution" % dflt_mean_counts_distr,
                        choices=("POINT", "EXPONENTIAL"))

    dflt_as_counts = 20.0
    parser.add_argument("--as_counts", default=dflt_as_counts, type=float,
                        help="expected number of allele-specific read counts "
                        "per regions (default=%.2f)" % dflt_as_counts)

    dflt_gene_disp = 0.01
    parser.add_argument("--gene_disp", default=dflt_gene_disp, type=float,
                        help="per-gene overdispersion parameter for "
                        "beta-negative binomial (default=%.2f)" %
                        dflt_gene_disp)

    dflt_gene_disp_distr = "POINT"
    parser.add_argument("--gene_disp_distr", default=dflt_gene_disp_distr,
                        help="distribution for sampling per-gene "
                        "overdispersion from (default=%s). If EXPONENTIAL "
                        "is specified, the value of gene_disp is used as "
                        "the scale parameter (mean) of the distribution"
                        % dflt_gene_disp, choices=("POINT", "EXPONENTIAL"))

    dflt_ind_disp = "100.0"
    parser.add_argument("--ind_disp", default=dflt_ind_disp,
                        help="per individual overdispersion parameter(s) for "
                        "beta-negative binomial. Can either provide a single value "
                        "that is used for all individuals or a comma-delimited list "
                        "of values (one per individual) (default=%s)." % dflt_ind_disp)
    
    dflt_as_disp = 0.2
    parser.add_argument("--as_disp", default=dflt_as_disp, type=float,
                        help="per individual allele-specific overdispersion "
                        "parameter for beta-binomial (default=%.2f)" % dflt_as_disp)

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

    dflt_sim_hom_as = False
    parser.add_argument("--sim_hom_as", action="store_true", dest="sim_hom_as",
                        help="simulate allele specific counts "
                        "at homozygous test SNPs (default=False)", default=False)

    options = parser.parse_args()

    # split apart individual dispersion value string
    vals = [float(x) for x in options.ind_disp.split(",")]
    if len(vals) == 1:
        options.ind_disp = [vals[0]] * options.num_inds
    elif len(vals) != options.num_inds:
        raise ValueError("number of ind_disp values should be "
                         "equal to 1 or to num_ind (%d)\n" % options.num_inds)
    else:
        options.ind_disp = vals
        
    return options
    



def main():
    options = parse_options()
    
    out_files = []
    sys.stderr.write("creating output files:\n")
    file_list = open("%s_file_list.txt" % options.prefix, "w")
    for i in range(options.num_inds):
        out_filename = "%s_%d.txt" % (options.prefix, i+1)
        sys.stderr.write("  %s\n" % out_filename)

        out_files.append(open(out_filename, "w"))
        # write_options(out_files[i], options)
        write_header(out_files[i])
        file_list.write(out_filename + "\n")
    file_list.close()
    
    
    ASseq_Y_file = open("%s_Y.txt" % options.prefix, "w")
    ASseq_Y1_file = open("%s_Y1.txt" % options.prefix, "w")
    ASseq_Y2_file = open("%s_Y2.txt" % options.prefix, "w")
    ASseq_Z_file = open("%s_Z.txt" % options.prefix, "w")

    test = 1
    while test <= options.num_tests:
        if random() > options.true_positives:
            # simulate a site with no effect, this is not a positive
            effect = 0
            alt_expr = 1
            AS_frac = 0.5
        elif random() < 0.5:
            # simulate a site with effect and beta > alpha
            effect = 0
            alt_expr = 1.0 + options.effect_size
            AS_frac = 1.0 / (2.0 + options.additivity * options.effect_size)
        else:
            # simulate a site with effect and alpha > beta
            effect = 1
            alt_expr = 1 / (1 + options.effect_size)
            AS_frac = (1 + options.additivity * options.effect_size) / \
                (2.0 + options.additivity * options.effect_size)

        snps = []
        counts = []
        num_hets = 0

        
        if options.mean_counts_distr == "POINT":
            mean_counts = options.mean_counts
        elif options.mean_counts_distr == "EXPONENTIAL":
            mean_counts = np.random.exponential(options.mean_counts)
        else:
            raise ValueError("unknown distribution %s\n" %
                             options.mean_counts_distr)

        if options.gene_disp_distr == "POINT":
            gene_disp = options.gene_disp
        elif options.gene_disp_distr == "EXPONENTIAL":
            gene_disp = np.random.exponential(options.gene_disp)
            sys.stderr.write("gene_disp: %.2f\n" % gene_disp)
        else:
            sys.stderr.write("unknown distribution: %s\n" % gene_disp)

        
        for ind in range(options.num_inds):
            # Simulate the individual's haps=[0,0]
            # prob of each minor allele is MAF (minor allele freq)
            is_het = False

            n_minor = int(random() < options.maf) + int(random() < options.maf)
            if n_minor == 0:
                # no minor alleles
                haps = [0,0]
            elif n_minor == 1:
                # heterozygous
                haps = [0,1]
                num_hets += 1
                is_het = True
            else:
                # two minor alleles
                haps = [1,1]
            
            # Expected number of reads based on genotypes
            ind_mean_counts = mean_counts * ((2 - n_minor) + (n_minor * alt_expr))
            #sys.stderr.write("n_minor: %d alt_expr: %g mean_counts: %g " %
            #                 (n_minor, alt_expr, ind_mean_counts))
            
            sim_count = simulate_BNB(ind_mean_counts, gene_disp, options.ind_disp[ind])

            if is_het:
                if random() < options.het_error_rate:
                    # simulate a homozygous site that was miscalled
                    # as a heterozygote
                    if haps[0] == 0:
                        ref, alt = simulate_BB(options.as_counts,
                                               options.read_error_rate, options.as_disp)
                    else:
                        ref, alt = simulate_BB(options.as_counts, 1-options.read_error_rate,
                                               options.as_disp)
                else:
                    ref, alt = simulate_BB(options.as_counts, AS_frac,
                                           options.as_disp)
            else:
                if options.sim_hom_as:
                    # simulate allele-specific counts even when test SNP
                    # is homozygous
                    ref, alt = simulate_BB(options.as_counts, 0.5, AS_disp)
                else:
                    ref, alt = 0, 0
            snps.append(TestSNP(effect, test, haps, sim_count,
                                ref, alt, 1.0 - options.het_error_rate))
            counts.append(sim_count)

        mean_counts = np.mean(counts)
        Y=[]
        Y1=[]
        Y2=[]
        Z=[]
        
        if num_hets >= options.min_hets:
            for snp_indx in range(len(snps)):
                snps[snp_indx].set_total_counts(mean_counts)
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
                    
            ASseq_Y_file.write("\t".join(str(y) for y in Y)+"\n")
            ASseq_Y1_file.write("\t".join(str(y1) for y1 in Y1)+"\n")
            ASseq_Y2_file.write("\t".join(str(y2) for y2 in Y2)+"\n")    
            ASseq_Z_file.write("\t".join(str(z) for z in Z)+"\n")
            test+=1

class TestSNP:
    def __init__(self,effect,test_num,haps,count,as_ref,as_alt,hetp):
        self.chrm = effect
        self.pos = test_num
        self.ref_allele = "A"
        self.alt_allele = "T"
        self.count = count
        self.ref_count = as_ref
        self.alt_count = as_alt
        self.haps = haps
        if haps[0] != haps[1]:
            self.hetp=hetp
        else:
            self.hetp = 0
        self.genotype = sum(haps)
        self.tot_count = 0

    def set_total_counts(self, tot_count):
        self.tot_count = tot_count

    def print_snp(self):
        return("%i %i %i %s %s %i %i|%i %i %i %i %f %i %i %i %i %i %i\n" %
               (self.chrm, self.pos, self.pos, self.ref_allele, self.alt_allele,
                self.genotype, self.haps[0], self.haps[1], self.pos, self.pos+1,
                self.pos, self.hetp, 1, self.ref_count, self.alt_count, 0,
                self.count, self.tot_count))

def simulate_BNB(mean, sigma, n):
    # sys.stderr.write("%g %g %g\n" % (mean, sigma, n))
    mean_p = np.float64(n) / (n+mean)
    sigma = (1 / sigma)**2
    a = mean_p * (sigma)+1
    b = (1 - mean_p)*sigma
    
    p = beta(a, b)
    #sys.stderr.write("%f %f\n"%(n,p))
    counts = negative_binomial(n, p)
    return counts

def simulate_BB(tot, mean_p, sigma):
    a = mean_p * (1/sigma**2 - 1)
    b = (1-mean_p) * (1/sigma**2 - 1)

    p = beta(a,b)
    counts = binomial(tot,p)
    #sys.stderr.write("%f %f %i\n"%(mean_p,p,counts))
    return counts, (tot-counts)

main()
