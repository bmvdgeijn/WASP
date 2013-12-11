CHT -  Combined Haplotype Test
======================

This repository contains source code for running the Combined Haplotype Test. 

Introduction
----------

The test is designed to test for genetic associations with data from ChIP-seq or 
other experiments that generate read counts. The test uses both read depth from a 
test region and allelic imbalance of reads that overlap phased heterozygous SNPs in 
the same test region. The SNP that is tested for association can be within or 
nearby the test region. The SNP should not be too far away, because it is assumed 
to be in LD with the other SNPs (i.e. the phasing information is assumed to be 
correct). 

The test accounts for overdispersion in the data (excess variation that cannot be 
attributed to binomial sampling or genotype) with three dispersion parameters. Two 
of these parameters are global (i.e. fixed across test regions), while one is 
estimated for each region.

The test is described in detail in the the supplementary material for 
[our paper](http://dx.doi.org/10.1126/science.1242429): McVicker G\*, van de Geijn B\*, 
Degner JF, Cain CE, Banovich NE, Raj A, Lewellen N, Myrthil M, Gilad Y, Pritchard JK. 
"Identification of genetic variants that affect histone modifications in human cells" 
Science. 2013

Some components of the test that are described in the paper are not included in the 
current version of the test. This version does not: (1) update the heterozygote probabilities
based on observed reads across all datatypes; or (2) correct variation due to GC content and 
fraction of reads in peak regions in each individual. These parts of the test were removed to 
simplify the code and to remove several dependencies. In principle these corrections could be 
applied to the input files prior to running the test. For example a script could update the 
heterozygote probabilities and read count totals using information about read counts and GC 
content etc.


Input Format
----------

The combined haplotype test script takes a set of input files with
read count and genotype information for each individual. 

Each file is a space-delimited text file that contains a single header line.

The columns in the input file are as follows (example values are given in [*brackets*]):

1. CHROM - name of chromosome [*chr1*]
2. TEST.SNP.POS - position of test SNP on chromosome (first base is numbered 1) [*963531*]
3. TEST.SNP.ID - identifier for test SNP [*rs2488994*]
4. TEST.SNP.REF.ALLELE - reference allele [*A*]
5. TEST.SNP.ALT.ALLELE - alternate allele [*G*]
6. TEST.SNP.GENOTYPE - sum of genotype probabilties (near 0 = homo ref, near 1 = het, near 2 = homo alt) [*1.01*] 
7. TEST.SNP.HAPLOTYPE - gives phasing of alleles [*0|1*]
8. REGION.START - start of region tested for association (where read counts are from) [*963155*]
9. REGION.END - end of region tested for association [*965155*]
10. REGION.SNP.POS - positions of linked SNPs in test region [*963199;963240;963531;964027;964043;964062;964088;964159;964216;964218;964219;964357;964433;964525;964996*]
11. REGION.SNP.HET.PROB - heterozygote probabilities for linked SNPs [*0.99;0.99;0.99;1.00;0.99;0.99;0.99;0.99;1.00;0.99;0.99;0.99;0.99;0.99;0.99*]
12. REGION.SNP.LINKAGE.PROB - linkage probabilities for linked SNPs [*1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00*]
13. REGION.SNP.REF.HAP.COUNT - number of allele-specific reference-matching reads at each linked SNP [*0;0;0;0;0;0;0;0;0;0;0;0;0;1;0*]
14. REGION.SNP.ALT.HAP.COUNT - number of allele-specific alternate-matching reads at each linked SNP [*0;1;0;0;0;1;0;1;1;0;0;3;1;1;3*]
15. REGION.SNP.OTHER.HAP.COUNT - number of reads matching neither ref nor alt at each linked SNP [*0;0;0;0;0;0;0;0;0;0;0;0;0;0;0*]
16. REGION.READ.COUNT - number of reads in region [*72*]
17. GENOMEWIDE.READ.COUNT - total number of mapped, filtered reads for this individual [*26219310*]

Note: We have written a set of scripts to generate these files,
however these have many file and source code dependencies that may
make them difficult for other people to use. That said, if you are
interested in obtaining the scripts, please contact us, and we can
make them available.

Some of the input files that we used for our paper can be downloaded from 
[here](http://eqtl.uchicago.edu/histone_mods/haplotype_read_counts/). 


Running the combined test
---------------------

    python combined_test.py [options] input_files.txt output_file.txt
        input_files.txt - A file containing a list of the input files to be tested
        output_file.txt - File to write output to

    Options:
        -a perform the allele specific portion of the test only
        -d perform the total read depth portion of the test only
        -s randomly permute the genotypes and allele specific counts
        -d <float> choose a total read depth overdispersion parameter (default=0.00001)
        -o <float> choose an allele specific overdispersion parameter (default=0.00001)
        -m <int>   choose a threshold for minimum number of allele specific reads (default 0)
        -e <float> choose a sequencing error rate in the reads (default 0.005)


Output
------

The script writes a tab delimited text file with the following columns:

1. chromosome - location of tested SNP
2. position - position of tested SNP
3. chi-square value - likelihood ratio test statistic
4. best alpha parameter - maximum likelihood estimate of alpha parameter
5. best beta parameter - maximum likelihood estimate of beta parameter
6. overdispersion parameter - maximum  likelihood estimate of per-region dispersion parameter
7. number of AS reads - number  of of allele specific reads in tested region



Choosing dispersion parameters
------

One of the dispersion parameters is estimated for each region using maximum likelhood. 
The other two parameters are fixed across regions and must be specified on the command 
line. In principle it should be possible to write code to estimate these parameters from 
the data, however we have not had time to implement this. In practice we run the code using 
several combinations of dispersion parameters, and choose the values that appear to give 
well-calibrated (uniformly distributed) p-values when the test is run on permuted genotypes
(i.e. p-values that follow the expectation line in a QQ-plot). One of the command line options
(-s) runs the test on permuted genotypes. 


Contact
------

For questions about the combined haplotype test, please contact Graham McVicker 
(mcvicker@jimmy.harvard.edu) or Bryce van de Geijn (bmvdgeijn@uchicago.edu).

