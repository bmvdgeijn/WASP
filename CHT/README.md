CHT -  Combined Haplotype Test
======================

This directory contains source code for the Combined Haplotype Test (CHT) portion of the WASP pipeline.

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

The test is described in  
[our paper](http://biorxiv.org/content/early/2014/11/07/011221): van de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci"


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
13. REGION.SNP.REF.HAP.COUNT - number of reads at each linked SNP that match the allele that is linked to the test SNP reference allele [*0;0;0;0;0;0;0;0;0;0;0;0;0;1;0*]
14. REGION.SNP.ALT.HAP.COUNT - number of reads at each linked SNP that match the allele that is linked to the test SNP alternative allele [*0;1;0;0;0;1;0;1;1;0;0;3;1;1;3*]
15. REGION.SNP.OTHER.HAP.COUNT - number of reads matching neither ref nor alt at each linked SNP [*0;0;0;0;0;0;0;0;0;0;0;0;0;0;0*]
16. REGION.READ.COUNT - number of reads in region [*72*]
17. GENOMEWIDE.READ.COUNT - total number of mapped, filtered reads for this individual [*26219310*]


These input files can be generated using the script `extract_haplotype_read_counts.py`.

An example workflow is provided in the file `../exapmle_workflow.sh`.

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
        -b <float> choose a total read depth overdispersion parameter (default=0.00001)
        -o <float> choose an allele specific overdispersion parameter (default=0.00001)
        -m <int>   choose a threshold for minimum number of allele specific reads (default 0)
        -e <float> choose a sequencing error rate in the reads (default 0.005)

An example input file named H3K27ac\_example\_input\_file.txt is included in this repo. This can be used to run the combined haplotype test on H3K27ac data from 10 individuals. The input file points to files that contain read count and haplotype information for dsQTLs, which can be downloaded from [here](http://eqtl.uchicago.edu/histone_mods/haplotype_read_counts/dsQTLs/).


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


Updating total read depths and heterozygous probabilities
-----

Prior to running the CHT, we recommend updating the heterozygote probabilities for the target SNPs using all of the available reads for each individual. This is done using the script ``update_het_probs.py``.

We also recommend correcting read depths based on the fraction of reads that fall within peaks for each sample and the GC content of each region. This is done using the script ``update_total_depth.py``.


Choosing dispersion parameters
------

The CHT includes three parameters that model the dispersion of the read count data. One dispersion parameter is estimated separately for each region when the CHT is run. The other two parameters are estimated across all regions, and can be considered hyperparameters. They are estimated before running the CHT using the scripts `fit_as_coefficients.py` and `fit_bnb_coefficients.py`. 

We suggest running the CHT on permuted data and visualizing the results with a quantile-quantile plot to ensure that the test is working properly. If the permutations do not follow the null distribution, this may indicate that the dispersion parameters have not been accurately estimated. In this case, you may manually set the overdispersion parameters or adjust the p-values according to the permuted distribution. One of the command line options (`-s`) runs the test on permuted genotypes. 


Contact
------

For questions about the combined haplotype test, please contact Graham McVicker 
(gpm@stanford.edu) or Bryce van de Geijn (bmvdgeijn@uchicago.edu).

