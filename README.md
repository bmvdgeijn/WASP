CHT -  Combined Haplotype Test
======================

This repository contains source code for running the Combined Haplotype Test. The 
test is described in the supplementary material for [our paper](http://dx.doi.org/10.1126/science.1242429): McVicker G\*, van de Geijn B\*, Degner JF, Cain CE, Banovich NE, Raj A, 
Lewellen N, Myrthil M, Gilad Y, Pritchard JK. "Identification of genetic variants that 
affect histone modifications in human cells" Science. 2013


Input Format
----------

The combined haplotype test script takes a set of input files with
read count and genotype information for each individual. 

Each file is a space-delimited text file that contains a single header line.

The columns in the input file are as follows:

1. CHROM - name of chromosome, e.g. chr1
2. TEST.SNP.POS - position of test SNP on chromosome (first base is numbered 1),  e.g. 963531
3. TEST.SNP.ID - identifier for test SNP, rs2488994
4. TEST.SNP.REF.ALLELE - reference allele, e.g. A
5. TEST.SNP.ALT.ALLELE - alternate allele, e.g. G
6. TEST.SNP.GENOTYPE - sum of genotype probabilties e.g. 1.01 (near 0 = homo ref, near 1 = het, near 2 = homo alt)
7. TEST.SNP.HAPLOTYPE - gives phasing of alleles e.g. 0|1
8. REGION.START - start of region tested for association (where read counts are from) e.g. 963155
9. REGION.END - end of region tested for association e.g. 965155
10. REGION.SNP.POS - positions of linked SNPs in test region e.g. 963199;963240;963531;964027;964043;964062;964088;964159;964216;964218;964219;964357;964433;964525;964996
11. REGION.SNP.HET.PROB - heterozygote probabilities for linked SNPs, e.g. 0.99;0.99;0.99;1.00;0.99;0.99;0.99;0.99;1.00;0.99;0.99;0.99;0.99;0.99;0.99
12. REGION.SNP.LINKAGE.PROB - linkage probabilities for linked SNPs, e.g. 1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00
13. REGION.SNP.REF.HAP.COUNT - number of allele-specific reference-matching reads at each linked SNP: e.g. 0;0;0;0;0;0;0;0;0;0;0;0;0;1;0
14. REGION.SNP.ALT.HAP.COUNT - number of allele-specific alternate-matching reads at each linked SNP: e.g. 0;1;0;0;0;1;0;1;1;0;0;3;1;1;3
15. REGION.SNP.OTHER.HAP.COUNT - number of reads matching neither ref nor alt at each linked SNP, e.g. 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0
16. REGION.READ.COUNT - numer of reads in region e.g. 72
17. GENOMEWIDE.READ.COUNT - total number of mapped, filtered reads for this individual, e.g. 26219310

Note: We have written a set of scripts to generate these files,
however these have many file and source code dependencies that may
make them difficult for other people to use. That said, if you are
interested in obtaining the scripts, please contact us, and we can
make them available.




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



