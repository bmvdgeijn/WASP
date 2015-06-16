This directory contains a set of example input and output datafiles
for the combined haplotype test pipeline. The example data is H3K27ac
ChIP-seq data and genotypes for 10 lymphoblastoid cell lines that are
part of the hapmap and 1000 genomes projects. To keep the datafiles
small only data from chromosome 22 are provided.

Several of the datafiles used by the pipeline are in HDF5 format. HDF5
is a compressed binary file format that can be used efficiently within
python and other programming languages. We provide programs to convert
file formats such as IMPUTE2, VCF and FASTA to HDF5.

Descriptions of the example data files are provided below.


## Genotype files

### genotypes/YRI_samples.txt

List of sample identifiers for genotype files. This example file
contains sample identifiers for 60 Yoruba lymphoblastoid cell lines.
The ordering of samples in this file corresponds to the ordering of
genotypes and haplotypes in other genotype files. This file is
provided by the user.

###  genotypes/chr22.hg19.impute2.gz

SNP information and genotypes in  IMPUTE2 output format.
Data in this format can be converted to HDF5 format by the snp2h5 program.
If desired, VCF files can be used instead of IMPUTE2 files. This
file is provided by the user.

### genotypes/chr22.hg19.impute2_haps.gz

Haplotype phasing information for SNPs in IMPUTE2 format.  Data in
this format can be convered to HDF5 format by the snp2h5 prgram.  If
desired, VCF files can be used instead of IMPUTE2 files. This file
is provided by the user.

### snp_tab.h5

HDF5 file containing table of SNP information (identifier, chromosome,
position, alleles, identifier). This file is created by the snp2h5
program.

### geno_probs.h5

HDF5 file containing genotype probabilities for all SNPS and
individuals. This file is createrd by the snp2h5 program.

### haps.h5

HDF5 file containing haplotype phasing and genotype calls for all SNPs
and individuals. Thi file is created by the snp2h5 program.

### snp_index.h5

HDF5 file which maps chromosomal positions to SNP row indices for the
geno_probs.h5, haps.h5, and snp_tab.h5 files. This file is created by
the snp2h5 file.


## Genome files

### ./seq.h5

HDF5 file containing the DNA sequence for the whole genome.  This file
is not provided in the repository because of its large size.  This
file can be created from fasta files, using the provided fasta2h5
program.

### ./chromInfo.hg19.txt

Text file containing names and lengths of all chromosomes in the
assembly. chromInfo.txt files can be downloaded from the UCSC genome
browser. For example, a chromInfo.txt.gz file for hg19 can be
downloaded from
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

## H3K27ac data files

### H3K27ac/samples.txt

List of samples for which H3K27ac ChIP-seq has been performed.  This
file should be provided by the user.

### H3K27ac/\*.chr22.keep.rmdup.bam

BAM files containing mapped reads for each individual. These files
should be sorted and indexed. They should also be filtered for
duplicates and allele-specific biases using the WASP mapping pipeline.

### H3K27ac/ref_as_counts.\*.h5

HDF5 file containing counts of allele-specific reads that match the
reference allele at every position in the genome. Allele-specific
counts are stored at the position of the SNP. This file is created by
the bam2h5.py script.

### H3K27ac/alt_as_counts.\*.h5

HDF5 file containing counts of allele-specific reads that match the
alternate allele at every position in the genome. Allele-specific
counts are stored at the position of the SNP. This file is created by
the bam2h5.py script.

### H3K27ac/other_as_counts.\*.h5

HDF5 file containing counts of allele-specific reads that do not match
the alternative or reference allele. Allele-specific counts are stored
at the position of the SNP. This file is created by the bam2h5.py
script.

### H3K27ac/read_counts.\*.h5

HDF5 containing counts of all mapped reads, regardless of whether they
overlap a SNP. Read counts are stored at the left-most position of the
mapped read. This file is created by the bam2h5.py script.


## Combined Haplotype Test files

### H3K27ac/chr22.peaks.txt.gz

BED-like file containing coordinates of target region(s) and test SNP
for the Combined Haplotype Test. The columns in the file are as follows:

1. name of chromosome [*chr22*]
2. test SNP start position  (uses BED convention that first base numbered 0)  [*16050408*]
3. test SNP end position (uses BED convention of half-open coords) [*16050409*]
4. test SNP reference allele [*T*]
5. test SNP alternate allele [*C*]
6. strand (always +)  [*+*] 
7. name of target region [*chr22.16050409*]
8. target region start(s), ';'-delimited list of start coordinates for
target region(s) with BED-style numbering [*16049408*]
9. target region end(s), ';'-delimited list of end coordinates for
target region(s) with BED-style numbering [*16051408*]

This file is provided by the user.


### H3K27ac/haplotype_read_counts.\*.txt.gz

Combined Haplotype Test input files, which contain information about
read counts, allele-specific reads, genotype probabilities etc. One
file per individual. These files are created by the
extract_haplotype_tead_counts.py script.

The files are tab (or space) delimited  and contain the following
columns  (example values are given in [*brackets*]):

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


### H3K27ac/haplotype_read_counts.*.adjusted.txt.gz

Combined Haplotype Test input files with adjusted genome-wide read
counts values. The values are adjusted to account for GC content
biases and the fraction of reads from peaks for each sample. The
adjusted total read counts are used to calculate an expected number of
reads for each genotype.

The read count adjustment is performed by the update_total_depth.py
script.


### H3K27ac/haplotype_read_counts.*.adjusted.hetp.txt.gz

Combined Haplotype Test input files with adjusted genome-wide read
counts and adjusted heterozygote probabilities. The heterozygote
probabilities are adjusted using the script update_het_probs.py.


### H3K27ac/cht_as_coef.txt

File containing estimates for allele-specific (Beta-Binomial)
dispersion parameters. The file contains one parameter estimate
per individual, each on its own line. This file is generated by the
fit_as_coefficients.py script.

### H3K27ac/cht_bnb_coef.txt

File containing estimates for genome-wide read depth (Beta Negative Binomial)
dispersion parameters. One parameter estimate for each individual,
each written on its own line. This file is generated by the
fit_bnb_coefficients.py script.

### H3K27ac/cht_results.txt

A tab-delimited text file containing results from running the combined
haplotype test. The columns in the file are:

1.  SNP chromosome
2.  SNP position
3.  Log Likelihood under null model
4.  Log likelihood under alternative model
5.  Chi-square test statistic (one degree of freedom)
6.  P-value
7.  ML estimate of alpha parameter (expected reads from reference chromosome)
8.  ML estimate of beta parameter (expected reads from alternate chromosome)
9.  ML estimate of BNB dispersion parameter for region
10.  total count of allele-specific reads for this target region
11. total count of reads for this target region

This file is generated by the script combined_test.py
