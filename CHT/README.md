# CHT -  Combined Haplotype Test

This directory contains the Combined Haplotype Test (CHT) portion of
the WASP pipeline.

## Introduction

The CHT is designed to test for genetic associations with
quantitatitve traits that can be measured with next-generation
sequencing reads. For example experiments such as ChIP-seq or RNA-seq
that generate read counts can be used. The test uses both read depth
from a test region and allelic imbalance of reads that overlap phased
heterozygous SNPs in the same test region. The test SNP whose
genotypes are tested for association can be within or nearby the test
region. (The SNP should not be too far away, because phasing
information between SNPs is currently assumed to be correct.)

The test accounts for overdispersion in the data (excess variation
that cannot be attributed to binomial sampling or genotype) with three
dispersion parameters. Two of these parameters are globally estimated
for each individual (i.e. fixed across test regions), while one is
estimated for each region.

The test is described in
[our paper](http://biorxiv.org/content/early/2014/11/07/011221): van
de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP:
allele-specific software for robust discovery of molecular
quantitative trait loci"


## Input Format

The combined haplotype test script takes a set of input files with
read count and genotype information for each individual.

Each file is a space-delimited text file that contains a single header
line.

The columns in the input file are as follows (example values are given
in [*brackets*]):

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

These input files can be generated using the following workflow or
created by the user.


## Obtaining phased genotype data

Step 3 of the CHT workflow requires *phased*
genotype data.  The [snp2h5](../snp2h5/README.md) program provided
with WASP can read phased genotypes from IMPUTE2 or VCF
formatted files.

Several methods are available for genotype phasing:

* [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html) can
  perform both imputation and genotype phasing. (Phasing information
  is written to a file specified with with -o option). Output
  files written by IMPUTE2 can be read directly by snp2h5.

* [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
  is a fast method for genotype phasing. Output files from SHAPEIT can
  be converted to IMPUTE2's .hap format or to VCF format using the
  `shapeit -convert`  command described in the SHAPEIT documentation.
  
* [BEAGLE](https://faculty.washington.edu/browning/beagle/b3.html) can
  be used for both imputation and genotype phasing. BEAGLE output
  files can be converted to VCF format using the
  [beagle2vcf](http://faculty.washington.edu/browning/beagle_utilities/utilities.html#beagle2vcf)
  utility program provided with BEAGLE.

If you are using samples that are part of the 1000 Genomes project,
genotypes that have been phased using SHAPEIT can be downloaded from
the [1000 Genomes website](http://www.1000genomes.org/data#DataAccess).


##  Workflow

An example workflow is provided in the file
`../example_workflow.sh`. This workflow uses example data under the
directory `../example_data`.

Some of the input files that we used for our paper can be downloaded from 
[here](http://eqtl.uchicago.edu/histone_mods/haplotype_read_counts/). 

The following steps can be used to generate input files and run the
Combined Haplotype Test. The examples given below use the example
data files that are provided in the `../example_data` directory.


### Step 1

Map reads to the genome, and filter them to correct for mapping
bias. This procedure is described by the README in the mapping
directory.

### Step 2

Convert SNP data to HDF5 format using the program snp2h5.  HDF5 files
are an efficient binary data format used by WASP scripts. The snp2h5
program can take VCF or IMPUTE2 files as input.

	./snp2h5/snp2h5 --chrom example_data/chromInfo.hg19.txt \
	      --format impute \
	      --geno_prob example_data/geno_probs.h5 \
	      --snp_index example_data/snp_index.h5 \
	      --snp_tab example_data/snp_tab.h5 \
	      --haplotype example_data/haps.h5 \
	      example_data/genotypes/chr*.hg19.impute2.gz \
	      example_data/genotypes/chr*.hg19.impute2_haps.gz


### Step 3

Convert FASTA files to HDF5 format. Note the HDF5 sequence files are
only used for GC content correction part of CHT. This step can be
ommitted if GC-content correction is not used.

	./snp2h5/fasta2h5 --chrom example_data/chromInfo.hg19.txt \
		--seq example_data/seq.h5 \
		/data/external_public/reference_genomes/hg19/chr*.fa.gz

### Step 4

Extract read counts from BAM files (from Step 1) and write them to
HDF5 files using the [bam2h5.py](README.bam2h5.md) program.  This must
be done for each sample/individual in the dataset.

    python CHT/bam2h5.py --chrom example_data/chromInfo.hg19.txt \
	      --snp_index example_data/snp_index.h5 \
	      --snp_tab example_data/snp_tab.h5 \
	      --haplotype example_data/haps.h5 \
	      --samples $ALL_SAMPLES_FILE \
	      --individual $INDIVIDUAL \
	      --ref_as_counts example_data/H3K27ac/ref_as_counts.$INDIVIDUAL.h5 \
	      --alt_as_counts example_data/H3K27ac/alt_as_counts.$INDIVIDUAL.h5 \
	      --other_as_counts example_data/H3K27ac/other_as_counts.$INDIVIDUAL.h5 \
	      --read_counts example_data/H3K27ac/read_counts.$INDIVIDUAL.h5 \
	      example_data/H3K27ac/$INDIVIDUAL.chr*.keep.rmdup.bam

### Step 5

Create a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-like input
file that defines both the "test SNP" and the "target regions" that
are to be tested for association.

For example, if the goal is to identify histone-mark QTLs, the target
regions should be ChIP-seq peaks, and the test SNPs should be SNPs
that are near-to or within the ChIP-seq peaks.

If the goal is to identify eQTLs, the target regions should be the
exons of genes, and the test SNPs could be SNPs within a specified
distance of the TSS.

Currently these input files must be generated by the user based on the
regions and SNPs that they want to use as input to the test.  We hope
to soon provide scripts that can be used to generate these input files
for different types of datasets.

We provide an example of a a test SNP / target region input file here:

	example_data/H3K27ac/chr22.peaks.txt.gz

The whitespace-delimited columns of the input file are as follows (example
values are given in [*brackets*])::

1. The chromosome of the test-SNP [*chr22*]
2. The position of the SNP. Unlike BED files, first base of chromosome
is 1. [*16050408*]
3. The end position of the SNP. *Currently ignored*. [*16050408*]
4. The reference allele of the SNP [*T*]
5. The non-reference allele of the SNP [ *C*]
6. The strand of the SNP.  *Currently ignored and assumed to be '+' *. [*+*]
7. The name of the test. *Currently ignored*. [*chr22. 16050408*]
8. A single start coordinate or a ';'-delimited list of start coordinates for
the target region. Multiple starts can be given if the target region is discontiguous,
such as the exons of a gene. [*16049408*]
9. A single end coordinate or a ';'-delimited list of end coordinate
   for the target region. [*16051408*]



### Step 6

Create a CHT input file for each individual using the
extract_haplotype_read_counts.py script. This script reads the test
SNPs and target regions from the input file created in Step 5.

    python CHT/extract_haplotype_read_counts.py \
       --chrom example_data/chromInfo.hg19.txt \
       --snp_index example_data/snp_index.h5 \
       --snp_tab example_data/snp_tab.h5 \
       --geno_prob example_data/geno_probs.h5 \
       --haplotype example_data/haps.h5 \
       --samples $ALL_SAMPLES_FILE \
       --individual $INDIVIDUAL \
       --ref_as_counts example_data/H3K27ac/ref_as_counts.$INDIVIDUAL.h5 \
       --alt_as_counts example_data/H3K27ac/alt_as_counts.$INDIVIDUAL.h5 \
       --other_as_counts example_data/H3K27ac/other_as_counts.$INDIVIDUAL.h5 \
       --read_counts example_data/H3K27ac/read_counts.$INDIVIDUAL.h5 \
       example_data/H3K27ac/chr22.peaks.txt.gz \
       | gzip > example_data/H3K27ac/haplotype_read_counts.$INDIVIDUAL.txt.gz


### Step 7

Adjust read counts in CHT input files by modeling relationship between read
depth and GC content & peakiness in each sample.

	IN_FILE=example_data/H3K27ac/input_files.txt
	OUT_FILE=example_data/H3K27ac/output_files.txt
	ls example_data/H3K27ac/haplotype_read_counts* | grep -v adjusted > $IN_FILE
	cat $IN_FILE | sed 's/.txt/.adjusted.txt/' >  $OUT_FILE

	python CHT/update_total_depth.py --seq example_data/seq.h5 $IN_FILE $OUT_FILE


### Step 8

Adjust heterozygote probabilities in CHT input files to account for possible
genotyping errors. Total counts of reference and alternative alleles
are used to adjust the probability. In this example we just provide
the same H3K27ac read counts, however you could also use read counts
combined across many different experiments or (perhaps ideally) from
DNA sequencing.

This must be done for each individual in the dataset.

    python CHT/update_het_probs.py \
	   --ref_as_counts example_data/H3K27ac/ref_as_counts.$INDIVIDUAL.h5  \
	   --alt_as_counts example_data/H3K27ac/alt_as_counts.$INDIVIDUAL.h5 \
	   example_data/H3K27ac/haplotype_read_counts.$INDIVIDUAL.adjusted.txt.gz
	   example_data/H3K27ac/haplotype_read_counts.$INDIVIDUAL.adjusted.hetp.txt.gz


### Step 9

Estimate overdispersion parameters for the allele-specific test (beta
binomial):

	python CHT/fit_as_coefficients.py \
		example_data/H3K27ac/cht_input_files.txt \
		example_data/H3K27ac/cht_as_coef.txt

and for the for association test (beta-negative binomial):

	python CHT/fit_bnb_coefficients.py \
		--min_counts 50 \
		--min_as_counts 10 \
		example_data/H3K27ac/cht_input_files.txt \
		example_data/H3K27ac/cht_bnb_coef.txt

### Step 10

Run the combined haplotype test:

	python CHT/combined_test.py --min_as_counts 10 \
		--bnb_disp example_data/H3K27ac/cht_bnb_coef.txt \
		--as_disp example_data/H3K27ac/cht_as_coef.txt \
		example_data/H3K27ac/cht_input_files.txt \
		example_data/H3K27ac/cht_results.txt


## Output

The combined haplotype test outputs a tab delimited text file with the
following columns:

1. chromosome - location of tested SNP 
2. position - position of tested SNP
3. log likelihood under null model
4. log likelihood under alternative model
5. chi-squared value - likelihood ratio test statistic
6. p-value - p-value from likelihood ratio test
7. best alpha parameter - maximum likelihood estimate of alpha parameter (expression level of reference allele)
8. best beta parameter - maximum likelihood estimate of beta parameter (expression level of non-reference allele)
9. overdispersion parameter - maximum  likelihood estimate of phi parameter (per-region beta-negative-binomial dispersion parameter)
10. number of AS reads - number of allele specific reads in tested regio,n summed across individuals
11. total reads - number of mapped reads in tested region, summed across individuals



## Updating total read depths and heterozygous probabilities

Prior to running the CHT, we recommend updating the heterozygote
probabilities for the target SNPs using all of the available reads for
each individual. This is done using the script
``update_het_probs.py`` as described in Step 8 above.

We also recommend correcting read depths based on the fraction of
reads that fall within peaks for each sample and the GC content of
each region. This is done using the script ``update_total_depth.py``
as described in Step 7 above.


## Choosing dispersion parameters

The CHT includes three parameters that model the dispersion of the
read count data. One dispersion parameter is estimated separately for
each region when the CHT is run. The other two parameters are
estimated across all regions, and can be considered
hyperparameters. They are estimated before running the CHT using the
scripts `fit_as_coefficients.py` and `fit_bnb_coefficients.py` as
described in Step 9 above.

We suggest running the CHT on permuted data and visualizing the
results with a quantile-quantile plot to ensure that the test is
working properly. If the permutations do not follow the null
distribution, this may indicate that the dispersion parameters have
not been accurately estimated. In this case, you may manually set the
overdispersion parameters or adjust the p-values according to the
permuted distribution. One of the command line options (`-s`) runs the
test on permuted genotypes.


## Contact

For questions about the combined haplotype test, please contact Graham McVicker 
(gpm@stanford.edu) or Bryce van de Geijn (bmvdgeijn@uchicago.edu).

