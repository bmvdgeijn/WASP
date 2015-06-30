snp2h5
======
This directory contains code for *snp2h5*, a program to convert
text files containing polymorphism data into HDF5 files. Currently,
[VCF files](http://faculty.washington.edu/browning/beagle/intro-to-vcf.html) and
[IMPUTE files](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html)
are supported.

HDF5 files are compressed platform-independent binary files. Data in
HDF5 files can be efficiently accessed using libraries written in C,
Python, R and other languages.


## Dependencies
snp2h5 depends on the HDF5 library (version 1.6 or higher). Pytables
or h5py can be used to access the HDF5 files from within python.

The easiest way to reliably install the HDF5 library (and Pytables) is 
to download and install [Anaconda](http://continuum.io/downloads).


## Compiling
snp2h5 is written in C. To compile, first install HDF5 (for example by installing
Anaconda). Make sure that the HDF5 library is in your library path. For example 
on Linux you could add the following to your .bashrc or .profile (replace
$HOME/anaconda/lib with the relevant directory):

    export LD_LIBRARY_PATH=$HOME/anaconda/lib:$LD_LIBRARY_PATH
       
Next modify the Makefile so that HDF_INSTALL points to the
HDF5 installation directory (replace $(HOME)/anaconda with the
appropriate path):

    HDF_INSTALL = $(HOME)/anaconda

Now compile snp2h5 using make:

    make

## Running snp2h5

    snp2h5 OPTIONS INPUT_FILE1 [INPUT_FILE2 ...]

### Input Files:
A separate input file must be provided for each chromosome. The
filename must contain the name of the chromosome. Chromosome names
should match those in the CHROM_FILE, which is provided by the --chrom option.

### Input Options:
* --chrom CHROM_FILE [required]

     Path to chromInfo.txt file (may be gzipped) with list of chromosomes
     for the relevant genome assembly. Each line in file should
     contain tab-separated chromosome name and chromosome length (in
     basepairs). chromInfo.txt files can be downloaded from the UCSC
     genome browser. For example, a chromInfo.txt.gz file for hg19 can
     be downloaded from
     http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

* --format vcf|impute [required]

     Specifies the format of the input files. Currently vcf and impute
     can be specified.


### Output Options:
*  --geno_prob GENO_PROB_OUTPUT_FILE [optional]
    
     Path to HDF5 file to write genotype probabilities to.  This option can
     only be used if the input VCF files provide genotype likelihoods
     (GL in the VCF FORMAT specifier) or input impute files contain
	 genotype probabilities.

*  --haplotype HAPLOTYPE_OUTPUT_FILE [optional]

     Path to HDF5 file to write haplotypes to.  This option can only be
     used if the input VCF files provide genotypes (GT in the VCF
     FORMAT specifier) or the impute files contain haplotypes. The
     genotypes should be phased (| haplotype separator in GT
     field). If the genotypes are unphased (/ haplotype separator
	 in GT field)  a warning is printed.

*  --snp_index SNP_INDEX_OUTPUT_FILE [optional]

     Path to HDF5 file to write SNP index to. The SNP index can
     be used to convert the genomic position of a SNP to its
     corresponding row in the geno_prob, haplotype and snp_tab
     HDF5 files.

*  --snp_tab SNP_TABLE_OUTPUT_FILE [optional]

     Path to HDF5 file to write SNP table to. Each row of SNP
     table contains SNP name (rs_id), position, allele1, allele2.

### Examples:

    # read 1000 genomes VCF files and write haplotype, snp_index
    # and snp_tab to HDF5 files
    snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \
	       --format vcf \
           --haplotype haplotypes.h5 \
           --snp_index snp_index.h5 \
           --snp_tab   snp_tab.h5 \
           1000G/ALL.chr*.vcf.gz


    # read 1000 genomes VCF files that contain genotype likelihoods
    # and write genotype probabilties to HDF5 file
    snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \
	       --format vcf \
           --geno_prob geno_probs.h5 \
           1000G/supporting/genotype_likelihoods/shapeit2/ALL.chr*.gl.vcf.gz

    # read IMPUTE-formatted files and write genotype probabilties,
	# haplotypes, snp_index and snp_tab to HDF5 files
    snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \
	       --format impute \
           --geno_prob geno_probs.h5 \
		   --snp_index snp_index.h5 \
		   --snp_tab snp_tab.h5 \
		   --haplotype haplotypes.h5 \
           data/*.hg19.impute2.gz data/*.hg19.impute2_haps.gz



## Output Files

### snp_tab
An HDF5 file containing one table of SNPs per chromosome.  Each row
provides a description of a single SNP with the following fields:

    name = String[16]
    pos = Int64
    allele1 = String[32]
    allele2 = String[32]

Note that maximum length of alleles is 32bp--long indel alleles are truncated.

### snp_index
An HDF5 file containing one array per chromosome. The array is used to lookup
SNPs by genomic position.  The value at each index in the array is -1 (no
SNP) or the row index of a SNP that can be used to obtain information
from the snp\_tab, geno\_probs and haplotype files.

### geno_probs
An HDF5 file containing one matrix of (float32) genotype probabilities
per chromosome. Each row corresponds to a SNP in the snp_tab.  There
are 3M columns, where M is number of samples.

Columns 1-3 are for the first sample, columns 4-6 are for the second
sample, etc.  The columns for each sample are "homozygous reference
probability", "heterozygous probability", and "homozygous
non-reference probability" and should sum to 1.0.

### haplotypes
An HDF5 file containing a matrix of phased (int8) genotypes. Each row
corresponds to a SNP in the snp_tab.  There are 2M columns, where M is
number of samples. Each column represents a haplotype and contains 0s
(for reference allele) and 1s (for non-reference allele).

Columns 1-2 are for haplotypes 1 and 2 for sample 1, columns 3-4 are
for haplotypes 1 and 2 for the second individuals, etc.



