Typical workflow for converting data to HDF5 format

# convert chromInfo.txt to an HDF5 file containing table
# of chromosome information
python convert_chrom.py  --chrom_tab data/chromosomes.h5  chromInfo.txt.gz

# parse VCF files and create two HDF5 files which contain:
#   [1] a SNP table for each chromosome. Each row has information for
#        each SNP including name, position, alleles
#   [2] a snp_index vector for each chromosome that allows SNPs to be looked up
#        by genome position
python convert_snps.py \
    --chrom_tab data/chromosomes.h5 \
    --snp_tab data/snp_tab.h5  \
    --snp_index data/snp_index.h5 \
    1000genomes/ALL.chr*genotypes.vcf.gz

# parse VCF file containing genotype likelihoods and create an HDF5
# file with genotype probabilities for every SNP. Also write a file
# containing identifiers of genotyped samples.
python convert_genos.py \
     --chrom_tab data/chromosomes.h5 \
     --snp_tab data/snp_tab.h5 \
     --genos data/geno_probs.h5 \
     --samples data/samples.txt \
     1000genomes/supporting/genotype_likelihoods/shapeit2/ALL.chr*.vcf.gz

# parse VCF file with phased genotypes and create an HDF5
# file with haplotype information for every SNP. Also write a file
# containing identifiers of genotyped samples.
python convert_haps.py \
    --chrom_tab data/chromosomes.h5 \
    --snp_tab data/snp_tab.h5 \
    --haps data/haplotypes.h5 \
    --samples data/samples.txt \
    1000genomes/shapeit2/ALL.chr*.vcf.gz
 


-------------------------

We store SNP data in HDF5 files, which can be read efficiently from within
python, C, or R.

SNP input data in multiple formats can be converted into HDF5.  For example,
we provide ways to convert IMPUTE2 and VCF files into HDF5 format.

The SNP data are stored in a directory of the users choice. The HDF5 files are
formatted as follows.

samples.txt
---------
List of sample names, one sample per row.

chr.h5
-----
Table of chromosomes for a genome assembly (e.g. hg19 )

    idnum = tables.Int32Col()
    name = tables.StringCol(32)
    length = tables.Int32Col()
    is_sex = tables.BoolCol(dflt=False)
    is_auto = tables.BoolCol(dflt=True)
    is_rand = tables.BoolCol(dflt=False)
    is_hap = tables.BoolCol(dflt=False)
    is_mito = tables.BoolCol(dflt=False)
    is_y = tables.BoolCol(dflt=False)
    is_x = tables.BoolCol(dflt=False)


snp_tab.h5
---------
One table of SNPs per chromosome containing descriptive information
including, name, chromosome position, and alleles.

    name = tables.StringCol(16)
    pos = tables.Int32Col()
    allele1 = tables.StringCol(32)
    allele2 = tables.StringCol(32)

Note that maximum length of alleles is 32bp--long indel alleles are truncated.


snp_index.h5
----------
One array per chromosome that is used as a positional lookup for SNPs.
Each index in the array gives a -1 (no SNP) or the index of a SNP in the
snp_tab.


geno_probs.h5
------------
One matrix per chromosome, containing genotype probabilities.  Each row corresponds to a SNP in the snp_tab.  There are 3M columns, where M is number of samples. The samples must be ordered in the same way as in the samples.txt file.

Columns 1-3 are for the first individual, columns 4-6 are for the second individual, etc.  The columns for each individual are "homozygous reference probability", "heterozygous probability", and "homozygous non-reference probability" and should sum to 1.0.


haplotypes.h5
------------
One matrix per chromosome that provides haplotype phasing. Each row corresponds to a SNP in the snp_tab.  There are 2M columns, where M is number of samples. The samples must be ordered in the same way as in the samples.txt file.

Columns 1-2 are for the first individual, columns 3-4 are for the second individual, etc. 


