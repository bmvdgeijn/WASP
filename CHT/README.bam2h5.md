## bam2h5.py

bam2h5.py is a python program for converting BAM
files to read count HDF5 files. Specifically, bam2h5.py counts the number
of reads that match the alternate and reference allele at every SNP in the provided
SNP HDF5 data files. The read counts are stored in specified HDF5 output
files.

Additionally counts of all reads are stored in another track (at the 
left-most chromosome position of the reads).

This program does not perform filtering of reads based on mappability.
It is assumed that the inpute BAM files are filtered appropriately prior to 
calling this script.

Reads that overlap known indels are not included in allele-specific
counts.

## Dependencies
bam2h5 requires the [pysam python library](https://github.com/pysam-developers/pysam)


## Running bam2h5.py

    bam2h5.py OPTIONS BAM_FILE1 [BAM_FILE2 ...]

### BAM Files:

Aligned reads are read from one or more BAM files. The provided
BAM files must be sorted and indexed.

### Input Options:
* --chrom CHROM_TXT_FILE [required]

    Path to chromInfo.txt file (may be gzipped) with list of
	chromosomes for the relevant genome assembly. Each line
	in file should contain tab-separated chromosome name and
	chromosome length (in basepairs). chromInfo.txt files can
	be downloaded from the UCSC genome browser. For example,
	a chromInfo.txt.gz file for hg19 can be downloaded from
	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

* --snp_index SNP_INDEX_H5_FILE [required]

    Path to HDF5 file containing SNP index. The SNP index is
    used to convert the genomic position of a SNP to its
    corresponding row in the haplotype and snp_tab
    HDF5 files.

* --snp_tab SNP_TABLE_H5_FILE [required]

    Path to HDF5 file to read SNP information from. Each row of SNP
    table contains SNP name (rs_id), position, allele1, allele2.

* --haplotype HAPLOTYPE_H5_FILE [optional]

    Path to HDF5 file to read phased haplotypes from.
    If supplied, when read overlaps multiple SNPs counts are randomly
    assigned to ONE of the overlapping HETEROZYGOUS SNPs; if not supplied 
    counts are randomly assigned to ONE of overlapping SNPs (regardless of 
    their genotype).

* --samples SAMPLES_TXT_FILE [optional]

    Path to text file containing a list of individual identifiers. The
    ordering of individuals must be consistent with the haplotype
    file. The samples file is assumed to have one identifier per line
    in the first column (other columns are ignored).

* --individual INDIVIDUAL [optional]

    Identifier for individual, used to determine which
    SNPs are heterozygous. Must be provided
    if --haplotype argument is provided and must match one of the
    individuals in the file provided with --samples argument.

### Output Options:
* --data_type uint8|uint16

    Data type of stored counts; uint8 takes up less disk
    space but has a maximum value of 255 (default=uint8).

* --ref_as_counts REF_AS_COUNT_H5_FILE [required]
	 
     Path to HDF5 file to write counts of reads that match reference allele.
     Allele-specific counts are stored at the position of the SNP.

* --alt_as_counts ALT_AS_COUNT_H5_FILE [required]

    Path to HDF5 file to write counts of reads that match alternate allele.
    Allele-specific counts are stored at the position of the SNP.

* --other_as_counts OTHER_AS_COUNT_H5_FILE [required]

    Path to HDF5 file to write counts of reads that match neither reference
    nor alternate allele. Allele-specific counts are stored at the position
    of the SNP.

* --read_counts READ_COUNT_H5_FILE [required]

    Path to HDF5 file to write counts of all reads, regardless of whether
    they overlap a SNP. Read counts are stored at the left-most position
    of the mapped read.


### Examples:

    # write HDF5 read count files for individual 18505
	INDIVIDUAL=18505
	
    ./bam2h5/bam2h5.py --chrom chromInfo.hg19.txt \
	      --snp_index snp_index.h5 \
	      --snp_tab snp_tab.h5 \
	      --haplotype haps.h5 \
	      --samples H3K27ac/samples.txt \
	      --individual $INDIVIDUAL \
	      --ref_as_counts ref_as_counts.$INDIVIDUAL.h5 \
	      --alt_as_counts alt_as_counts.$INDIVIDUAL.h5 \
	      --other_as_counts other_as_counts.$INDIVIDUAL.h5 \
	      --read_counts read_counts.$INDIVIDUAL.h5 \
	      H3K27ac/$INDIVIDUAL.chr*.keep.rmdup.bam

