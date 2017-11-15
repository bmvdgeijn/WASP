Version 0.2.2 - November 15, 2017
-----------

This release incorporates several modifications/fixes that have been made to
the master branch over the past year.

Changes include:
* check first 2 bytes of input files to decide if they are gzipped rather than relying on .gz filename extension
* fix bug in which unmapped reads in BAM were causing find_intersecting_snps.py to crash
* fix bug in which some cigar codes caused find_intersecting_snps to crash because of typo in snptable
* fix bug related to soft-clipped reads (switch to use pysam 'query_sequence' instead of 'query' attribute)
* added check for pysam version to rmdup_pe and find_intersecting_snps.py
* added check for pytables version to find_intersecting_snps.py
* rescale count totals in CHT so that alpha and beta estimates stay in reasonable range
* allow empty header lines in VCF parsing
* change max SNP identifier length from 16 to 255
* include much more information in output table from combined haplotype test such as SNP info, total numbers of AS reads, etc.
* warn and gracefully handle situation where there are no matching samples on a chromosome


Version 0.2.1 - September 5, 2016
-----------
This minor release corrects an issue where the cht_data module was not included in the v0.2 release. This module is required by the fit_bnb_coefficients.py script.


Version 0.2 - September 3, 2016
-----------

Version 0.2 of WASP is a major update to the code,
especially the mapping code. It fixes several bugs related
to how paired-end reads are handled. For this reason it is
strongly recommended that users switch to this version
of the pipline.

Changes include:
* re-wrote mapping scripts to make simpler and more modular
* re-wrote mapping test scripts and added of many tests
* fixed several mapping pipeline bugs related to paired-end reads
* find_intersecting_snps.py window size no longer required (is now
	unlimited)
* find_intersecting_snps.py can now take HDF5 files as input
* find_intersecting_snps.py can now consider only haplotypes
	present in samples, rather than all possible allelic combinations
	of SNPs overlapping reads.
* added get_as_counts.py script that outputs allele-specific read
	counts at all polymorphic SNPs. 
* snp2h5 now records sample info in output HDF5 files
* improved speed of many CHT pipeline steps
* improved stability of CHT dispersion parameter estimation
* added Snakemake workflows for both mapping and CHT pipelines
* added qqplot.R script to CHT workflow


Version 0.1
-----------
Initial version of WASP
