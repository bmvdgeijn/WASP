WASP: allele-specific pipeline for unbiased read mapping and molecular QTL discovery
----

WASP is a suite of tools for unbiased allele-specific read mapping and
discovery of molecular QTLs

WASP has two parts, which can be used independently of each
other: 

1. Read filtering tools that correct for biases in allele-specific
   mapping. 

2. A Combined Haplotype Test (CHT) that tests for genetic association
   with a molecular trait using counts of mapped and allele-specific
   reads.

The following directories and files are included with WASP.
Each directory contains its own README file:

* CHT - Code for running the Combined Haplotype Test

* mapping -Tools for correcting mapping biases

* snp2h5 - A program for converting common SNP data formats
  (IMPUTE and VCF) to an efficient binary format, HDF5.

* example_data - Example data files that can be used to try out the
  Combined Haplotype Test.

* example_workflow.sh - A script illustrating how each step of the
Combined Haplotype Test workflow can be run.

