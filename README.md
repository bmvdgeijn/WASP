# WASP: allele-specific pipeline for unbiased read mapping and molecular QTL discovery

## Introduction

WASP is a suite of tools for unbiased allele-specific read mapping and
discovery of molecular QTLs

WASP is described in
[our paper](http://biorxiv.org/content/early/2014/11/07/011221): van
de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP:
allele-specific software for robust discovery of molecular
quantitative trait loci"

WASP has two parts, which can be used independently of each
other: 

1. Read filtering tools that correct for biases in allele-specific
   mapping. 

2. A Combined Haplotype Test (CHT) that tests for genetic association
   with a molecular trait using counts of mapped and allele-specific
   reads.

The following directories and files are included with WASP.
Each directory contains its own README file:

* [CHT](./CHT) - Code for running the Combined Haplotype Test

* [mapping](./mapping) - Mappability filtering pipeline for correcting allelic mapping biases

* [snp2h5](./snp2h5) - Contains snp2h5 and fasta2h5:  programs for converting
  common SNP and sequence data formats (IMPUTE, VCF and FASTA)
  to an efficient binary format, HDF5.

* [examples](./examples) - Example data files that can be used to try out the
  Combined Haplotype Test.

* [example_mapping_workflow.sh](./examples/example_mapping_workflow.sh) - A script illustrating
  how each step of the Mappability Filtering Pipeline can be run. 

* [example_cht_workflow.sh](./examples/example_cht_workflow.sh) - A script illustrating
  how each step of the Combined Haplotype Test workflow can be run.


## Dependencies

WASP is written in C and python and uses an efficient binary file format
known as [HDF5](https://www.hdfgroup.org/HDF5/).

WASP depends on the following:

* python version 3.x

* [numpy](http://www.numpy.org)

* [scipy](http://www.scipy.org).

* [pysam](https://github.com/pysam-developers/pysam) version 0.8.4 or higher.

* the HDF5 C library version 1.8 or higher

* [PyTables](http://www.pytables.org/) version 3.x


The easiest way to install [HDF5](https://www.hdfgroup.org/HDF5/),
[numpy](http://www.numpy.org), [scipy](http://scipy.org) and
[Pytables](http://www.pytables.org/) is to download and install
[Anaconda](http://continuum.io/downloads). *Installing Anaconda is
highly recommended.* After installing Anaconda, configure [Bioconda] (https://bioconda.github.io/)
and do `conda install pysam`, or download and install 
[pysam](https://github.com/pysam-developers/pysam) directly.


## Installation

1. Download and install [Anaconda](http://continuum.io/downloads),
(or download and install Numpy, Scipy, HDF5, and Pytables separately).

2. Configure Bioconda and install pysam:

		conda config --add channels r
		conda config --add channels bioconda
		conda install pysam
  Alternatively, download and install [pysam](https://github.com/pysam-developers/pysam)
  yourself.

3. Make sure that the HDF5 library is in your library path. For example 
on Linux or OSX you can add the following to your .bashrc or .profile (replace
$CONDA_PREFIX with your Anaconda installation directory if CONDA_PREFIX is not defined):

		export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

4. Clone or download the WASP repository from github:

		# clone the WASP repository
		git clone https://github.com/bmvdgeijn/WASP.git

		# Alternatively download the respository instead:
		wget https://github.com/bmvdgeijn/WASP/archive/master.zip

5. Compile snp2h5 (optional: only needs to be done if you plan to use
snp2h5 or fasta2h5). First modify the snp2h5/Makefile to point to the
Anaconda (or HDF5) installation directory. For example open
[snp2h5/Makefile](snp2h5/Makefile) with a text editor and modify the
HDF_INSTALL variable to point to your Anaconda installation directory:

		HDF_INSTALL = $(CONDA_PREFIX)

	Now compile snp2h5 using make:

		cd WASP/snp2h5
		make


