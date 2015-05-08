# WASP: allele-specific pipeline for unbiased read mapping and molecular QTL discovery

## Introduction

WASP is a suite of tools for unbiased allele-specific read mapping and
discovery of molecular QTLs

WASP is described in
[our paper](http://biorxiv.org/content/early/2014/11/07/011221): van de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci"

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

* snp2h5 - Contains snp2h5 and fasta2h5:  programs for converting
  common SNP and sequence data formats (IMPUTE, VCF and FASTA)
  to an efficient binary format, HDF5.

* example_data - Example data files that can be used to try out the
  Combined Haplotype Test.

* example_workflow.sh - A script illustrating how each step of the
  Combined Haplotype Test workflow can be run.


## Dependencies

WASP is written in C and python and depends on both [numpy](http://www.numpy.org) and
[scipy](http://www.scipy.org).  The code also depends on
[argparse](https://code.google.com/p/argparse/), which is included by default in newer versions of python (>= 2.7).

Some scripts depend on the [pysam python library](https://github.com/pysam-developers/pysam)

The combined haplotype test uses
[HDF5](https://www.hdfgroup.org/HDF5/), an efficient
compressed binary file format.  For this reason, the Combined
Haplotype Test requires the HDF5 library
(version 1.6 or higher) and [PyTables](http://www.pytables.org/).

The easiest way to install HDF5, numpy, scipy and
Pytables is to download and install
[Anaconda](http://continuum.io/downloads).

*Installing Anaconda is highly recommended.*  After installing
Anaconda, the only other dependency that must be downloaded
and installed is [pysam](https://github.com/pysam-developers/pysam).


## Installation

1. Download and install [Anaconda](http://continuum.io/downloads),
(or download and install Numpy, Scipy, HDF5, and Pytables separately).

2. Download and install [pysam](https://github.com/pysam-developers/pysam)

3. Make sure that the HDF5 library is in your library path. For example 
on Linux or OSX you can add the following to your .bashrc or .profile (replace
$HOME/anaconda with your Anaconda installation directory):

    export LD_LIBRARY_PATH=$HOME/anaconda/lib:$LD_LIBRARY_PATH

4. Clone the WASP repository from github:

    git clone https://github.com/bmvdgeijn/WASP.git

(Alternatively download the respository instead):
	wget https://github.com/bmvdgeijn/WASP/archive/master.zip

5. Compile snp2h5 (optional: only needs to be done if you plan to use
snp2h5 or fasta2h5). First modify the snp2h5/Makefile to point to the
Anaconda (or HDF5) installation directory. For example:
    HDF_INSTALL = $(HOME)/anaconda

Now compile snp2h5 using make:
    cd snp2h5
    make


