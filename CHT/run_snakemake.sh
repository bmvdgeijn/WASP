#!/bin/bash

# Run snakemake on an SGE cluster (that uses qsub command for submission)
# Run at most 30 jobs at a time
# To run on an LSF cluster, change "qsub -V" to "bsub"
snakemake --cluster "qsub -V" --jobs 30
