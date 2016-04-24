#!/bin/bash

snakemake --cluster "qsub -V" --jobs 5
