#!/bin/bash
#
# This script takes VCF files and generates files that can be used
# by the find_intersecting_snps.py script. 
# The script takes an INPUT directory and an OUTPUT directory.
#
# The INPUT directory is expected to contain files ending with .vcf.gz
# and to contain the name of the chromosome (like chr22 or chr1).
#
# OUTPUT files are named $CHR.snps.txt.gz (where $CHR is the name of the
# chromosome). OUTPUT files contain <position> <allele1> <allele2> on
# each line.
#

INPUT_DIR=$1
OUTPUT_DIR=$2

if [ ! $INPUT_DIR ]; then
    echo "usage: extract_vcf_snps.sh <input_dir> <output_dir>" >&2
    exit 2
fi

if [ ! $OUTPUT_DIR ]; then
    echo "usage: extract_vcf_snps.sh <input_dir> <output_dir>" >&2
    exit 2
fi

mkdir -p $OUTPUT_DIR

for FILE in $INPUT_DIR/*vcf.gz; do
     echo $FILE >&2
     CHR=`echo $FILE | sed -n 's/^.*\(chr[0-9A-Z]*\).*.vcf.gz$/\1/p'`
     echo $CHR >&2
     OUTPUT_FILE=$OUTPUT_DIR/$CHR.snps.txt.gz
     gunzip -c $FILE | egrep -v "^#" | awk '{print $2,$4,$5}' | gzip > $OUTPUT_FILE
done

