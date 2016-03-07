
# Set these environment vars to point to
# your local installation of WASP
WASP=$HOME/proj/WASP
DATA_DIR=$WASP/example_data
SNP_DIR=$DATA_DIR/genotypes/snps

# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2
INDEX=/data/external_public/reference_genomes/hg19/hg19
BOWTIE=$HOME/bowtie2-2.2.5


# First make files containing lists of SNPs. Each file should be
# named chr<#>.snps.txt.gz  (e.g. chr22.snps.txt.gz) and contain 3 columns:
# position, ref allele, alt allele.
#
# Here is an example of creating SNP files from impute2 output files:
# mkdir -p $SNP_DIR
# for FILE in /data/internal/Yoruba/IMPUTE/CEU/hg19/*impute2.gz; do
#     echo $FILE >&2
#     CHR=`echo $FILE | sed -n 's/^.*\(chr[0-9A-Z]*\).*.impute2.gz$/\1/p'`
#     echo $CHR >&2
#     OUTPUT_FILE=$SNP_DIR/$CHR.snps.txt.gz
#     gunzip -c $FILE | awk '{print $3,$4,$5}' | gzip > $OUTPUT_FILE
# done


# Map reads using bowtie2 (or another mapping tool of your choice)
$BOWTIE/bowtie2 -x $INDEX -1 $DATA_DIR/sim_pe_reads1.fastq.gz \
		-2 $DATA_DIR/sim_pe_reads2.fastq.gz \
    | samtools view -S -b -q 10 - > $DATA_DIR/sim_pe_reads.bam

# Pull out reads that need to be remapped to check for bias
# Use the -p option for paired-end reads.
python $WASP/mapping/find_intersecting_snps.py -p $DATA_DIR/sim_pe_reads.bam $SNP_DIR

# Remap the reads, using same the program and options as before.
# NOTE: If you use an option in the first mapping step that modifies the
# reads (e.g. the -5 read trimming option to bowtie2) you should omit this
# option during the second mapping step here (otherwise the reads will be modified
# twice)!
$BOWTIE/bowtie2 -x $INDEX -1 $DATA_DIR/sim_pe_reads.remap.fq1.gz \
		-2 $DATA_DIR/sim_pe_reads.remap.fq2.gz \
    | samtools view -S -b -q 10 - >  $DATA_DIR/sim_pe_reads.remap.bam

# Use filter_remapped_reads.py to create filtered list of reads that correctly
# remap to same position
python $WASP/mapping/filter_remapped_reads.py -p $DATA_DIR/sim_pe_reads.to.remap.bam \
       $DATA_DIR/sim_pe_reads.remap.bam $DATA_DIR/sim_pe_reads.remap.keep.bam \
       $DATA_DIR/sim_pe_reads.to.remap.num.gz

# Create a merged BAM containing [1] reads that did
# not need remapping [2] filtered remapped reads
samtools merge $DATA_DIR/sim_pe_reads.keep.merged.bam \
	 $DATA_DIR/sim_pe_reads.keep.bam $DATA_DIR/sim_pe_reads.remap.keep.bam

# Sort and index the bam file
samtools sort $DATA_DIR/sim_pe_reads.keep.merged.bam \
	 $DATA_DIR/sim_pe_reads.keep.merged.sorted
	 
samtools index $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam

# Filter out duplicate reads. Use rmdup_pe.py for paired-end reads,
# rmdup.py for single-end reads.
python $WASP/mapping/rmdup_pe.py $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
       $DATA_DIR/sim_pe_reads.keep.rmdup.merged.sorted.bam


