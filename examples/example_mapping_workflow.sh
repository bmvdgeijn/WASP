
# Set these environment vars to point to
# your local installation of WASP
WASP=$HOME/proj/WASP
DATA_DIR=$WASP/examples/example_data

# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2
INDEX=$HOME/data1/external/GRC37/combined/bowtie2_index/hg37
BOWTIE=$HOME/anaconda2/bin/bowtie2


$WASP/snp2h5/snp2h5 --chrom $DATA_DIR/chromInfo.hg19.txt \
		    --format impute \
		    --snp_index $DATA_DIR/genotypes/snp_index.h5 \
		    --geno_prob $DATA_DIR/genotypes/geno_probs.h5 \
		    --snp_tab $DATA_DIR/genotypes/snp_tab.h5 \
		    --haplotype $DATA_DIR/genotypes/haps.h5 \
		    --samples $DATA_DIR/genotypes/YRI_samples.txt \
		    $DATA_DIR/genotypes/chr*.hg19.impute2.gz \
		    $DATA_DIR/genotypes/chr*.hg19.impute2_haps.gz


# Map reads using bowtie2 (or another mapping tool of your choice)
$BOWTIE -x $INDEX -1 $DATA_DIR/sim_pe_reads1.fastq.gz \
	-2 $DATA_DIR/sim_pe_reads2.fastq.gz \
    | samtools view -S -b -q 10 - > $DATA_DIR/sim_pe_reads.bam

# Pull out reads that need to be remapped to check for bias
# Use the -p option for paired-end reads.
python $WASP/mapping/find_intersecting_snps.py \
       --is_paired_end \
       --output_dir $DATA_DIR  \
       --snp_index $DATA_DIR/genotypes/snp_index.h5 \
       --snp_tab $DATA_DIR/genotypes/snp_tab.h5 \
       --haplotype $DATA_DIR/genotypes/haps.h5 \
       --samples $DATA_DIR/H3K27ac/samples.txt \
       $DATA_DIR/sim_pe_reads.bam

# Remap the reads, using same the program and options as before.
# NOTE: If you use an option in the first mapping step that modifies the
# reads (e.g. the -5 read trimming option to bowtie2) you should omit this
# option during the second mapping step here (otherwise the reads will be modified
# twice)!
$BOWTIE -x $INDEX -1 $DATA_DIR/sim_pe_reads.remap.fq1.gz \
	          -2 $DATA_DIR/sim_pe_reads.remap.fq2.gz \
    | samtools view -S -b -q 10 - >  $DATA_DIR/sim_pe_reads.remap.bam

# Use filter_remapped_reads.py to create filtered list of reads that correctly
# remap to same position
python $WASP/mapping/filter_remapped_reads.py \
       $DATA_DIR/sim_pe_reads.to.remap.bam \
       $DATA_DIR/sim_pe_reads.remap.bam \
       $DATA_DIR/sim_pe_reads.remap.keep.bam

# Create a merged BAM containing [1] reads that did
# not need remapping [2] filtered remapped reads
samtools merge $DATA_DIR/sim_pe_reads.keep.merged.bam \
	 $DATA_DIR/sim_pe_reads.keep.bam $DATA_DIR/sim_pe_reads.remap.keep.bam

# Sort and index the bam file
samtools sort $DATA_DIR/sim_pe_reads.keep.merged.bam \
	 -o $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam
	 
samtools index $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam

# Filter out duplicate reads. Use rmdup_pe.py for paired-end reads,
# rmdup.py for single-end reads.
python $WASP/mapping/rmdup_pe.py $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
       $DATA_DIR/sim_pe_reads.keep.rmdup.merged.sorted.bam


