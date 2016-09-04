Pipeline for mappability filtering
==================================

This directory contains scripts that can be used to eliminate mapping
bias from mapped allele-specific reads.  First, reads are mapped
normally using a mapper chosen by the user (must output BAM or SAM
format).  Then mapped reads that overlap single nucleotide
polymorphisms (SNPs) are identified. For each read that overlaps a
SNP, its genotype is swapped with that of the other allele and the
read is re-mapped. Re-mapped reads that fail to map to exactly the
same location in the genome are discarded.


SnakeMake
---------

We now provide a Snakemake workflow that can be used to run the entire
mappability filtering pipeline. For more information see the
[Snakemake README](README.snakemake.md)

Step 1:
-------

Create input files containing genetic variants (SNPs). The WASP
pipeline can take SNP input files in HDF5 format or in a text-based format.
The HDF5 format is preferred as it contains phase information that improves
the performance of the mappability filtering pipeline. (If you do not
have phasing information, then the text-based format should be used.)

### Creating HDF5 SNP files

Convert SNP data to HDF5 format using the program snp2h5.  HDF5 files
are an efficient binary data format used by WASP scripts. The snp2h5
program can take VCF or IMPUTE2 files as input.

If VCF files are used, the sample names are read directly from the
header lines of the VCF file. If IMPUTE2 files are used, then
the names of the samples should be provided in a separate text file.

    	# using an IMPUTE input file:
       ./snp2h5/snp2h5 --chrom example_data/chromInfo.hg19.txt \
              --format impute \
              --geno_prob example_data/geno_probs.h5 \
              --snp_index example_data/snp_index.h5 \
              --snp_tab example_data/snp_tab.h5 \
              --haplotype example_data/haps.h5 \
              --samples samples_names.txt
              example_data/genotypes/chr*.hg19.impute2.gz \
              example_data/genotypes/chr*.hg19.impute2_haps.gz

       # using VCF files:
       ./snp2h5/snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \
             --format vcf \
             --haplotype haplotypes.h5 \
             --snp_index snp_index.h5 \
             --snp_tab   snp_tab.h5 \
             data/1000G/ALL.chr*.vcf.gz

       # read separate VCF files that contain genotype likelihoods (GL)
       # or genotype probabilities (GP) if these were not in the main
       # VCF file
       ./snp2h5/snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz \
            --format vcf \
            --geno_prob geno_probs.h5 \
            1000G/supporting/genotype_likelihoods/shapeit2/ALL.chr*.gl.vcf.gz

### Creating text-based SNP files

The text-based input files have three space-delimited columns
(position, ref_allele, alt_allele), and one input file per chromosome.
The filenames must contain the name of the chromosome (e.g. chr2).

We provide example scripts that can be used to create these files
from IMPUTE or VCF files:

       # get SNPs from IMPUTE files:
       ./mapping/extract_impute_snps.sh example_data/genotypes/ output_snp_dir/

       # get SNPs from VCF files:
       ./mapping/extract_vcf_snps.sh example_data/genotypes/ output_snp_dir/



Step 2: 
-------

Map the fastq files using your favorite mapper/options and filter for
quality using a cutoff of your choice

### Example:
       bowtie2 -x bowtie2_index/hg37 -1 ${SAMPLE_NAME}_1.fq.gz \
               -2 ${SAMPLE_NAME}_2.fq.gz \
               | samtools view -b -q 10 - > map1/${SAMPLE_NAME}.bam
       samtools sort -o map1/${SAMPLE_NAME}.sort.bam map1/${SAMPLE_NAME}.bam
       samtools index map1/${SAMPLE_NAME}.sort.bam


Step 3:
-------

Use find_intersecting_snps.py to identify reads that may have mapping biases

### Example:
       python mapping/find_intersecting_snps.py \
              --is_paired_end \
              --is_sorted \
              --output_dir find_intersecting_snps \
              --snp_tab snp_tab.h5 \
              --snp_index snp_index.h5 \
              --haplotype haplotype.h5 \
              --samples my_samples.txt \
              map1/${SAMPLE_NAME}.sort.bam

#### Usage:
       positional arguments:
            bam_filename          Coordinate-sorted input BAM file containing
                                  mapped reads.

       optional arguments:
             -h, --help            show this help message and exit
             --is_paired_end, -p   Indicates that reads are paired-end (default
                                   is single).
             --is_sorted, -s       Indicates that the input BAM file is
	                           coordinate-sorted (default is False).
             --max_seqs MAX_SEQS   The maximum number of sequences with 
                                   different allelic combinations to consider
                                   remapping (default=64). Read pairs wi
                                   more allelic combinations than MAX_SEQs are
                                   discarded
             --max_snps MAX_SNPS   The maximum number of SNPs allowed to
                                   overlap a read before discarding the read.
                                   Allowing higher numbers will decrease speed
                                   and increase memory usage (default=6).
             --output_dir OUT_DIR  Directory to write output files to. If not
                                   specified, output files are written to the
                                   same directory as the input BAM file.
             --snp_dir SNP_DIR     Directory containing SNP text files This
                                   directory should contain one file per
                                   chromosome named like chr<#>.snps.txt
                                   Each file should contain 3 columns: position
                                   RefAllele AltAllele. This option should
                                   only be used if --snp_tab, --snp_index,
                                   and --haplotype arguments are not used.
                                   If this argument is provided, all possible
                                   allelic combinations are used (rather
                                   than set of observed haplotypes)
             --snp_tab SNP_TABLE_H5_FILE
                                   Path to HDF5 file to read SNP info from. Each
                                   row of SNP table contains SNP name, position,
                                   allele1, allele2.
             --snp_index SNP_INDEX_H5_FILE
                                   Path to HDF5 file containing SNP index.
                                   The SNP index is used to convert the
                                   genomic position of a SNP to
                                   its corresponding row in the haplotype and
                                   snp_tab HDF5 files.
             --haplotype HAPLOTYPE_H5_FILE
                                   Path to HDF5 file to read phased haplotypes
                                   from. When generating alternative reads use
                                   known haplotypes from this file rather than
                                   all possible allelic combinations.
             --samples SAMPLES     Use only haplotypes and SNPs that are
                                   polymorphic in these samples. SAMPLES can
                                   either be a comma-delimited string of sample
                                   names or a path to a file with one
                                   sample name per line (file is assumed to be
                                   whitespace-delimited and first column is
                                   assumed to be sample name). Sample names
                                   should match those present in the
                                   --haplotype file. Samples are ignored if no
                                   haplotype file is provided.


#### Output:
         PREFIX.keep.bam - bamfile with reads that did not intersect SNPs
                          or indels that can be kept without remapping
         PREFIX.to.remap.bam - bamfile with original reads that overlapped SNPs
                          that need to be remapped
         PREFIX.remap.fq.gz - fastq file containing the reads with flipped
                          alleles to remap. If paired-end option is used
                          two files ending with .fq1.gz and .fq2.gz are output.
         (PREFIX is the name of the input file, excluding the trailing .bam)
	
         Note: Reads that overlap indels are currently excluded and
         will not be present in any of the 'remap' files or the
         input.keep.bam file. For this reason the total number of reads
         will not add up to the number of reads provided in the
         input.sort.bam file.


Step 4
-----
Map the PREFIX.remap.fq.gz using the same mapping arguments used in
Step 1. Note that the arguments should be exactly the same as those in
Step 1 EXCEPT for arguments that directly modify the reads that are
used by the aligner. For example the read trimming arguments to bowtie
(-3 and -5 arguments) should be used in Step 1 ONLY, because they
modify the reads that are output by bowtie.

### Example:
         bowtie2 -x bowtie2_index/hg37 \
                   -1 find_intersecting_snps/${SAMPLE_NAME}_1.remap.fq.gz \
                   -2 find_intersecting_snps/${SAMPLE_NAME}_2.remap.fq.gz \
               | samtools view -b -q 10 - > map2/${SAMPLE_NAME}.bam
         samtools sort -o map2/${SAMPLE_NAME}.sort.bam map2/${SAMPLE_NAME}.bam
         samtools index map2/${SAMPLE_NAME}.sort.bam


Step 5
------
Use filter_remapped_reads.py to filter out reads where one or more
of the allelic versions of the reads fail to map back to the same
location as the original read.

#### Usage:
         filter_remapped_reads.py [-h] to_remap_bam remap_bam keep_bam
       
         positional arguments:
           to_remap_bam  input BAM file containing original set of reads that
           needed to be remapped after having their alleles
			 flipped. This file is output by the
			 find_intersecting_snps.py script.
           remap_bam     input BAM file containing remapped reads (with flipped
                         alleles)
           keep_bam      output BAM file to write filtered set of reads to

#### Example:
         python mapping/filter_remapped_reads.py \
           find_intersection_snps/${SAMPLE_NAME}.to.remap.bam \
           map2/${SAMPLE_NAME}.sort.bam \
           filter_remapped_reads/${SAMPLE_NAME}.keep.bam


Step 6
------

Merge  ${SAMPLE_NAME}.keep.bam and ${SAMPLE_NAME}.remap.keep.bam 
can be merged for a complete set of mappability-filtered aligned reads.
The merged file should then be sorted and indexed:

#### Example:
         samtools merge merge/${SAMPLE_NAME}.keep.merge.bam \
                  filter_remapped_reads/${SAMPLE_NAME}.keep.bam  \
                  find_intersecting_snps/${SAMPLE_NAME}.keep.bam
         samtools sort -o  merge/${SAMPLE_NAME}.keep.merge.sort.bam \
                  merge/${SAMPLE_NAME}.keep.merge.bam 
         samtools index ${SAMPLE_NAME}.keep.merged.sort.bam


Step 7
------

Filter duplicate reads. Programs such as samtools rmdup introduce bias
when they filter duplicate reads because they retain the read with the
highest score (which usually matches the reference). We provide a
script rmdup.py which performs unbiased removal of duplicate
reads. The script discards duplicate reads at random (independent of
their score). The input BAM or SAM file must be sorted.

#### Usage:
         # for single end reads:
         python rmdup.py <sorted.input.bam> <output.bam>
         # for paired-end reads:
         python rmdup_pe.py <sorted.input.bam> <output.bam>
	
## Testing

To run the tests, execute `py.test` from within the mapping directory.
The tests currently require bowtie2 and samtools to be in the PATH.
