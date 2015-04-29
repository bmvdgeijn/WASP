#
# convert snp files to HDF5 format
#
./snp2h5/snp2h5 --chrom test_data/chromInfo.hg19.txt \
	      --format impute \
	      --geno_prob test_data/geno_probs.h5 \
	      --snp_index test_data/snp_index.h5 \
	      --snp_tab test_data/snp_tab.h5 \
	      --haplotype test_data/haps.h5 \
	      test_data/genotypes/chr*.hg19.impute2.gz \
	      test_data/genotypes/chr*.hg19.impute2_haps.gz

#
# convert FASTA files to HDF5 format
# Note the HDF5 sequence files are only used for GC content
# correction part of CHT. This step can be ommitted if
# GC-content correction is not used.
#
./snp2h5/fasta2h5 --chrom test_data/chromInfo.hg19.txt \
	--seq test_data/seq.h5 \
	/data/external_public/reference_genomes/hg19/chr*.fa.gz


INDIVIDUAL=18505
#
# read BAM files for individual 18505 and write read counts to HDF5 files
#
./bam2h5/bam2h5.py --chrom test_data/chromInfo.hg19.txt \
	      --snp_index test_data/snp_index.h5 \
	      --snp_tab test_data/snp_tab.h5 \
	      --haplotype test_data/haps.h5 \
	      --samples test_data/H3K27ac/samples.txt \
	      --individual $INDIVIDUAL \
	      --ref_as_counts test_data/ref_as_counts.$INDIVIDUAL.h5 \
	      --alt_as_counts test_data/alt_as_counts.$INDIVIDUAL.h5 \
	      --other_as_counts test_data/other_as_counts.$INDIVIDUAL.h5 \
	      --read_counts test_data/read_counts.$INDIVIDUAL.h5 \
	      test_data/H3K27ac/$INDIVIDUAL.chr*.keep.rmdup.bam


#
# create combined haplotype test input file for individual 18505
#
python CHT/extract_haplotype_read_counts.py \
       --chrom test_data/chromInfo.hg19.txt \
       --snp_index test_data/snp_index.h5 \
       --snp_tab test_data/snp_tab.h5 \
       --geno_prob test_data/geno_probs.h5 \
       --haplotype test_data/haps.h5 \
       --samples test_data/H3K27ac/samples.txt \
       --individual $INDIVIDUAL \
       --ref_as_counts test_data/ref_as_counts.$INDIVIDUAL.h5 \
       --alt_as_counts test_data/alt_as_counts.$INDIVIDUAL.h5 \
       --other_as_counts test_data/other_as_counts.$INDIVIDUAL.h5 \
       --read_counts test_data/read_counts.$INDIVIDUAL.h5 \
       test_data/H3K27ac/chr22.peaks.txt \
       > test_data/H3K27ac/haplotype_read_counts.$INDIVIDUAL.h5




