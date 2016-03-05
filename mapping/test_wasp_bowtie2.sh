
# INDEX=$HOME/data1/external/GRC38/combined/bowtie2_index/hg38

INDEX=../example_data/test_genome
bowtie2-build ../example_data/test_genome.fa $INDEX

MAP_CMD1="bowtie2 -x $INDEX -1 FASTQ1 -2 FASTQ2 | samtools view -S -b -q 10 - > OUTPUT_BAM"

./wasp --map_command "$MAP_CMD1" --fastq1 ../example_data/test_reads1.fq --fastq2 ../example_data/test_reads2.fq --output_dir ../test ../example_data/test_*snps.txt




