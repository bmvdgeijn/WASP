CHT
===

Combined Haplotype Test

===

Input Format



===

Using the test script

python combined_test.py [options] input_files.txt output_file.txt
       input_file.txt-
		A file containing a list of the input files to be tested
       output_file.txt-
		Each line is the output from a tested SNP with the format:

chromosome | position | chi-square value | best alpha parameter | best beta parameter | overdispersion parameter | # of allele specific reads

Options:
	-a	perform the allele specific portion of the test only
	-d	perform the total read depth portion of the test only
	-s	randomly permute the genotypes and allele specific counts
	-d <float> choose a total read depth overdispersion parameter (default=0.00001)
	-o <float> choose an allele specific overdispersion parameter (default=0.00001)
	-m <int>   choose a threshold for minimum number of allele specific reads (default 0)
	-e <float> choose a sequencing error rate in the reads (default 0.005)