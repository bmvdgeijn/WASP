## Snakemake Mappability Filtering Pipeline

[Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) is a
workflow management system, designed to streamline the execution of
software pipelines. We now provide a Snakemake rule file that can be
used to run the entire Mappability Filtering Pipeline.

For a more complete description of Snakemake see the
[Snakemake tutorial](http://snakemake.bitbucket.org/snakemake-tutorial.html).

## Installing Snakemake

Snakemake requires python3, however the CHT pipeline requires
python2. For this reason, if you are using
[Anaconda](https://www.continuum.io/downloads), it is recommended that
you create a [python3
environment](http://conda.pydata.org/docs/py2or3.html#create-a-python-3-5-environment). For example you can create a python3.5 Anaconda environment with the following shell command (this only needs to be done once):

        conda create -n py35 python=3.5 anaconda

You can then activate the py35 environment, and install the latest version of
Snakemake with the following commands:

        source activate py35
        conda install snakemake

Then when you want to switch back to your default (e.g. python2) environment
do the following:

        source deactivate



## Configuring the Mappability Filtering Pipeline

The rules for the Snakemake tasks are defined in the [Snakefile](Snakefile).

Configuration parameters for this Snakefile are read from the YAML file
[snake_conf.yaml](snake_conf.yaml).

Before running Snakemake, edit this file to specify the location
of all of the input directories and files that will be used by the pipeline.
This includes locations of the impute2 or VCF SNP files, input BAM files etc.

Importantly you must set `wasp_dir` to point to the location of WASP
on your system, and set `py2` to setup the environment
for python (e.g. by modifying your PATH) and call the
appropriate interpreter.  This is necessary because Snakemake is run
using python3, but most of the scripts require python2.


## Running the Mappability Filtering pipeline

Snakemake can be run as a single process or on a compute cluster with
multiple jobs running simultaneuously. To run Snakemake on a single node
you could do something like the following:

        source activate py35
        cd $WASP_DIR/CHT
        snakemake

We provide a script [run_snakemake.sh](run_snakemake.sh) to run Snakemake
on a SGE compute cluster. You must be in a python3 environment to run this
script, and the script must be run from a job submission host.

        source activate py35
        cd $WASP_DIR/CHT
        ./run_snakemake.sh

It should be possible to make simple modifications to this script to
run on queue management systems other than SGE (e.g. LSF or Slurm).


You should Snakemake from within a [Screen](https://www.gnu.org/software/screen/) virtual terminal or using [nohup](https://en.wikipedia.org/wiki/Nohup) so
that if you are disconnected from the cluster, Snakemake will continue to run.

## Output from the Mappability Filtering Pipeline

The Snakemake pipeline creates the following directories under the
output_dir specified in snake_conf.yaml. The rmdup and as_counts directories
contain the final output.

* map1 - results from first mapping of reads to genome
* map1_sort - sorted BAM files from first mapping
* find_intersecting_snps - results from find_intersecting_snps.py
* map2 - results from second mapping of reads to genome
* map2_sort - sorted BAM files from second mapping
* filter_remapped_reads - results from filter_remapped_reads.py
* merge - merged reads that should be kept because they did not
        overlap a SNP or because they overlapped a SNP and all
	alleles remapped to the same location.
* rmdup - Final BAM files with duplicate reads removed. Only *sort* files
        need to be kept.
* as_counts - allele specific read counts at polymorphic SNPs

Once the pipeline is complete, most of the files can be removed to save space:

       # remove unsorted files from rmdup dir:
       ls rmdup/ | grep -v sort | xargs rm
       # remove intermediate files and directories:
       rm -rf map1 map1_sort find_intersecting_snps map2 map2_sort filter_remapped_reads merge


## Debugging the Snakemake jobs

By default Snakemake will write an output and error file for each job
to your home directory. These files will be named like `snakejob.<rulename>.<job_num>.sh.{e|o}<sge_jobid>`. For example:

   	# contains error output for extract_haplotype_read_counts rule:
   	snakejob.find_intersecting_snps_paired_end.13.sh.e4507125

If a rule fails, you should check the appropriate output file to see what
error occurred. A major benefit of Snakemake is that if you re-run snakemake
after a job fails it will pickup where it left off.

You should make sure that all of the existing jobs have been killed before
re-starting the pipeline.

