#!/bin/bash
#
#$ -q iblm.q
#
#$ -o $HOME/sge/out
#$ -e $HOME/sge/err
#
#$ -t 1-126
#$ -V
#

WASP_DIR=$HOME/proj/WASP/

DATA_DIR=/iblm/netapp/home/gmcvicker/data1/gmcvicker/data/WASP_power_calcs

TASK_LIST=$DATA_DIR/power_sim_params.txt

# add one to task id to skip header line
LINE=$((SGE_TASK_ID + 1))

# read simulation parameters from task file
REF_ALLELE_FREQ=`sed -n ${LINE}p $TASK_LIST | awk '{print $2}'`
EFFECT_SIZE=`sed -n ${LINE}p $TASK_LIST | awk '{print $3}'`
SAMPLE_SIZE=`sed -n ${LINE}p $TASK_LIST | awk '{print $4}'`

# create output directory for simulation results
OUT_DIR=$DATA_DIR/sim_af${REF_ALLELE_FREQ}_es${EFFECT_SIZE}_ss${SAMPLE_SIZE}
mkdir -p $OUT_DIR

# run simulation
Rscript --vanilla $WASP_DIR/CHT/power_analysis/sim_haplotype_read_counts.R \
  --n.sample $SAMPLE_SIZE \
  --allele.freq $REF_ALLELE_FREQ \
  --effect.size $EFFECT_SIZE \
  --output.dir $OUT_DIR


# make list of input files for combined haplotype test
ls $OUT_DIR/sim_hap_read_counts*.txt > $OUT_DIR/input_files.txt

# estimate dispersion parameters from simulation data
python $WASP_DIR/CHT/fit_as_coefficients.py $OUT_DIR/input_files.txt $OUT_DIR/as_coef.txt

python $WASP_DIR/CHT/fit_bnb_coefficients.py --sample 1000 $OUT_DIR/input_files.txt $OUT_DIR/bnb_coef.txt

# run combined haplotype test
# python $WASP_DIR/CHT/combined_test.py \
#         --bnb_disp $OUT_DIR/bnb_coef.txt \
#         --as_disp $OUT_DIR/as_coef.txt \
# 	--as_only \
#         $OUT_DIR/input_files.txt \
#         $OUT_DIR/cht_results.as_only.txt

# python $WASP_DIR/CHT/combined_test.py \
#         --bnb_disp $OUT_DIR/bnb_coef.txt \
#         --as_disp $OUT_DIR/as_coef.txt \
# 	--bnb_only \
#         $OUT_DIR/input_files.txt \
#         $OUT_DIR/cht_results.bnb_only.txt

python $WASP_DIR/CHT/combined_test.py \
        --bnb_disp $OUT_DIR/bnb_coef.txt \
        --as_disp $OUT_DIR/as_coef.txt \
        $OUT_DIR/input_files.txt \
        $OUT_DIR/cht_results.txt
