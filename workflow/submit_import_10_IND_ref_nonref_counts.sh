#!/bin/sh
#
#
# This is an example of an SGE cluster submission script
# for creating HDF5 tracks containing HDF5 allele-specific
# read count files. It will need to be updated for your specific purpose.
# It uses the script import_bam_ref_nonref_counts.py
# and a 'taskfile' containing a list of data types and individuals.  
#
#
# preserve shell environment vars:
#$ -V
#
# submit a job-array consisting of following tasks
#$ -t 1-50
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 8GB of mem, use some IO
#$ -l h_vmem=8g,bigio=1

#
# Version 1 - bwa mapped files
#
SCRIPT=$HOME/proj/script/python/10_IND/import_bam_ref_nonref_counts.py
TASKFILE=$HOME/data/10_IND_nb/10_IND_TASKLIST.txt
DATA_TYPE=`head -$SGE_TASK_ID $TASKFILE | tail -1 | awk '{print $2}'`
INDIVIDUAL=`head -$SGE_TASK_ID $TASKFILE | tail -1 | awk '{print $3}'`
INPUT_FILE=/data/share/10_IND/Final_BAMs/${DATA_TYPE}/$INDIVIDUAL.quality.sort.rmdup.bam
IND_STR=

#
# Version 2 - mcmapped files
#
#TASKFILE=/data/share/10_IND_v2/10_IND_tasklist.txt
#SCRIPT=$HOME/proj/script/python/10_IND/import_ref_nonref_counts.py
# DATA_TYPE=`head -$SGE_TASK_ID $TASKFILE | tail -1 | awk '{print $2}'`
# INDIVIDUAL=`head -$SGE_TASK_ID $TASKFILE | tail -1 | awk '{print $3}'`
#INPUT_DIR=/data/share/10_IND_v2/filtered/mappability/$DATA_TYPE/$INDIVIDUAL
#INPUT_FILE=$INPUT_DIR/$INDIVIDUAL.chr*.txt.gz
#IND_STR=$INDIVIDUAL

STAT_SCRIPT=$HOME/proj/genome/python/script/db/set_track_stats.py

echo "$HOSTNAME $DATA_TYPE $INDIVIDUAL" 1>&2

REF_COUNT_TRACK=10_IND/$DATA_TYPE/read_counts/${DATA_TYPE}_${INDIVIDUAL}_AS_ref_count
ALT_COUNT_TRACK=10_IND/$DATA_TYPE/read_counts/${DATA_TYPE}_${INDIVIDUAL}_AS_alt_count
OTHER_COUNT_TRACK=10_IND/$DATA_TYPE/read_counts/${DATA_TYPE}_${INDIVIDUAL}_AS_other_count

READ_COUNT_TRACK=10_IND/$DATA_TYPE/read_counts/${DATA_TYPE}_${INDIVIDUAL}_read_start_count


echo "CMD: python $SCRIPT $IND_STR $REF_COUNT_TRACK $ALT_COUNT_TRACK $OTHER_COUNT_TRACK $READ_COUNT_TRACK $INPUT_FILE"

python $SCRIPT $IND_STR $REF_COUNT_TRACK $ALT_COUNT_TRACK $OTHER_COUNT_TRACK $READ_COUNT_TRACK $INPUT_FILE


if [ "$?" -ne "0" ]; then
    echo "import of read counts failed" 1>&2
    exit 1
fi

python $STAT_SCRIPT $REF_COUNT_TRACK

if [ "$?" -ne "0" ]; then
    echo "set stats of AS ref count track failed" 1>&2
    exit 1
fi


python $STAT_SCRIPT $ALT_COUNT_TRACK

if [ "$?" -ne "0" ]; then
    echo "set stats of AS alt count track failed" 1>&2
    exit 1
fi


python $STAT_SCRIPT $OTHER_COUNT_TRACK

if [ "$?" -ne "0" ]; then
    echo "set stats of AS other count track failed" 1>&2
    exit 1
fi


python $STAT_SCRIPT $READ_COUNT_TRACK

if [ "$?" -ne "0" ]; then
    echo "set stats of read count track failed" 1>&2
    exit 1
fi

echo "done" 1>&2
