#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

echo $'\n'"["`date`"]: Job started."

## Add Modules ########################
module load apps/samtools/1.3
module load apps/perl/5.18.2
#######################################

SAMPLE_ID=$1
SAMPLE_PATH=$2
SCRIPTPATH=$3
REFDIR=$4
TARGETS=$5
COV_OUT_DIR_NAME=$6
DUP_FREE_BAM_DIR_NAME=$7
OUT_PATH=$8
REF_FILE=$9
Library_ID=${10}

REF2="${REFDIR}/${REF_FILE}"

COV_DIR="${OUT_PATH}"
DUP_FREE_BAM_DIR="${OUT_PATH}"

DUP_FREE_BAMhg38="${DUP_FREE_BAM_DIR}/${SAMPLE_ID}"
PILEUP_NODUPS_FILEhg38="${COV_DIR}/${SAMPLE_ID}_nodups_Q25.pileup"

samtools mpileup -q 25 -Q 25 -f $REF2 $DUP_FREE_BAMhg38 > $PILEUP_NODUPS_FILEhg38
perl "${SCRIPTPATH}/Coverage_from_pileup_combinedscript.pl" --inPath1 $COV_DIR --inPath2 $COV_DIR --inFile $PILEUP_NODUPS_FILEhg38 --batchID ${Library_ID} --targets $TARGETS

#clean pileup
rm $PILEUP_NODUPS_FILEhg38

echo $'\n'"["`date`"]: Coverage from pileup complete!!"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
