#! /bin/bash

#/************************* Paras need to be adjusted for diff samples
SAMPLE_PATH="/sharedlustre/IGM/Rach_SC"
OUT_PATH="/sharedlustre/IGM/Rach_SC"

SCRIPT_PATH="/sharedlustre/IGM/Rach_SC"
TARGETS="/sharedlustre/IGM/Rach_SC/hg38/CCDS_20160421.txt"
Library_ID="test.24Aug2016" #e.g. {PI name}.{Batch Number} to identify the library
REFDIR="/sharedlustre/IGM/Rach_SC/hg38"

REF_FILE="hg38.fa"
COV_OUT_DIR_NAME="coverage" #coverage file output directory name. Just a name, it will be placed under the sample directory
DUP_FREE_BAM_DIR_NAME="dup_free_bam" # folder to ouput bam file without any GATK process, this can be input of freebayes
#**************************/

## Submitting jobs ##
INALL=( "${SAMPLE_PATH}"/*bam )

for SampleFile in "${INALL[@]}"; do
	SAMPLE_ID="${SampleFile##*/}"
	 	JOB_IDc="cov_${SAMPLE_ID}"
	    arr2=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPT_PATH}" "$REFDIR" "$TARGETS" "${COV_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${OUT_PATH}" "${REF_FILE}" "${Library_ID}")
        echo "${SAMPLE_ID}"
		qsub -l "h_vmem=50G" -M "helen.griffin@ncl.ac.uk" -N "${JOB_IDc}" ${SCRIPT_PATH}/Coverage_from_pileup_submit_perl3.sh "${arr2[@]}"
done
