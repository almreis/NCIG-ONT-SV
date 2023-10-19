#!/bin/bash
#PBS -N SLOW5_GUPPY
#PBS -q gpuvolta
#PBS -l ncpus=48
#PBS -l ngpus=4
#PBS -l mem=384GB
#PBS -l walltime=48:00:00
#PBS -l wd

# Set error handling
set -e
set -o pipefail

# Load required modules
module load /path/to/modules/slow5-guppy/6.0.1
module load /path/to/modules/slow5tools/0.6.0

# Define variables
SAMPLE_DATA_DIR=/path/to/${SAMPLE_ID}
TMP_LOG=hac_basecall_${SAMPLE_ID}.log
BLOW5_FILE=${SAMPLE_DATA_DIR}/${SAMPLE_ID}_merged.blow5
BLOW5_FILE_FALLBACK=${SAMPLE_DATA_DIR}/${SAMPLE_ID}_reads.blow5
BLOW5_LOG=${SAMPLE_DATA_DIR}/slow5_${SAMPLE_ID}.log
num_threads=${PBS_NCPUS}
GUPPY_VERSION=`guppy_basecaller --version | grep "Version" | sed 's/.*Version ' | sed 's/+.*' | awk '{print "guppy_"$1}'`
DATE=`date "+%Y%m%d"`
MODEL=hac
OUTDIR=/path/to/${GUPPY_VERSION}/${MODEL}/${SAMPLE_ID}/${DATE}
LOG=${OUTDIR}/hac_basecall_${SAMPLE_ID}.log
TMP_DIR=/path/to/${SAMPLE_ID}_${GUPPY_VERSION}_${MODEL}

# Function for script termination
die() {
  echo "$1" >&2
  echo
  exit 1
}

# Function for sanity check
sanity_check() {
  # Check and compare the number of reads in FASTQ and SLOW5
  if grep "SLOW5_read_count" ${BLOW5_LOG}; then
    NUM_SLOW5_READS=$(grep "SLOW5_read_count" ${BLOW5_LOG} | awk '{print $2}')
  elif grep -q "estimated raw reads in FAST5" ${BLOW5_LOG}; then
    NUM_SLOW5_READS=$(grep "estimated raw reads in FAST5" ${BLOW5_LOG} | awk '{print $1}')
  else
    die "Cannot deduce the number of records using ${BLOW5_LOG}"
  fi

  NUM_READS_PASS=$(echo "$(cat basecalls/pass/*.fastq | wc -l)/4" | bc)
  NUM_READS_FAIL=$(echo "$(cat basecalls/fail/*.fastq | wc -l)/4" | bc)
  NUM_READS=$(echo "$NUM_READS_PASS+$NUM_READS_FAIL" | bc)

  # Compare the number of reads in FASTQ and SLOW5
  if [ ${NUM_READS} -ne ${NUM_SLOW5_READS} ]; then
    echo "WARNING: Sanity check failed. $NUM_READS in FASTQ, but $NUM_SLOW5_READS reads in SLOW5"
    FASTQ_PERCENT=$(echo "(((($NUM_READS))/($NUM_SLOW5_READS))*100)" | bc -l)
    FASTQ_PERCENTINT=$(echo "$FASTQ_PERCENT/1" | bc)
    if [ $FASTQ_PERCENTINT -lt 99 ]; then
      die "ERROR: Sanity check failed. Only ${FASTQ_PERCENTINT}% of reads are in FASTQ."
    else
      echo "${FASTQ_PERCENTINT}% reads are in FASTQ."
    fi
  else
    echo "$NUM_READS in FASTQ, $NUM_SLOW5_READS reads in SLOW5"
  fi

  echo "Reads passed QC: $NUM_READS_PASS"
  echo "Reads failed QC: $NUM_READS_FAIL"

  PASS_PERCENT=$(echo "(((($NUM_READS_PASS))/($NUM_READS))*100)" | bc -l)
  PASS_PERCENTINT=$(echo "$PASS_PERCENT/1" | bc)
  
  if [ $PASS_PERCENTINT -lt 80 ]; then
    echo "WARNING: Only ${PASS_PERCENT}% reads have passed QC."
  else
    echo "${PASS_PERCENT}% reads have passed QC."
  fi
}

# Check for required variable
[ -z "${SAMPLE_ID}" ] && die "Usage: qsub -v SAMPLE_ID=PBXXXXX ./slow5_guppy.pbs.sh"

# Check for log file and temporary directory existence
test -e ${LOG} && die "Log file ${LOG} already exists, indicating that high accuracy basecalling was already run. Please remove that first, if you want to re basecall."
test -d ${TMP_DIR} && die "Temporary directory ${TMP_DIR} already exists, indicating unfinished business. Please delete it first. Exiting."

# Check for the existence of BLOW5 files
if test -e ${BLOW5_FILE}; then
  test -e ${BLOW5_LOG} || die "BLOW5 conversion log ${BLOW5_LOG} does not exist. Please convert FAST5 to BLOW5 first."
else
  echo "${BLOW5_FILE} not found. Checking for fallback ${BLOW5_FILE_FALLBACK}"
  BLOW5_FILE=${BLOW5_FILE_FALLBACK}
  test -e ${BLOW5_FILE} || die "Merged BLOW5 file ${BLOW5_FILE} does not exist. Please convert FAST5 to BLOW5 first."
  NUM_SLOW5_READS=$(slow5tools skim ${BLOW5_FILE} --rid | wc -l)
  echo "SLOW5_read_count $NUM_SLOW5_READS" | tee -a ${BLOW5_LOG}
fi

# Create a temporary directory and move to it
mkdir ${TMP_DIR} || die "Creating ${TMP_DIR} failed. Exiting."
cd ${TMP_DIR} || die "cd to ${TMP_DIR} failed. Exiting."

# Basecalling with Guppy
echo "basecalling sample ${SAMPLE_ID}" | tee -a ${TMP_LOG}
echo -n "Guppy version: " | tee -a ${TMP_LOG}
guppy_basecaller --version 2>&1 | tee -a ${TMP_LOG}
echo "Basecalling" | tee -a ${TMP_LOG}
# R10.4
( /usr/bin/time -v guppy_basecaller -c dna_r10.4_e8.1_hac_prom.cfg -i ${BLOW5_FILE} -s basecalls/ -r --device cuda:all --slow5_threads ${PBS_NCPUS} || die "Basecalling failed. Exiting." ) 2>&1 | tee -a ${TMP_LOG}

# Generate FASTQ files
echo "Generating ${SAMPLE_ID}_pass.fastq" | tee -a ${TMP_LOG}
cat basecalls/pass/*.fastq | pigz -p ${num_threads} > ${SAMPLE_ID}_pass.fastq.gz || die "Generating ${SAMPLE_ID}_pass.fastq.gz failed. Exiting."
echo "Generating ${SAMPLE_ID}_fail.fastq" | tee -a ${TMP_LOG}
cat basecalls/fail/*.fastq | pig

