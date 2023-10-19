#!/bin/bash
#PBS -N ALIGN_MINIMAP2_CHM13
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l jobfs=400GB
#PBS -l walltime=48:00:00
#PBS -l wd

# Load required modules
module load minimap2/2.22
module load samtools/1.12

# Define reference paths and directories
FEMALE_REF=/path/to/chm13_reference_dir/chm13.draft_v1.1.fasta
MALE_REF=/path/to/chm13_reference_dir/chm13_v1.1_with_HG002_chrY.masked_PARs.fasta
SAMPLE_GENDER=/path/to/sex_info.tab
SAMPLE_DATA_DIR=`ls -tr -d /path/to/basecalls/*/*/${SAMPLE_ID}/* | tail -n1`
TMP_DIR=/path/to/${SAMPLE_ID}_minimap2
MINIMAP2_VERSION=`minimap2 --version | awk '{print "v"$1}'`
DATE=`date "+%Y%m%d"`
OUT_DIR=/path/to/minimap2/${MINIMAP2_VERSION}/${SAMPLE_ID}/${DATE}
FASTQ_PASS=`ls -tr /path/to/basecalls/*/*/${SAMPLE_ID}/*/*pass.fastq.gz | tail -n1`
TMP_LOG=align_${SAMPLE_ID}.minimap2.log
LOG=${OUT_DIR}/align_${SAMPLE_ID}.minimap2.log
num_threads=${PBS_NCPUS}

# Terminate script
die() {
  echo "$1" >&2
  echo
  exit 1
}

# Check for required variable
[ -z "${SAMPLE_ID}" ] && die "Usage: qsub -v SAMPLE_ID=PBXXXXX ./align.minimap2.sh"

# Check for log file and temporary directory existence
test -e ${LOG} && die "Log file ${LOG} already exists, indicating that aligning was already run. Please remove that first, if you want to re-align."
test -d ${TMP_DIR} && die "Temporary directory ${TMP_DIR} already exists, indicating an unfinished business. Please delete it first. Exiting."

# Check for the existence of sample gender file
test -e ${SAMPLE_GENDER} || die "Sample gender file is not present at ${SAMPLE_GENDER}. Exiting."

# Create temporary directories and local scratch storage
mkdir ${TMP_DIR} || die "Creating ${TMP_DIR} failed. Exiting."
cd ${TMP_DIR} || die "cd to ${TMP_DIR} failed. Exiting."
mkdir fastq || die "Creating tmp fastq dir failed. Exiting."

# Deduce sample gender
GENDER=$(grep ${SAMPLE_ID} ${SAMPLE_GENDER} | head -n1 | awk '{print $3}')
if [ "${GENDER}" = "XX" ]; then
  REFERENCE=${FEMALE_REF}
  echo "Sample ${SAMPLE_ID} deduced to be female" | tee -a ${TMP_LOG}
elif [ "${GENDER}" = "XY" ]; then
  REFERENCE=${MALE_REF}
  echo "Sample ${SAMPLE_ID} deduced to be male" | tee -a ${TMP_LOG}
else
  die "Could not deduce gender. Deduced value: ${GENDER}"
fi

# Perform alignment
echo "Aligning sample ${SAMPLE_ID}" | tee -a ${TMP_LOG}
echo -n "Minimap2 version: " | tee -a ${TMP_LOG}
minimap2 --version 2>&1 | tee -a ${TMP_LOG}
echo -n "Samtools version: " | tee -a ${TMP_LOG}
samtools --version 2>&1 | tee -a ${TMP_LOG}

# Copy the reference and FASTQ file to the local scratch storage
REF_LOCAL=$PBS_JOBFS/ref.fa
FASTQ_LOCAL=$PBS_JOBFS/reads.fastq.gz
echo "Copying reference $REFERENCE to $REF_LOCAL" | tee -a ${TMP_LOG}
cp $REFERENCE $REF_LOCAL || die "Copying $REFERENCE to $REF_LOCAL failed"
echo "Copying ${SAMPLE_ID}_pass.fastq to $FASTQ_LOCAL" | tee -a ${TMP_LOG}
cp $FASTQ_PASS $FASTQ_LOCAL || die "Copying ${FASTQ_PASS} to $FASTQ_LOCAL failed"

# Minimap2 Alignment
echo "Minimap2 Alignment" | tee -a ${TMP_LOG}
( /usr/bin/time -v minimap2 -x map-ont -a -t$num_threads --secondary=no --MD $REF_LOCAL $FASTQ_LOCAL > ${SAMPLE_ID}_pass.sam || die "Minimap2 failed. Exiting." ) 2>&1 | tee -a ${TMP_LOG}

# Samtools sort and index
echo "Samtools sort" | tee -a ${TMP_LOG}
( /usr/bin/time -v samtools sort -@$num_threads -m 3G ${SAMPLE_ID}_pass.sam > ${SAMPLE_ID}_pass.bam || die "Samtools sort failed. Exiting." ) 2>&1 | tee -a ${TMP_LOG}
echo "Samtools index" | tee -a ${TMP_LOG}
( /usr/bin/time -v samtools index ${SAMPLE_ID}_pass.bam || die "Samtools index failed. Exiting." ) 2>&1 | tee -a ${TMP_LOG}

# Create output directories and move data
mkdir -p ${OUT_DIR}
echo "Moving data" | tee -a ${TMP_LOG}
mv ${SAMPLE_ID}_pass.bam ${SAMPLE_ID}_pass.bam.bai ${OUT_DIR}/ || die "Moving ${SAMPLE_ID}_pass.bam ${SAMPLE_ID}_pass.bam.bai to ${OUT_DIR} failed. Exiting."

echo "Completed successfully" | tee -a ${TMP_LOG}
mv ${TMP_LOG} ${LOG} || die "Moving ${TMP_LOG} to ${LOG} failed. Exiting."
rm -r ${TMP_DIR}

