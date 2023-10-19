#!/bin/bash
#PBS -N CNVPYTOR 
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -l wd

###################################################################

module load pythonlib/3.9.2
module load kentutils/0.0

export PATH=$PATH:/home/913/ar5564/.local/bin

# Define the path to the BAM file
BAM=$(ls -1tr /path/to/alignments/minimap2/v2.22-r1101/${SAMPLE_ID}/*/${SAMPLE_ID}_pass.bam | tail -n1)

# Set the current date for creating an output directory
DATE=`date "+%Y%m%d"`
OUT_DIR=/path/to/sv_calls/cnvpytor/1.3.1/${SAMPLE_ID}/${DATE}

# Define the log file path
LOG=${OUT_DIR}/${SAMPLE_ID}_cnvpytor.log

# Define the Pytor file path
PYTOR_FILE=${OUT_DIR}/${SAMPLE_ID}.pytor

###################################################################

# Create the output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Create a temporary directory (if needed)
mkdir ${TMP_DIR}

# Check if the input BAM file exists
test -e ${BAM} || die "Sample ${SAMPLE_ID} does not have an alignment file. Exiting."

# Run cnvpytor with different options
( /usr/bin/time -v /home/913/ar5564/.local/bin/cnvpytor -root ${PYTOR_FILE} -rd ${BAM} || die "cnvpytor failed. Exiting." ) 2>&1 | tee -a ${LOG}  

( /usr/bin/time -v /home/913/ar5564/.local/bin/cnvpytor -root ${PYTOR_FILE} -his 1000 10000 100000 || die "cnvpytor failed. Exiting." ) 2>&1 | tee -a ${LOG}

( /usr/bin/time -v /home/913/ar5564/.local/bin/cnvpytor -root ${PYTOR_FILE} -partition 1000 10000 100000 || die "cnvpytor failed. Exiting." ) 2>&1 | tee -a ${LOG}

# Set the bin size
BIN=10000

# Generate CNV calls and save to a TSV file
/home/913/ar5564/.local/bin/cnvpytor -root ${PYTOR_FILE} -call ${BIN} > ${OUT_DIR}/${SAMPLE_ID}.calls.${BIN}.tsv

# Navigate to the output directory
cd ${OUT_DIR}

# Export data for JBrowse visualization
/home/913/ar5564/.local/bin/cnvpytor -root ${SAMPLE_ID}.pytor -export jbrowse

# Convert bigWig to Wig format
bigWigToWig jbrowse_${SAMPLE_ID}/bw/${SAMPLE_ID}/his_rd_p_${BIN}.bw jbrowse_${SAMPLE_ID}/bw/${SAMPLE_ID}/his_rd_p_${BIN}.wig

# Convert Wig to Bed format and compress it
wig2bed --do-not-sort < jbrowse_${SAMPLE_ID}/bw/${SAMPLE_ID}/his_rd_p_${BIN}.wig | gzip > ${SAMPLE_ID}.bed.gz

