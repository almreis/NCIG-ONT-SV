#!/bin/bash
#PBS -P te53
#PBS -N CENTRIFUGE_CHM13
#PBS -q normal
#PBS -l ncpus=8
#PBS -l mem=192GB
#PBS -l jobfs=100GB
#PBS -l walltime=24:00:00
#PBS -l wd
#PBS -l storage=gdata/te53

###################################################################

export CENTRIFUGE_HOME=/path/to/Figures/Review/microbial_contamination/centrifuge

INDEX_DIR=/path/to/Figures/Review/microbial_contamination/database/

FASTQ_PASS=`ls -tr /path/to/data/basecalls/*/*/${SAMPLE_ID}/*/*pass.fastq.gz | tail -n1`

num_threads=${PBS_NCPUS}

OUTDIR=$(ls -1dtr /path/to/data/alignments/minimap2/v2.22-r1101/${SAMPLE_ID}/* | tail -n1)

REPORT=${OUTDIR}/${SAMPLE_ID}.centrifuge_report.tsv

CLASSIFICATION=${OUTDIR}/${SAMPLE_ID}.centrifuge_classification.tsv

###################################################################

$CENTRIFUGE_HOME/centrifuge -p ${num_threads} -q -x ${INDEX_DIR}/p_compressed+h+v -U ${FASTQ_PASS} -S ${CLASSIFICATION} --report-file ${REPORT}

