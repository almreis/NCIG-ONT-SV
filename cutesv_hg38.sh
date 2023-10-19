#!/bin/bash
#PBS -N CUTESV 
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=128GB
#PBS -l walltime=16:00:00
#PBS -l wd

###################################################################

module load cuteSV/1.0.13

# Define BAM file path
BAM=$(ls -1tr /path/to/data/alignments/minimap2/v2.22-r1101/hg38/${SAMPLE_ID}/*/${SAMPLE_ID}_pass.bam | tail -n1)

num_threads=${PBS_NCPUS}

# Define temporary directory and log file
TMP_DIR=/path/to/${SAMPLE_ID}_${ALIGNER}_cutesv
TMP_LOG=${TMP_DIR}/sv_call.${SAMPLE_ID}_pass.vcf.log

CUTESV_VERSION=`cuteSV --version | awk '{print "v"$2}'`
DATE=`date "+%Y%m%d"`
OUT_DIR=/path/to/sv_calls/cutesv/${CUTESV_VERSION}/hg38/${SAMPLE_ID}/${DATE}

# Define VCF file and main log file
VCF=${OUT_DIR}/${SAMPLE_ID}_pass.vcf
LOG=${OUT_DIR}/sv_call.${SAMPLE_ID}_pass.vcf.log

# Define reference paths
FEMALE_REF=/path/to/references/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xo.no_alt.no_decoy.fa
MALE_REF=/path/to/references/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.no_alt.no_decoy.fa
REFERENCE=/path/to/references/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.no_alt.no_decoy.fa

###################################################################

mkdir -p ${OUT_DIR}

mkdir ${TMP_DIR}

test -e ${BAM} || die "Sample ${SAMPLE_ID} does not have an alignment file. Exiting."

# Run cuteSV with specified parameters
( /usr/bin/time -v cuteSV --threads ${num_threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --sample ${SAMPLE_ID} --report_readid --min_support 5 --min_size 20 --max_size 1000000 --genotype ${BAM} ${REFERENCE} ${VCF} ${TMP_DIR} || die "cuteSV failed. Exiting." ) 2>&1 | tee -a ${TMP_LOG}  

echo "Completed successfully" | tee -a ${TMP_LOG}
mv ${TMP_LOG} ${LOG} || die "moving ${TMP_LOG} to ${LOG} failed. Exiting." 
rm -r ${TMP_DIR}

