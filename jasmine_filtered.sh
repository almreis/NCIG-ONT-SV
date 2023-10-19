#!/bin/bash
#PBS -N JASMINE_FILTERED_CHROM
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=30:00:00
#PBS -l wd

###################################################################

module load jasminesv/1.1.4
module load minimap2/2.22
module load samtools/1.12

# Define paths to files and directories
GENOME=/path/to/chm13_reference_dir/chm13_v1.1_with_HG002_chrY.masked_PARs.fasta
VCF_LIST=/path/to/sv_calls/cutesv/v1.0.13/vcf_list.filtered.txt
BAM_LIST=/path/to/alignments/minimap2/v2.22-r1101/bam_list.txt
OUT_DIR=/path/to/sv_calls/cutesv/v1.0.13/all.filtered.joint_call.v2/${CHROM}
WORKING_DIR=${OUT_DIR}/jasmine_workdir_filtered
OUT=${OUT_DIR}/${CHROM}.all.filtered.joint_call.vcf
LOG=${OUT_DIR}/${CHROM}.all.filtered.joint_call.log
THREADS=${PBS_NCPUS}

mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

# Define per-chromosome VCF and BAM lists
VCF_LIST_CHROM=/path/to/sv_calls/cutesv/v1.0.13/all.filtered.joint_call.v2/${CHROM}/${CHROM}.vcf_list.filtered.txt
rm -f ${VCF_LIST_CHROM}
while read VCF
do
	SAMPLE_ID=$(echo ${VCF} | tr "/" "\t" | cut -f10)
	VCF_CHROM=${OUT_DIR}/${SAMPLE_ID}.${CHROM}.vcf
	grep "^#" ${VCF} > ${VCF_CHROM}
	grep -v  "^#" ${VCF} | grep -w "^${CHROM}" | grep -v "SVTYPE=BND" >> ${VCF_CHROM} 
	realpath ${VCF_CHROM} >> ${VCF_LIST_CHROM}
done < ${VCF_LIST}

BAM_LIST_CHROM=/path/to/sv_calls/cutesv/v1.0.13/all.filtered.joint_call.v2/${CHROM}/${CHROM}.bam_list.filtered.txt
rm -f ${BAM_LIST_CHROM}
while read BAM
do
	SAMPLE_ID=$(echo ${BAM} | tr "/" "\t" | cut -f10)
	BAM_CHROM=${OUT_DIR}/${SAMPLE_ID}.${CHROM}.bam
	samtools view -b ${BAM} ${CHROM} > ${BAM_CHROM}
	samtools index ${BAM_CHROM}
	realpath ${BAM_CHROM} >> ${BAM_LIST_CHROM}
done < ${BAM_LIST}

echo "Run Jasmine" | tee -a ${LOG}
( /usr/bin/time -v jasmine threads=${THREADS} out_dir=${WORKING_DIR} genome_file=${GENOME} file_list=${VCF_LIST_CHROM} bam_list=${BAM_LIST_CHROM} out_file=${OUT} min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --allow_intrasample --clique_merging --dup_to_ins --normalize_type --run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants || die "Jasmine failed. Exiting." ) 2>&1 | tee -a ${LOG}

