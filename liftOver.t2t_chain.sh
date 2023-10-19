#!/bin/bash
#PBS -N LIFTOVER 
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l jobfs=300GB
#PBS -l walltime=10:00:00
#PBS -l wd

###################################################################

# T2T chain
CHAIN=/path/to/references/chm13_reference_dir/liftover_files/chm13v2-grch38.chain

/path/to/scripts/liftOver ${BED} ${CHAIN} ${LIFTED} ${UNMAPPED}
