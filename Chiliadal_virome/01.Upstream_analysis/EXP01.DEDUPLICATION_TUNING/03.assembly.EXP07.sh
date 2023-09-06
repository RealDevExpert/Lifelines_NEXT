#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=48gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load BBMap

mkdir -p ../${SAMPLE_ID}/assembly/EXP07

module load SPAdes/3.15.3-GCC-11.2.0

### Assembly of paired-processed reads only
spades.py \
        --sc \
        --only-assembler \
        -k 21,33,55,77,99,127 \
	-1 ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_1.fastq \
	-2 ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_2.fastq \
	-s ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_singletons.fastq \
        -o ../${SAMPLE_ID}/assembly/EXP07 \
	-m $((${SLURM_MEM_PER_NODE} / 1000)) \
        --threads ${SLURM_CPUS_PER_TASK}


module list

module purge

