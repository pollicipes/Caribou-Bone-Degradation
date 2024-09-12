#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=10:00:00 
#SBATCH --array=1-1%1

set -Eeuo pipefail;

module purge;
module load bowtie2 samtools/1.6;
year=$1;
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ~/data/caribou/fastq/${year}/fastq.list | awk '{print $1}');

OUT="/home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs/"${year}"/"
mkdir -p ${OUT};
DATA="/home/ngp704/data/caribou/anvio_dirs/runs/"${year}"/01_QC/"
cd ${DATA};
DB="/home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs/MAGs_concatenated_for_metaDMG/ALL_MAGs_single_contig_per_genome"

bowtie2 --threads 8 \
	-x ${DB} \
        -1 ${sample}"-FILTERED_R1.fastq.gz" \
        -2 ${sample}"-FILTERED_R2.fastq.gz" \
        --no-unal | samtools view -@ 8 -b -F 4 - | samtools sort -@ 8 -O bam - > ${OUT}"/"${sample}".MAGs_sorted.bam"

module purge;
module load jdk picard;
cd ${OUT};

picard MarkDuplicates I=${sample}".MAGs_sorted.bam" O=${sample}".MAGs_sorted.rmdup.bam" M=${sample}"_picard.rmdup.log" AS=True REMOVE_DUPLICATES=True;

exit;
