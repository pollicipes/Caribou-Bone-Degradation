#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=20        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=350G
#SBATCH --time=120:00:00          # total run time limit (HH:MM:SS)
#SBATCH --array=1-17%17

module load krakenuniq/1.0.3;
set -Eeuo pipefail;

year=$1;
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ~/data/caribou/fastq/${year}/fastq.list | awk '{print $1}');

RAW="/home/ngp704/data/caribou/fastq/"${year}"/fastp_filtered_data_merged_collapsed/"
DB="/datasets/globe_databases/krakenuniq/standard/20200816/2020-08-16";
tipo="standard";
OUT="/home/ngp704/data/caribou/results/krakenuniq_collapsed/"${year}"/"

krakenuniq --report-file ${OUT}"/"${sample}"_kraken_report_"${tipo}"_collapsed" --threads 20 --db ${DB} ${RAW}"/"${sample}"_collapsed+forwardread.fq.gz"

module load bracken/2.8;
bracken-build 
bracken -d $DB -i ${OUT}"/"${sample}"_kraken_report_MT_SE_dedup" -o ${OUT}"/"${sample}"_kraken_report_MT_SE_dedup.abundance"  -w ${OUT}"/"${sample}"_kraken_report_MT_SE_dedup_bracken.out" -l G # -r  <read length>

