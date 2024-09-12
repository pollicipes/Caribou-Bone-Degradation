#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 10
#SBATCH --mem=60G
#SBATCH --time=60:00:00 

set -Eeuo pipefail;

mkdir -p ./log;

# Remember to activate the environment before submitting.

###########################################
#### NOTE: IT IS metaDMG2, not metaDMG ####
###########################################
# conda activate /projects/mjolnir1/people/tmk528/.conda/envs/metaDMG2;

cd /home/ngp704/data/caribou/results/metadmg/FINAL_MAGs_mapping;

metaDMG config \
	/home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs/*/*filtered.bam /home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs/2021/SAMPLE_101.MAGs_sorted.rmdup.bam /home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs/2021/SAMPLE_116.MAGs_sorted.rmdup.bam \
	--damage-mode local \
	-m /home/tmk528/data/bin/metaDMG-cpp/metaDMG-cpp \
	--cores-per-sample 5 \
	--parallel-samples 2;

metaDMG compute;
metaDMG convert --output results_as_csv.csv;

# Interactive view in browser:
# ssh -L 8050:localhost:8050 ngp704@mjolnirhead01fl.unicph.domain
# metaDMG dashboard config.yaml --port 8050
