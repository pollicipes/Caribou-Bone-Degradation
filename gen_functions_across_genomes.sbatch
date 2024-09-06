#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 6
#SBATCH --mem=6G
#SBATCH --time=01:00:00 

set -Eeuo pipefail;

mkdir -p ./log;
module purge;

module load anvio;

annot=$1 # This elements one at the time, as string arguments: HMM_PF00082, COG20_PATHWAY, COG20_CATEGORY, HMM_PF01752, COG20_FUNCTION, HMM_PF16499
int_gen="/home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/internal_genomes.txt"
out="/home/ngp704/data/caribou/anvio_dirs/runs/metabolism/functions_across_genomes"

anvi-script-gen-function-matrix-across-genomes -i ${int_gen} \
					       -O ${out}/${annot}_matrix_functions_across_genomes \
					       --annotation-source ${annot} \
					       --also-report-accession-ids
