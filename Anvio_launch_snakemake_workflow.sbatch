#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=10:00:00

set -Eeuo pipefail;

module purge;
module load anvio/7.1;

cnf=$1
anvi-run-workflow -w metagenomics -c ${cnf} # json file as an argument
