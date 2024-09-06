#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=03:00:00 

set -Eeuo pipefail;

mkdir -p ./log;
module purge;
module load anvio/8.0;

path="/home/ngp704/data/caribou/anvio_dirs/runs/metabolism/";
ext="/home/ngp704/data/caribou/anvio_dirs/runs/external_genomes/external_genomes.txt";
modules="/home/ngp704/data/caribou/anvio_dirs/runs/bone_degradome_collagenases/BONE_DEGRADATION_METABOLISM/";

# Annotate the bone degradome  
anvi-estimate-metabolism -e ${ext} \
			 --user-modules ${modules} \
			 --only-user-modules \
			 -O $path"/metabolism_bone_degradome" \
			 --matrix-format;
