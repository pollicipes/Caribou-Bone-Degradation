#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 3
#SBATCH --mem=3G
#SBATCH --time=01:00:00 
#SBATCH --array=1-115%10

set -Eeuo pipefail;

mkdir -p ./log;
module purge;
module load anvio/8.0;

# Modify the path accordingly
list="/home/ngp704/data/caribou/anvio_dirs/runs/REDUNDANT_MAGs/CONTIGS_DB/list.contigs" 
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

path="/home/ngp704/data/caribou/anvio_dirs/runs/bone_degradome_collagenases/"
contigs="/home/ngp704/data/caribou/anvio_dirs/runs/REDUNDANT_MAGs/CONTIGS_DB/"

# To annotate our 3 HMMs 
anvi-run-hmms -c ${contigs}/${sample} -H ${path}/HMM_PF01752 --add-to-functions-table -T 4
anvi-run-hmms -c ${contigs}/${sample} -H ${path}/HMM_PF00082 --add-to-functions-table -T 4
anvi-run-hmms -c ${contigs}/${sample} -H ${path}/HMM_PF16499 --add-to-functions-table -T 4

# To annotate our MAGs with the COGs, we need forst to add them to the CONTIGS.db
anvi-run-ncbi-cogs -c ${contigs}/${sample} \
		   --cog-data-dir /projects/mjolnir1/data/databases/anvio-8.x/ncbi-cogs/20231006/ \
		   -T 4
