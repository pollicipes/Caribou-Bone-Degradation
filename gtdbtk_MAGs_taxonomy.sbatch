#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 20
#SBATCH --mem=150G
#SBATCH --time=30:00:00 

set -Eeuo pipefail;

mkdir -p ./log;
module purge;
module load gtdbtk/2.2.6

cd "/home/ngp704/data/caribou/anvio_dirs/runs/";
GTDBTK_DATA_PATH="/projects/mjolnir1/data/databases/GTDBTK_DB/release214/";

gtdbtk classify_wf \
       # FOLDER WITH THE MAGS IN FASTA
       --genome_dir /home/ngp704/data/caribou/anvio_dirs/runs/REDUNDANT_MAGs/ \ 
       --out_dir /home/ngp704/data/caribou/anvio_dirs/runs/REDUNDANT_MAGs_gtdbtk \
       --cpus 20 \
       --mash_db gtdb-tk-r214.msh \
       --extension "fa";
