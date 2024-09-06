#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 10
#SBATCH --mem=25G
#SBATCH --time=05:00:00 

set -Eeuo pipefail;

mkdir -p ./log;
module purge;

module load anvio/8.0;

anvi-dereplicate-genomes \
    -f /home/ngp704/data/caribou/anvio_dirs/runs/REDUNDANT_MAGs/fasta_paths.txt \
    --program fastANI \
    -o /home/ngp704/data/caribou/anvio_dirs/runs/dereplicate \
    --similarity-threshold 0.90 \
    -T 10

### USED MANUALLY TO DERREPLICATE AND CHOOSE THE GOOD GENOMES
# while read i;
# do
#     x=$(echo $i | awk '{print $4}') ;
#     cac=$(sed 's/,/\\|/g' <(echo ${x}));
#     echo ${i};
#     grep $cac ../REDUNDANT_MAGs/ALL_REDUNDANT_MAGs.stats;
# done < to_derreplicate_manually.txt > to_derreplicate_manually_wStats
# 
# awk '$2 > 1' CLUSTER_REPORT.txt  > to_derreplicate_manually.txt

### ADD ALL FINAL MAGS TOGETHER
# awk '$2 == 1' CLUSTER_REPORT.txt | awk '{print $3}' > ../FINAL_MAGs/FINAL_MAGs.list
# grep 'XXX' to_derreplicate_manually_wStats | cut -f1 | sed 's/-CONTIGS//g' >> ../FINAL_MAGs/FINAL_MAGs.list
