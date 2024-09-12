#!/bin/bash

#SBATCH --output=./log/%x_o.%A
#SBATCH --error=./log/%x_e.%A
#SBATCH -c 5
#SBATCH --mem=8G
#SBATCH --time=3:00:00 

set -Eeuo pipefail;

mkdir -p ./log;
module purge;

cd /home/ngp704/data/caribou/anvio_dirs/runs/FINAL_MAGs/bowtie2_indexed_MAGs;

# Iterate all the MAGs and concatenate all the contigs in a single contig.
while read l; do
    # get the mag name
    mag=`basename $l`;
    # the number of contigs
    num=$(grep '>' ${l} | wc -l);
    echo $mag $num;
    # merge the $num number of contigs into a single 1. 
    awk -v i=$num -f merge_contigs.awk ${l} | sed 1d > MAGs_concatenated_for_metaDMG/${mag};
    # make a header for the sequence
    h=${mag%.fa};
    # add the sequence to the fasta file
    sed  -i "1i >$h" MAGs_concatenated_for_metaDMG/${mag};
# Read the MAGs in fasta
done <  <(ls ../../FINAL_MAGs/FASTAs/*fa)


module load bowtie2;
bowtie2-build ALL_MAGs_single_contig_per_genome.fa ALL_MAGs_single_contig_per_genome;
