#!/bin/bash
# ------------------------------------------------------------------
# [Author] Title
#          Description
# ------------------------------------------------------------------

cut -f1,2 gtdbtk.ar53.summary.tsv | sed 's/;/\t/g' | sed 's/s__\b\|g__\b\|f__\b\|o__\b\|c__\b\|p__\b\|d__\b//g' | sed 1,2d > MAGs_taxonomy_ARCH.txt
cut -f1,2 gtdbtk.bac120.summary.tsv | sed 's/;/\t/g' | sed 's/s__\b\|g__\b\|f__\b\|o__\b\|c__\b\|p__\b\|d__\b//g' > MAGs_taxonomy_BACT.txt
cat MAGs_taxonomy_BACT.txt MAGs_taxonomy_ARCH.txt > MAGs_FULL_taxonomy.txt

