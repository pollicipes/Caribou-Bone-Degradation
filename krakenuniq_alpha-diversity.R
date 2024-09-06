###
# 0. SCRIPT USED TO ANALYZE THE RESULTS FROM KRAKENUNIQ ####
###

library(tidyverse)
library(phyloseq)
library(ape)
library(data.table)
library(microbiome)
library(RColorBrewer)
library(scales)
library(viridis)
library(microViz)
library(ggpubr)
library(ggrepel)

options(scipen = 10)
setwd("~/HPC_data/caribou/results/krakenuniq/")
ind = "counts" # "abund"

# 0. SETUP ####
# OTU count table
level = "Genus" # "Genus" # "Species"
# Read abundances bracken 
rawd = fread(paste0("bracken_abundances_",level,"_ALL")) %>% as.data.frame()
if(ind == "counts"){
  df0 = rawd[,c(1,seq(4,ncol(rawd),2))]
} else {
  df0 = rawd[,c(1,seq(5,ncol(rawd),2))]
}
sx0 = gsub("SAMPLE_|\\.", "", colnames(df0))
sx1 = gsub("\\D.*", "", sx0)
df1 = df0[,-1]
colnames(df1) = c(sx1[-1])

# Metadata read 
mtd = fread("metadata_layer4.txt") # Use the metadata.txt if we want the original metadata and not layers 4 and 4B merged.
mtd[] = lapply(mtd, function(x) gsub(",", ".", x))
rownames(mtd) = mtd$ID
# convert to numeric type all numbers.
mtd[] = lapply(mtd, type.convert, as.is = TRUE)

# Adjust DNA % bins:
# dna_thr = 9
# mtd = mtd %>% mutate(classif_caribou_dna = case_when(caribou_dna <= dna_thr ~ "Low", caribou_dna > dna_thr ~ "Optimal"))

mtd$date = factor(mtd$date,levels=c("1978","2021"))
# mtd$layer = factor(mtd$layer,levels=unique(mtd$layer))
# mtd$profile = factor(mtd$profile,levels=unique(mtd$profile))

# Taxonomy read
taxonomy_tab = as.matrix(df0$name)
colnames(taxonomy_tab) = level
# assign rownames to the tax level chosen
rownames(taxonomy_tab) = taxonomy_tab[,1] # [,tax_lvl_int]
rownames(df1) = rownames(taxonomy_tab)

# CLEAN CROSSCHECK NAMES
inters_names = intersect(colnames(df1), rownames(mtd))
otu_cleaned = df1[, inters_names]
metadata_cleaned = mtd[match(inters_names, rownames(mtd)),]
mx1 = sample_data(metadata_cleaned, errorIfNULL = FALSE)
mx1$caribou_dna = as.numeric(mx1$caribou_dna)
sample_names(mx1) = metadata_cleaned$ID

# Read into phyloseq
p1 = phyloseq(otu_table(as.matrix(otu_cleaned),taxa_are_rows = TRUE),
              tax_table(taxonomy_tab),mx1)
dim(otu_table(p1))

# ALPHA DIVERSITY BY DATE ####
# Correct by median sequencing depth
total = median(sample_sums(p1))
standf = function(x, t = total) round(t * (x / sum(x))) 
pad = transform_sample_counts(p1, standf)
ad = prune_taxa(taxa_names(pad) != "Homo sapiens", pad) # Remove human comtamination
rich = estimate_richness(ad)

# Plots
theme_set(theme_bw())
palpha_year = plot_richness(ad, x="date", measures="Shannon", color = "date") +
  geom_boxplot(alpha = 0.8,size=1) +
  scale_color_manual(breaks = c("1978","2021"), values = c("darkgray","orange")) +
  scale_x_discrete(limits=c("1978","2021"),labels=c("1978","2021")) +
  stat_compare_means(label.x = 1.35) + 
  xlab("Year") +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  guides(color="none")

palpha_year

# png(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_YEAR_KrakenUniq.png",
#    width = 19, height = 19, units ="cm", res = 600)
# palpha_year
# dev.off()
# ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_YEAR_KrakenUniq.pdf",palpha_year)

# ALPHA DIV PER LAYER
# Reorder layers to plot x-axis on reverse.
sample_data(ad)$layer = factor(sample_data(ad)$layer, 
                               levels = c("Layer 4","Layer 3","Layer 2","Layer 1"))

palpha_layer = plot_richness(ad, x = "layer", measures = "Shannon", shape = NULL) +
  geom_dotplot(aes(fill = date), binaxis='y', stackdir='center',dotsize = 0.5) +
  #geom_boxplot(alpha =  0.8, size = 1) +
  #geom_point(size = 3) +
  #geom_jitter(aes(group=date)) +
  stat_summary(fun=mean, geom="line",lty = "dashed", lwd = 1.2, aes(group = date, color = date)) +
  scale_colour_manual(values = c("darkgray","orange")) + 
  xlab("Layer") +
  # stat_compare_means(aes(group=date),label.y = c(3.5,3.5,3.5,3.5),vjust = -1.5) +
  labs(fill = "Year") +
  guides(color="none") +
  scale_fill_manual(labels = c("1978","2021"),breaks = c("1978","2021"),values = c("darkgray","orange")) +
  theme(text = element_text(size = 20)) +
  coord_flip()

palpha_layer

# png(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_LAYER_KrakenUniq.png",
#     width = 19, height = 19, units ="cm", res = 600)
# palpha_layer
# dev.off()
# 
# ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_LAYER_KrakenUniq.pdf",palpha_layer)

# END #### 