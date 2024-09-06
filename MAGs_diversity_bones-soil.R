###
## THIS IS THE ANALYSIS DONE FROM THE ALPHA/BETA DIVERSITY FOR THE SOIL SAMPLES + 1978 + 2021 batch using the MAGs
###

message("READ ABOVE!!")

library(tidyverse)
library(RColorBrewer)
library(microbiome)
library(data.table)
library(patchwork)
library(phyloseq)
library(microViz)
library(viridis)
library(ggrepel)
library(scales)
library(ggpubr)
library(vegan)
library(ape)
theme_set(theme_classic())

# Options:
options(scipen = 10)
options(rstudio.help.showDataPreview = FALSE)

# 1. READ AND CLEAN COVERAGE FILE FOR THE MAGs AND SOIL ####
setwd("~/HPC_data/caribou/") # "results/alpha_beta_diversity_analysis_MAGs/"
cvg_years = "anvio_dirs/runs/soil/summary_mean_coverage_1978_2021.txt"
cvg_soil = "anvio_dirs/runs/soil/summary_mean_coverage_soil.txt"

df_years = fread(cvg_years, header=T)
df_soil = fread(cvg_soil, header=T,drop = "bins")
df = as.matrix(cbind(df_years,df_soil))
rownames(df) = df[,1]
mx1 = df[,-c(1)]
class(mx1) = "numeric"

# 2. METADATA ####
mtd = data.frame(fread("/Users/jrodriguez/Projects/caribou/Sediment_samples/soil_mtd_NEW.txt", dec = ","))
# Tidy up md
mtd = as.matrix(mtd)
rownames(mtd) = mtd[,1]
mtd = as.data.frame(mtd)
mtd$year_date = factor(mtd$year_date, levels = c("year_1978","year_2021","soil"))
mtd$layer = factor(mtd$layer, levels = unique(mtd$layer))

# 3. TAXONOMY ####

mag_mtd = as.data.frame(fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/MAGs_graph/MAGs_metadata_F.txt", fill=T,na.strings = c("",NA),dec="."))
# genome size
gsize = data.frame(MAG = mag_mtd$derrep_MAGs, Length = mag_mtd$Length)

tax = mag_mtd[,c(6:12)]
tax = as.data.frame(tax) # class(tax)
# Tidy up tax
tax = as.matrix(tax)
rownames(tax) = mag_mtd$derrep_MAGs
tax = tax[,-c(1)]

# 4. NORMALIZATION TPM ####
NormalizeTPM <- function(data, sce = NULL, tr_length = NULL, log = FALSE, scale = 1e+06, pseudo.count = 1) {
  
  # Median length of all transcripts for a given gene
  med_length <- stats::aggregate(x = tr_length$Length,
                                 by = list(hgnc_symbol = gsize$MAG), FUN = stats::median)
  common <- intersect(rownames(data), med_length$hgnc_symbol)
  if (length(common) < 2)
    stop("No enough overlap between data rownames and transcript length.\n")
  data <- data[common, ]
  med_length <- med_length[match(common, med_length$hgnc_symbol), ]
  # divide by length of transcript in kb - reads per kilobase (RPK)
  data1 <- sweep(data, 1, STATS = med_length$x/1000, FUN = "/")
  # divide by million RPK in sample - transcripts per million (TPM)
  data2 <- sweep(data1, 2, STATS = colSums(data1)/(10^6), FUN = "/")
  data3 <- data2/scale
  if (log) {
    if (pseudo.count == 0)
      warning("Using 0 pseudocount: Inf may be generated.\n")
    data <- log2(data + pseudo.count)
  }
  if(is.null(sce)){ return(data3) } else{ return(sce) }
}
mx2 = NormalizeTPM(data = mx1, tr_length = gsize, scale = 1e+06)

# 5. !!! THIS SHOULD BE APPLIED IF WE WANT ONLY THE DIFFERENTIAL ABUNDANCE SOIL-2021 COMPARISON ####
# Uncomment this and go to section 8.
## s78 = mtd[mtd$year_date == "year_1978",]$ID
## mx2 = mx2[,!colnames(mx2) %in% s78]

# 6. ALPHA DIVERSITY ####

mx.counts = round(mx2*1000000) # Multiply the data TPM normalized by 1M to get back the counts
colSums(mx.counts)

# Read to phyloseq 
p.rich.0 = phyloseq(otu_table(as.matrix(mx.counts), taxa_are_rows = TRUE), 
                  tax_table(tax), 
                  sample_data(mtd))

# Plot alpha diversity
palpha_year = plot_richness(p.rich.0, x = "year_date", measures = "Shannon", color = "year_date") +
  scale_color_manual(labels = c("1978","2021","soil"),breaks = c("year_1978","year_2021","soil"),values = c("darkgray","orange","blue")) +
  geom_boxplot(alpha = 0.6) +
  guides(color = "none") +
  xlab("Year") +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  stat_pwc(label.size = 6, label = "Wilcoxon, p = {p}")

palpha_year
# ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_YEAR_SOIL_MAGs.pdf",palpha_year)

# Alpha div per layer
# Reorder layers to plot x-axis on reverse.
sample_data(p.rich.0)$layer = factor(sample_data(p.rich.0)$layer, levels = c("Layer_4","Layer_3","Layer_2","Layer_1"))

palpha_layer = plot_richness(p.rich.0, x = "layer",   measures ="Shannon", shape = NULL) +
  geom_dotplot(aes(fill = year_date), binaxis='y',stackdir = "center",binwidth = 1/150,dotsize = 5) +
  stat_summary(fun=mean, geom="line",lty = "dashed", lwd = 1.2, aes(group = year_date, color = year_date)) +
  scale_colour_manual(values = c("darkgray","orange","blue")) + 
  xlab("Layer") +
  stat_compare_means(aes(group=year_date),label.y = c(2,2,2,2),vjust = -1.5) +
  labs(fill = "Year") +
  guides(color="none") +
  scale_fill_manual(labels = c("1978","2021","Soil"),breaks = c("year_1978","year_2021","soil"),values = c("darkgray","orange","blue")) +
  theme(text = element_text(size = 20)) +
  coord_flip()

palpha_layer

# png(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_LAYER_SOIL_MAGs.png",
#    width = 19, height = 19, units ="cm", res = 600)
# palpha_layer
# dev.off()

# 7. B-DIVERSITY ####
p.rich = microViz::tax_fix(p.rich.0) # COOL TO FIX TAXA TABLE NAMES, IF NEEDED!

# Filter 
pb1 = filter_taxa(p.rich, function(x) sum(x > 1000) > (0.10*length(x)),TRUE); pb1
# (otu_table(pb1))

dna1 = pb1 %>% 
  tax_transform(rank = "Genus", trans = "normalize") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc(method = "auto") %>%
  ord_plot(axes = c(1, 2),
           color = "year_date", 
           alpha = 0.5, size = 5) + 
  microViz::stat_chull(ggplot2::aes(colour = year_date)) +
  theme(text = element_text(size = 18)) +
  labs(color = "Origin") +
  scale_color_manual(labels = c("1978","2021","Soil"),breaks = c("year_1978", "year_2021", "soil"),values = c("gray34","orange","blue")) +
  geom_text_repel(aes(label = gsub("SAMPLE_","",ID)), size = 3,min.segment.length = 0.2)

dna1
#ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/MDS_with_NEW_soil_samples.pdf",dna1)
#ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/MDS_with_2021-SOIL_samples.pdf",dna1)


## 7.1 EUCLIDEAN DISTANCE between soil and each year. ####

df = data.frame(dna1$data$MDS1, dna1$data$MDS2, dna1$data$year_date) 
colnames(df) = c("MDS1","MDS2","year_date")
df_2021 <- df[df$year_date == "year_2021", c("MDS1", "MDS2")]
df_1978 <- df[df$year_date == "year_1978", c("MDS1", "MDS2")]
df_soil <- df[df$year_date == "soil", c("MDS1", "MDS2")]
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}
calculate_distances <- function(soil_point, year_points) {
  apply(year_points, 1, function(year_point) {
    euclidean_distance(soil_point[1], soil_point[2], year_point[1], year_point[2])
  })
}
distances_2021 <- apply(df_soil, 1, function(soil_point) calculate_distances(soil_point, df_2021))
distances_1978 <- apply(df_soil, 1, function(soil_point) calculate_distances(soil_point, df_1978))
avg_distances_2021 <- rowMeans(distances_2021)
avg_distances_1978 <- rowMeans(distances_1978)
# closer_to <- ifelse(avg_distances_2021 < avg_distances_1978, "year_2021", "year_1978")
overall_avg_2021 <- mean(avg_distances_2021)
overall_avg_1978 <- mean(avg_distances_1978)
print(paste("Average distance to 2021:", overall_avg_2021))
print(paste("Average distance to 1978:", overall_avg_1978))

if (overall_avg_2021 < overall_avg_1978) {
  print("On average, soil points are closer to year 2021")
} else {
  print("On average, soil points are closer to year 1978")
}

# 8. DESeq2 Differential between the 2021 and the soil samples. ####

# Filter 
pb1 = filter_taxa(p.rich, function(x) sum(x > 1000) > (0.10*length(x)),TRUE); pb1
# (otu_table(pb1))

library(DESeq2)
das.date = phyloseq_to_deseq2(pb1, ~ year_date)
diagdd.date = DESeq(das.date, test = "Wald", fitType = "parametric")
res.date = results(diagdd.date, cooksCutoff = FALSE)
alpha = 1
sig.date = res.date[which(res.date$padj < alpha),]
sig.date$ID = rownames(sig.date)
diff.deseq.date = merge(as.data.frame(sig.date), mag_mtd, by.x="ID", by.y="derrep_MAGs")

# Top-significant hits ####
head(diff.deseq.date[order(diff.deseq.date$log2FoldChange),])
# otu_table(pb1)

## Plot pvalue x log-fold change ####
labeling.signif = diff.deseq.date %>% mutate(colv = ifelse(abs(log2FoldChange) > 2,"coral1" ,"grey")) %>% filter(abs(log2FoldChange) > 2)
diff1 = diff.deseq.date %>% mutate(colv = ifelse(abs(log2FoldChange) > 2,"coral1" ,"grey")) %>%
  ggplot(aes(x = -log10(padj),y = log2FoldChange,color=colv)) + 
  geom_point(size = 2.5) + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=2,color="red",lty="dashed") +
  geom_hline(yintercept=-2,color="red",lty="dashed") +
  # geom_text_repel(aes(label = paste0(Order," | ",Family," | ",Genus)),size=3) +
  geom_text_repel(data = labeling.signif, aes(label = paste0(Phylum," | ",Order," | ",Family)), size = 4,min.segment.length = 0.1) +
  scale_color_identity() +
  geom_hline(yintercept = 0,col = "gray90",lty=2) +
  ylim(-8,8) +
  xlim(0,30) + 
  ggtitle("Differential taxa between samples\n logFold < 0: abundant in 2021; logFold > 0: abundant in the sediment\n Phylum | Order | Family")

diff1
# ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Differential_2021-SOIL_samples.pdf",diff1)

# END #### 
