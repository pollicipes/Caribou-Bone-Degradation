###
# THIS IS THE ANALYSIS DONE FROM THE ALPHA/BETA DIVERSITY FROM THE MAGS
###

### CLUSTER VERSION ###

message("READ ABOVE!!")

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
library(vegan)

# Options:
options(scipen = 10)
options(rstudio.help.showDataPreview = FALSE)

# 1. READ COVERAGE FILE FOR THE MAGs ####
setwd("~/HPC_data/caribou/")
df = fread("anvio_dirs/runs/ALL_FINAL_MAGs_MERGED-SUMMARY/bins_across_samples/mean_coverage.txt", header=T)
df = as.matrix(df)
rownames(df) = df[,1]
mx1 = df[,-c(1)]
class(mx1) = "numeric"

# 2. METADATA ####
mtd = data.frame(fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/SAMPLES_metadata.txt", dec = ","))
# Tidy up md
mtd = as.matrix(mtd)
rownames(mtd) = mtd[,1]
mtd = as.data.frame(mtd)
mtd$year_date = factor(mtd$year_date, levels = c("year_1978","year_2021"))
mtd$layer = factor(mtd$layer, levels = unique(mtd$layer))
mtd$profile = factor(mtd$profile, levels = unique(mtd$profile))

# 3. BINS ####
mag_mtd = as.data.frame(fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/MAGs_graph/MAGs_metadata_F.txt", fill=T,na.strings = c("",NA),dec="."))
# genome size
gsize = data.frame(MAG = mag_mtd$derrep_MAGs, Length = mag_mtd$Length)

# 4. TAXONOMY ####
tax = mag_mtd[,c(6:12)]
# tax = fread("MAG_Tax_ALL.txt")
tax = as.data.frame(tax)
# class(tax)
# Tidy up tax
tax = as.matrix(tax)
rownames(tax) = mag_mtd$derrep_MAGs
tax = tax[,-c(1)]

# 5. NORMALIZATION TPM ####
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
# colSums(mx2)

# 6. ALPHA DIVERSITY ####
mx.counts = round(mx2*1000000) # Multiply the data TPM normalized by 1M to get back the counts
p.rich = phyloseq(otu_table(as.matrix(mx.counts), taxa_are_rows = TRUE), 
              tax_table(tax), 
              sample_data(mtd))
# estimate_richness(p.rich)
palpha_year = plot_richness(p.rich, x="year_date", measures="Shannon", color = "year_date") +
  geom_boxplot(alpha = 0.8, size = 1) + 
  scale_color_manual(breaks = c("year_1978","year_2021"), values = c("darkgray","orange")) +
  scale_x_discrete(limits=c("year_1978","year_2021"),labels=c("1978","2021")) +
  stat_compare_means(label.x = 1.35) + 
  xlab("Sample") +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  guides(color="none")
palpha_year
# ggsave(filename = "~/Projects/caribou/FIGURES/Shannon_alpha_YEAR_MAGs.pdf",palpha_year)

## 6.1 Alpha div per layer ####

# Reorder layers to plot x-axis on reverse.
sample_data(p.rich)$layer = factor(sample_data(p.rich)$layer, 
                               levels = c("Layer_4","Layer_3","Layer_2","Layer_1"),labels = c("Layer 4","Layer 3","Layer 2","Layer 1"))


palpha_layer = plot_richness(p.rich, x = "layer",   measures ="Shannon", shape = NULL) +
  geom_dotplot(aes(fill = year_date), binaxis='y',stackdir = "center",binwidth = 1/150,dotsize = 5) +
  #geom_boxplot(alpha =  0.8, size = 1) +
  #geom_point(size = 3) +
  #geom_jitter(aes(group=date)) +
  stat_summary(fun=mean, geom="line",lty = "dashed", lwd = 1.2, aes(group = year_date, color = year_date)) +
  scale_colour_manual(values = c("darkgray","orange")) + 
  xlab("Layer") +
  stat_compare_means(aes(group=year_date),label.y = c(1.95,1.95,1.95,1.95),vjust = -1.5) +
  labs(fill = "Year") +
  guides(color="none") +
  scale_fill_manual(labels = c("1978","2021"),breaks = c("year_1978","year_2021"),values = c("darkgray","orange")) +
  theme(text = element_text(size = 20)) +
  coord_flip()
palpha_layer

# png(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_LAYER_MAGs.png",
#     width = 19, height = 19, units ="cm", res = 600)
# palpha_layer
# ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Shannon_alpha_LAYER.pdf", plot = palpha_layer)
# dev.off()

# 7. BETA-DIVERSITY (TESTING) ####
set.seed(42)
dist.met = "bray";
meth1 ="NMDS" # "MDS" "PCoA" #  ;
ord.nmds.bray = ordinate(p.rich, method = meth1, distance = dist.met)

## NMDS Plot #####
plot_ordination(p.rich,
                ord.nmds.bray,
                color = "year_date",
                # color = "layer",
                # color = "profile",
                # shape="year_date"
                ) +
  geom_point(size = 4) + 
  geom_label_repel(aes(label = ID),size=3)

# 8. DIFFERENTIAL ABUNDANCE ####
library("DESeq2")

## 8.1 By DATE #####
# Filter 
pb1 = filter_taxa(p.rich, function(x) sum(x > 10000) > (0.10*length(x)),TRUE); pb1
# otu_table(pb1)

# otu_table(pb1)
das.date = phyloseq_to_deseq2(pb1, ~ year_date)
diagdd.date = DESeq(das.date, test = "Wald", fitType = "parametric")
res.date = results(diagdd.date, cooksCutoff = FALSE)
alpha = 1
sig.date = res.date[which(res.date$padj < alpha),]
sig.date$ID = rownames(sig.date)
diff.deseq.date = merge(as.data.frame(sig.date), mag_mtd, by.x="ID", by.y="derrep_MAGs")
# fwrite(diff.deseq.date,"~/Projects/caribou/DESEQ2_diff_abundance_by_date.txt", sep="\t", dec = ",")

## Plot pvalue x log-fold change ####
labeling.signif = diff.deseq.date %>% mutate(colv = ifelse(abs(log2FoldChange) > 2,"red3" ,"grey1")) %>% filter(abs(log2FoldChange) > 2)

diff_years = diff.deseq.date %>% mutate(colv = ifelse(abs(log2FoldChange) > 2,"red3" ,"grey1")) %>%
  ggplot(aes(x = -log10(padj),y = log2FoldChange,color=colv)) + 
  geom_label_repel(data = labeling.signif, aes(label = paste0(Phylum,"\n",Order,"\n",Family)), size = 4, min.segment.length = 0.4) +
  geom_point(size=2.5) + 
  geom_hline(yintercept=2,color="red3",lty="dashed") +
  geom_hline(yintercept=-2,color="red3",lty="dashed") +
  # geom_text_repel(aes(label = paste0(Order," | ",Family," | ",Genus)),size=3) +
  scale_color_identity() +
  geom_hline(yintercept = 0,col = "gray90",lty=2) +
  ylim(-8,8) +
  xlim(0,25) +
  theme(text = element_text(size = 20)) + 
  annotate("text", x = 2, y = 7, label = "Abundant in 2021",size = 10,fontface = "italic") +
  annotate("text", x = 2, y = -7, label = "Abundant in 1978", size = 10, fontface = "italic") 
  # ggtitle("Differential taxa between years\n logFold > 0: abundant in 2021; logFold < 0: abundant in 1978\n Phylum | Order | Family")
diff_years

#pdf(file = "~/Projects/caribou/FIGURES_PAPER/Differential_2021-1978_samples.pdf",width = 15, height = 10) # , units ="cm", res = 600)
#diff_years
#dev.off()

stop("FINISH HERE!")

# END  ####