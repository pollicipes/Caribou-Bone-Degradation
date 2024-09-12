# MetaDMG analysis; plot the results.
# Following the script from https://github.com/miwipe/KapCopenhagen/blob/main/scripts/KapKmetaDMG6AnimalFamilies.R

# Get rid of the rstudio error trying to show preview of data frames
options(rstudio.help.showDataPreview = FALSE)

library(tidyverse)
library(data.table)
library(gplots)
library(scales)
library(reshape2)
library(vegan)
library(rioja)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)

# 1. READ AND FILTER INITIALLY ####
a = fread("~/HPC_data/caribou/results/metadmg/FINAL_MAGs_mapping/results_as_csv.csv",sep=",")
# View(table((unique(a[,c(1,2)]))))

# Filters
# MAP significance
signif = 1.5
# Minimum reads for parsing taxa
MinRead = 5
# Minimum average mean read length
MinLength = 30
# Minimum amount of damage (larger or smaller)
# DamMin = 0.02  ### > 0.025 

## 1.1 EXTRACT THE READ COUNTS ####
dw1 = a %>% filter(
    # MAP_damage < DamMin, 
    N_reads > MinRead, 
    mean_L > MinLength, 
    MAP_significance  > signif) %>% 
  as.data.frame() %>% 
  reshape2::dcast(tax_id ~ sample, value.var="N_reads") %>% # cast the data frame into a matrix with the read number per sample X MAGs
  column_to_rownames(var = "tax_id") # convert first column to rownames, which contain the MAG names, to get a numbers-only table
dim(a); dim(dw1)
# set NAs as zeros
dw1[is.na(dw1)] = 0 

# 2. FILTER MAGs OUT FOR LOW COUNTS ####
# We need to do this, because otherwise, as we are analyzing damage, we won't know whether that MAG is not damaged or not present

# If we want to remove samples based on low coverage:
# col = median.default(colSums(dw1))/2 
# drop = colSums(dw1) > col # applies the threshold1
# dd1 = dw1[c(drop)] # applies the threshold1
# View(dd1)

# Removes all taxa with less than half the median of the sum of reads to all taxa 
row.median2 = median.default(rowSums(dw1))/2
dw2 = as.data.frame(t(dw1)) # tranposes the dataframe  
drop2 = colSums(dw2) > row.median2 # applies the threshold2
dw3 = dw2[c(drop2)] # applies the threshold2 dim(dw3)
dw4 = as.data.frame(t(dw3)) # transposes the dataframe (df) back to original
message("We have ",dim(dw4)[1]," MAGs and ",dim(dw4)[2]," SAMPLES remaining")

# adding new column with sum of reads
b1 = cbind(dw4, NReads = rowSums(dw4))
# Number of columns with observations (excludes sum column). Equal to sample number. Should be.
cdat = ncol(b1) - 1
# Setting replicated e.g. how many times must a taxa appear in different samples to be considered 'replicated', considers all samples as different. 
b2 = cbind(b1, Nreplicated = rowSums(b1[,1:cdat] > 0))
b2 = b2[b2$Nreplicated >= 3,] # applies the replication threshold summary(b2$Nreplicated)
message("We have ",dim(b2)[1]," MAGs and ",dim(b2)[2]-2," SAMPLES remaining")

# 3. GETS THE NAME OF FILTERED MAGS FOR DAMAGE ####
# We select those passing filters to be kept for damage analyses.
mags.sel = rownames(b2)

mag.df = a %>% filter(
  # MAP_damage < DamMin, 
  N_reads > MinRead, 
  mean_L > MinLength, 
  MAP_significance  > signif) %>% 
  as.data.frame() %>% 
  reshape2::dcast(tax_id ~ sample, value.var="MAP_damage") %>%
  filter(tax_id %in% mags.sel) %>% 
  column_to_rownames(var = "tax_id") # convert first column to rownames, which contain the MAG names, to get a numbers-only table

# set NAs as zeros
mag.df[is.na(mag.df)] = 0 

## 3.1 GRAPH VIEW ####
col = colorRampPalette((brewer.pal(9, "PuRd")))(20)
ComplexHeatmap::Heatmap(as.matrix(mag.df),col = col, border = T)
# 3D bar heatmap; with MAG ID.
ComplexHeatmap::Heatmap3D(as.matrix(mag.df), col = col)

# Do 3D heatmap with the taxa names
mag.df$mag = rownames(mag.df)
mag_mtd = as.data.frame(fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/MAGs_graph/MAGs_metadata_F.txt", fill=T,na.strings = c("",NA),dec="."))
mag.df.tax = merge(mag.df, mag_mtd, by.x="mag", by.y="derrep_MAGs")
mag.df.tax$id = seq(1, nrow(mag.df.tax))
n1 = paste0(mag.df.tax$Order,"_",mag.df.tax$Family,"_",mag.df.tax$id)
# Select the columns for the heatmap
mag.df.tax = mag.df.tax[,c(2:34)] # Column 42 is Order
rownames(mag.df.tax) = n1
ComplexHeatmap::Heatmap3D(as.matrix(mag.df.tax),col = col, bar_rel_width = 0.4, show_column_dend = T,bar_angle = 45)

# EXPORT THE DATA FRAME WITH THE DAMAGE PER SAMPLE AND MAG.
fwrite(mag.df, file = "~/Projects/caribou/anvio_dirs/metadmg/samples_damage_MAGs.txt", sep="\t",row.names = T, col.names = T, dec=",")

message("Now use the anvio interactive to plot/view the damage!")

###
stop("FIRST PART FOR THE PAPER ENDS HERE. IF WE WANT TO PROCEED TO THE PLOTS DONE IN THE ORIGINAL SCRIPT FROM THE KAP KOBENHAVN GITHUB (Kjaer et al 2022, Nature), WE CAN CONTINUE.")
### 



# 4. GET THE READ NUMBERS MAPPED TO MAGS. ####
i = ncol(b2) - 1
b3 = as.matrix(b2[,seq(1,i)])  # change number of columns to the total of my_data

# IT SHOULD BE NORMALIZED BY SAMPLE, AS THE % ABUNDANCE PER SAMPLE IS ALREADY PROPORTIONAL
b4 = prop.table(data.matrix(b3), margin = 2)*100 # makes proportion table, needs 2 margins e.g. header and 1st row names
colSums(prop.table(b4, margin = 2)*100) # prints sum of column, should give 100 for each 
b4[is.nan(b4)] = 0 # View(b4)

### Filter again MAGs with not much data
# j = 1 # at least 1% across samples
# b5 = b4[apply(b4[,1:i], MARGIN = 1, function(x) any(x > j)), ] # applies the threshold
# # Save to file 
# now = format(Sys.time(), "%d-%m-%Y-%H-%M")
# wd = getwd()
# csvFileName = paste0(wd,"/results/metadmg/FINAL_MAGs_mapping/",now,"_DAMAGED_TAXAxSAMPLE_sig-",LR,"_Dam-",DamMin,".csv")
# fwrite(data.frame(b5), csvFileName,sep = "\t", row.names = TRUE,dec=",")
# pheatmap((b5[,-ncol(b5)]), scale = "row")
# pheatmap(mag.df.tax, scale = "row")

# 5. WEIGHT THE DAMAGE BY THE ABUNDANCE
# Damage will be weighted by multiplying i X abundance.
# We can scale all values in abundance table.
dmg.weigh = scale(mag.df[,1:33] * b4[,1:33])
### dmg.weigh.scale = t(apply(dmg.weigh,FUN = scale, MARGIN = 1)) ## Not need to re-scale BY ROW
fwrite(as.data.frame(dmg.weigh), file = "~/Projects/caribou/anvio_dirs/metadmg/samples_damage_MAGs_ABUNDANCE_WEIGHTED.txt", sep="\t",row.names = T, col.names = T, dec=",")
ComplexHeatmap::Heatmap(as.matrix(scale(dmg.weigh)),col = col, border = T)

## 5.1 INVERSE WEIGHTING DAMAGE ####

#' This helps representing what's truly active. 
#' We can weight abundances by the inverse of damage %, so that those with less damage (more active) 
#' get higher weights, which gets scaled after multiplying X abundance.

dmg.weigh.inv = scale(1/mag.df[,1:33] * b4[,1:33])
# fwrite(as.data.frame(dmg.weigh.inv), file = "~/Projects/caribou/anvio_dirs/metadmg/samples_damage_MAGs_ABUNDANCE_WEIGHTED_INVERSE.txt", sep="\t",row.names = T, col.names = T, dec=",")
ComplexHeatmap::Heatmap(as.matrix(scale(dmg.weigh.inv)),col = col, border = T)

### Just checking the scaled/non-scaled distributions...
# noscale = as.numeric(unlist(dmg.weigh))
# withscale = as.numeric(unlist(dmg.weigh))
# hist(noscale,breaks=30)
# plot(withscale,noscale)

stop("END PART 2!")

# EXTRA-PLOTS ####
# removing Nreads column
z = ncol(b5)
y = z - 1
# makes b5 to long table for ggplot heatmap below
b6 = melt(b5[,1:y])
# NMDS 
set.seed(42)
# b5t = t(b5)
nmds = vegan::metaMDS(b5[,1:y],trymax = 200)

# converting NMDS to long format for ggplot
data.scores = as.data.frame(vegan::scores(nmds,"site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site = rownames(data.scores)  # create a column of site names, from the rownames of data.scores
# data.scores$grp = colnames(nmds)  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores = as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species = rownames(species.scores)  # create a column of species, from the rownames of species.scores
vec.smp = as.numeric(unlist(lapply(strsplit(species.scores$species,split="_"),`[`,2)))
year = ifelse(vec.smp < 118,"2021","1978")
species.scores$year = year

# 1. NMDS of samples and taxa driving ordination ####
ggplot() + 
  geom_text_repel(data = data.scores, aes(x=NMDS1, y=NMDS2, label=site), size=3) +  # add the site labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, fill="grey"),size=0.5,alpha=0.5) +  # add the site labels
  geom_point(data = species.scores, aes(x = NMDS1, y=NMDS2, colour = year),size=3) + # add the point markers
  geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species, colour=year)) +  # add the species labels
  ggtitle("NMDS ordination of filtered samples (all plots following have filtered data)") +
  labs(col = "Samples", fill="Driving taxa") +
  theme_classic()

# filter only taxa and genus used and plot a damage plot prior to filtering 
# taxa1 = rownames(dd4)
# dm1 = df1[df1$tax_id %in% taxa1,]

## Filter from the main dataframe after main filters only the taxa we used (47 taxa; 27 samples)
taxa2 = unique(b6$Var1)
dm2 = a[a$tax_id %in% taxa2,]

# 2. TOTAL READS per SAMPLE ####
df77 = data.frame(Samples = colnames(dw1),Total_reads = colSums(dw1))
ggplot(data = df77, aes(x=Samples, y=Total_reads)) +
  geom_bar(stat="identity")  +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# 3. TAXA vs. MEAN READ LENGTH ####
ggplot(dm2, aes(x = mean_L, y=tax_id)) + geom_boxplot() + theme_test() +
  ylab("Taxa") + xlab("Mean Length (Bp)") + ggtitle("Taxa vs. mean read length")

# 4. TAXA vs. MAP_damage ####
# WE SEE SOME DAMAGE PATTERNS; SOME SAMPLES FROM 2021 HAVE MORE DAMAGE THAN AVERAGE.
ggplot(dm2, aes(x = MAP_damage, y = reorder(tax_id, MAP_damage, FUN=median, .desc=FALSE))) + geom_boxplot() + ggtitle("Taxa vs. damage") +
  ylab("Taxa") + xlab("MAP_damage")  + xlim(0,0.03) + theme_test()

# 5. HEATMAP PLOT WITH THE PROPORTION MAPPED TO EACH OF THE DEGRADED MAGS ####
DamMin = 0
p2 = ggplot(b6, aes(x = Var2, y = Var1, fill = log10(value))) +
  geom_tile(colour="lightgrey") + 
  theme_minimal() + 
  scale_fill_gradient(low="white", high="darkred") + 
  scale_y_discrete(limits=rev) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust =1)) + ggtitle(paste0("Damaged reads mapping to taxa/sample, minimum ",DamMin*100,"% damaged")) +
  ylab("Taxa") + 
  xlab("Samples") + labs(fill = "Logbase10 transformed \npercentage")
p2

stop("ENDS HERE!")