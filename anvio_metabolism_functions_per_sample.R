# Read the matrices with the functions across genomes and translates them to functions per sample.
# Reads also as well the mean coverage normalized per sample.

library(data.table)
library(ggplot2)

# 1. READ FUNCTIONS AND MAGS ####
basedir = "~/HPC_data/caribou/anvio_dirs/runs/metabolism/functions_across_genomes/"
acc.enz = fread("~/HPC_data/caribou/anvio_dirs/runs/bone_degradome_collagenases/enzymes_list.txt")
COGs = unique(grep("COG", acc.enz$enzyme, value=TRUE))
annot = c("COG20_FUNCTION", "HMM_PF00082", "HMM_PF01752", "HMM_PF16499")
mags = fread("~/HPC_data/caribou/anvio_dirs/runs/FINAL_MAGs/FINAL_MAGs.list", header=F)
colnames(mags) = c("mag")

message("Remember note below!")
### CUSTOM PROFILES BUILT WITH CUSTOM COLLAGENASE DATABASES ARE EXCLUDED. 
# THESE -> "non_M09A_subgroup_AA_HMM", "M09A_subgroup_AA_HMM"
### RESULTS LOOK TOO DIVERGENT AND NO TRUE COLLAGENASES SEEM TO BE RETURNED.
### NEED TO LEARN HOW TO BUILD MY OWN HMM PROFILE FROM FASTAS.

func_mags = data.frame()
for (i in annot){
  mats0 = as.data.frame(fread(paste0(basedir,"/",i,"_matrix_functions_across_genomes-FREQUENCY.txt")))
  if(i == "COG20_FUNCTION"){
    # get only the columns for our cogs
    mats.f = mats0[mats0$COG20_FUNCTION_accession %in% COGs,][,c(-1,-(ncol(mats0)-1))]
    rownames(mats.f) = mats.f$COG20_FUNCTION_accession
    mats.f = mats.f[,-ncol(mats.f)]
  } else {
    mats.f = mats0[,c(-1,-(ncol(mats0)-1),-ncol(mats0))]
    rownames(mats.f) = i
  }
  miss = setdiff(mags$mag,colnames(mats.f))
  mats.f[ , miss] = 0      
  if (ncol(mats.f) != ncol(mats.f)) {
    stop("something wrong with matrix columns. More columns removed than it should!")
  }
  mats.fs = mats.f[ ,order(names(mats.f))]
  func_mags = rbind(func_mags,mats.fs)
}

# View(func_mags)

# 2. READ MEAN COVERAGE ####
# SAMPLE METADATA
options(scipen = 10)
sample_mtd = fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/SAMPLES_metadata.txt", dec = ",")

# Read coverage
# ALL BONE SAMPLES : 
mean_cov = as.data.frame(fread("~/HPC_data/caribou/anvio_dirs/runs/ALL_FINAL_MAGs_MERGED-SUMMARY/bins_across_samples/mean_coverage.txt"))
mags_cov = mean_cov$bins
mean_cov = mean_cov[,!colnames(mean_cov) %in% c("bins")]
rownames(mean_cov) = mags_cov
# head(mean_cov[1:10,1:10])

mean_cov.t = transpose(mean_cov)
colnames(mean_cov.t) = mags_cov
rownames(mean_cov.t) = colnames(mean_cov)
# head(mean_cov.t[1:10,1:10])

# Normalize it
mag_sums = rowSums(mean_cov.t)
mean_cov_NORM = (mean_cov.t/mag_sums)*100

cov_final = transpose(mean_cov_NORM)
rownames(cov_final) = colnames(mean_cov_NORM)
colnames(cov_final) = colnames(mean_cov)

# 3. CALCULATE FUNCTIONS PER SAMPLE ####
head(cov_final[1:10,1:10])
head(func_mags[1:10,1:10])
# fwrite(x = cov_final, file = "~/Projects/caribou/anvio_dirs/metabolism/FUNCTIONS_PER_SAMPLE/cov_final_normalized.txt", row.names = T, col.names = T, sep = "\t",quote = F)
# fwrite(x = func_mags, file = "~/Projects/caribou/anvio_dirs/metabolism/FUNCTIONS_PER_SAMPLE/functions_per_MAG.txt", row.names = T, col.names = T, sep = "\t",quote = F)

# Cross tables to get the coverage-corrected functions per sample
fx_samples = list()
# iterates sample names
for (s in colnames(cov_final)){
  values_per_sample = c()
# iterates mags
  for (f in rownames(func_mags)){
# sums the values of the function per mag in sample "s"
    value_s = sum(cov_final[,s] * func_mags[f,])
# puts them all together in a vector
    values_per_sample = c(values_per_sample, value_s)
  }
  fx_samples[[s]] = values_per_sample
}

# Make a dataframe functions x samples
df_fx_samples = data.frame(fx_samples)
rownames(df_fx_samples) = rownames(func_mags)

# RAW VALUES #
# fwrite(x = df_fx_samples, file = "~/Projects/caribou/anvio_dirs/metabolism/FUNCTIONS_PER_SAMPLE/FUNCTIONS_PER_SAMPLE-MATRIX_RAW.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# Check if 1978 samples have lower COVERAGE values than 2021.
message("No difference in raw coverage for both\n")
wilcox.test(colSums(df_fx_samples[,1:16]), colSums(df_fx_samples[,17:33]))

# Testing that raw values for HMM_PF01752 are higher for 2021
acc.enz # Check column number for HMM_PF01752 peptidase M9
# sample_mtd
# As should be expected, because its peptidase M9 ### ROW NUMBER 11
summary(unname(unlist(df_fx_samples[11,1:16])))
summary(unname(unlist(df_fx_samples[11,17:33])))

# 3.1 Normalization options #####
# Normalize it by SAMPLE: Summing each sample and making the proportion
sx = colSums(df_fx_samples)
fx_s = sweep(df_fx_samples, 2, sx, `/`)*100
# fwrite(x = fx_s, file = "~/Projects/caribou/anvio_dirs/metabolism/FUNCTIONS_PER_SAMPLE/FUNCTIONS_PER_SAMPLE-MATRIX.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# 3.1.1 Scale values per columns or rows ######
fx_s_SCALED = data.frame(t(apply(fx_s, 1, scale)))
rownames(fx_s_SCALED) = rownames(fx_s)
colnames(fx_s_SCALED) = colnames(fx_s)
fwrite(x = fx_s_SCALED, file = "~/Projects/caribou/anvio_dirs/metabolism/FUNCTIONS_PER_SAMPLE/FUNCTIONS_PER_SAMPLE-MATRIX_scaled_by_row.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# Statistical test for groups, by rows
raw_values_pvalues = sapply(1:nrow(fx_s_SCALED), function(i) wilcox.test(unname(unlist(fx_s_SCALED[i,1:16])), 
                                                      unname(unlist(fx_s_SCALED[i,17:33])))[c("p.value")])
pv_enzs0 = data.frame(rownames(fx_s_SCALED),as.numeric(raw_values_pvalues))
colnames(pv_enzs0) = c("ID","pv")
# Get the correct ID from the accessions
acc.enz$ID = c(paste0("HMM_",acc.enz[1:3,1]$enzyme),acc.enz[4:12,1]$enzyme)
pv_enzs = merge(pv_enzs0,acc.enz,by = "ID")[,c(5,1,2)]

# View heatmap
pheatmap::pheatmap(fx_s_SCALED, scale = "column")
sample_mtd

## END ##
