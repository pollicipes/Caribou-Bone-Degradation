# MetaDMG analysis; statistics about damage. Linear models, etc.

# Get rid of the RStudio error trying to show preview of data frames
options(rstudio.help.showDataPreview = FALSE)

library(tidyverse);
library(ggpubr);
library(dplyr);
library(ggplot2);
library(car);
library(pwr2);
library(data.table)
library(sjstats)
library(effectsize)

# Read samples group metadata
# Use dec="," (comma) to tell R that the decimals are represented by comma
mtd = data.frame(fread("/Users/jrodriguez/Projects/caribou/anvio_dirs/metabolism/09_SAMPLES_metadata.txt", dec = ","))

# Read mag damage per sample
in_dmg = fread(file = "~/Projects/caribou/anvio_dirs/metadmg/samples_damage_MAGs.txt", header = T, dec = ",")
mag.df = in_dmg %>% column_to_rownames(var = "V1") %>% t() %>% as.data.frame
# class(mag.df)

# ADD ID TO MERGE
mag.df$ID = rownames(mag.df)
mag.m =  merge(mag.df,mtd[,c(1:4)],by="ID")
mag.f = reshape2::melt(mag.m)

# Explicit factors to the df
mag.f$layer = factor(mag.f$layer, levels = unique(mtd$layer))
mag.f$year_date = factor(mag.f$year_date, levels = unique(mtd$year_date))
mag.f$profile = factor(mag.f$profile, levels = unique(mtd$profile))

# Test distribution fit:
library(fitdistrplus)
# Just selecting the vector with the numbers
year_1 = mag.f[mag.f$year_date == "year_1978",]$value
year_2 = mag.f[mag.f$year_date == "year_2021",]$value
descdist(year_1, discrete = FALSE) 

# LOOKS LIKE IT IS NOT NORMALLY DIST. Non parametric tests will be used.
hist(year_1,breaks=100);hist(year_2,breaks=100)
shapiro.test(year_1); shapiro.test(year_2)
summary(year_1); summary(year_2)

# 1. ANOVA TEST EFFECTS OF VARIABLES ####
# PLOTS differences:
# By years: #####
library(ggpubr)
mag.f %>% filter(year_date=="year_2021") %>% dim

plt.year = 
  ggboxplot(mag.f, x = "year_date", y = "value", color = "year_date", palette = c("darkgray","orange"), size = 1) + 
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("1978\nn=1056","2021\nn=1122")) +
  xlab("") + 
  ylab("Proportion deamination 1st mapped read base pair\n(Damage)") +
  ggtitle("Average damage per year" ,subtitle = "Each sample from group mapped to 66 MAGs") +
  guides(fill="none", color="none") + 
  stat_compare_means(label.x = 1.25, label.y = 0.07,size = 5)
plt.year
ggsave(filename = "~/Projects/caribou/FIGURES_PAPER//Damage_per_year_NEW.pdf",plot = plt.year)
pv_year = wilcox.test(year_1,year_2)
message("There are differences in damage % between years, wilcox pv: ", pv_year$p.value)

# Average values
wilcox.test((mag.f[mag.f$year_date == "year_1978",]$value), (mag.f[mag.f$year_date == "year_2021",]$value))[["p.value"]]

# By layer: #####
plt.layer = ggboxplot(mag.f, x = "layer", y = "value", color = "year_date",size = 1) + 
  theme(text = element_text(size = 20)) +
  xlab("") + 
  scale_color_manual(labels = c("1978","2021"), values = c("darkgray","orange")) +
  labs(color = NULL) +
  ylab("Proportion deamination 1st mapped read base pair\n(Damage)") +  # legend
  ggtitle("Damage per year/layer") + #,subtitle = "Interaction analyses by layer") +
  stat_compare_means(ref.group = ".all.", hide.ns = FALSE)
plt.layer
ggsave(filename = "~/Projects/caribou/FIGURES_PAPER/Damage_per_year_interaction_layer.pdf",plot = plt.layer)

lay.aov = aov(value ~ year_date * layer, data = mag.f)
summary(lay.aov)
effectsize::omega_squared(lay.aov, partial = FALSE)
# Interpret OMEGA squared:
# https://cran.r-project.org/web/packages/effectsize/vignettes/interpret.html
message("Looks like year_date has a medium effect size in explaining damage to DNA.")
message("Layer has a small effect.")

# # By excavation profile: #####
# plt.profile = ggboxplot(mag.f, x = "profile", y = "value", color = "year_date") +
#   theme(text = element_text(size = 20)) +
#   rotate_x_text(angle = 45) +
#   xlab("") + 
#   scale_color_manual(labels = c("Year 1978","Year 2021"), values = c("#E7B800","#00AFBB")) +
#   labs(color = NULL) +
#   ylab("Damage") + 
#   ggtitle("Damage per profile",subtitle = "Interaction analyses by profile")
# ggsave(filename = "~/Projects/caribou/REPORT_MAGS_Dec2023/PLOTS_R/Damage_per_year_interaction_profile.pdf",plot = plt.profile)
# 
# plt.profile
# pro.aov = aov(value ~ year_date + profile, data = mag.f) # Additive effect
# summary(pro.aov)
# effectsize::omega_squared(pro.aov, partial = FALSE)

# 2. MULTIPLE PAIRED TESTS FOR FINDING SIGNIFICANT VALUES FOR MICROBES. ####
# 2.1 DIFFERENCE BY YEAR: #####

species_list = as.character(unique(mag.f$variable))

results = data.frame()
# species = "year_2021_Bin_44_MAG_00002"
for (species in species_list) {
  subset_data = mag.f[mag.f$variable == species, ]
  xy_aov = aov(subset_data$value ~ subset_data$year_date, data = subset_data)
  omega2 = effectsize::omega_squared(xy_aov, partial = FALSE)$Omega2
  cohens_d = cohens_d(subset_data$value ~ subset_data$year_date)$Cohens_d
  wilcox_test_result = wilcox.test(subset_data$value ~ subset_data$year_date)
  # Store the results in a data frame
  results = rbind(results, data.frame(
    Species = species,
    p_value = wilcox_test_result$p.value,
    Omega2 = omega2,
    Cohens_d = cohens_d,
    mean_1978 = mean(subset_data[subset_data$year_date == "year_1978",]$value),
    mean_2021 = mean(subset_data[subset_data$year_date == "year_2021",]$value)
  ))
}

# PVALUE ADJUSTED Benjamini-Hochberg
fdrs = p.adjust(results$p_value, method = "BH")
results$padj = fdrs
results$logp = -log10(results$padj)
# Display the results
View(results)
library(data.table)
fwrite(results, file = "~/Projects/caribou/anvio_dirs/metadmg/DAMAGE_BY_SAMPLES_and_SPECIES.stats", sep="\t", col.names = T, row.names = FALSE, dec=",")

stop("ALL DONE!")
