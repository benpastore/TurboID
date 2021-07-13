library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readr)
library(lemon)
library(scales)
library(ggpubr)
library(ggrepel)
library(ggpmisc)
library(GGally)
library(RColorBrewer)

setwd('C:/Users/Ian/OneDrive - The Ohio State University/work/eggd_data/mutant_intensity_quantificiation/intensity_quant')


int <- read_csv('intensity_rats_xl.csv')

int$genotype <- factor(int$genotype, c('wt','e1','e2','e1_2'))

comarisons <- list(c('wt','e1'),c('wt','e2'),c('wt','e1_2'),c('e1','e2'),c('e1','e1_2'),c('e2','e1_2'))


int

ggplot(data = int, aes(x = genotype, y = X9, fill = genotype)) + 
  geom_boxplot(width = 0.7) +    
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20),
        axis.text = element_text(size = 20)) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") +
  geom_jitter() + 
  stat_compare_means(method = 't.test', comparisons = comarisons)


