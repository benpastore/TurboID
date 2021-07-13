setwd('C:/Users/Ian/OneDrive - The Ohio State University/work/eggd_data/brood')

brood <- read_csv('brood_size.csv')


brood

brood$Strain <- factor(brood$Strain, levels = c('N2','eggd-1','eggd-2','eggd-1/eggd-2'))
brood  

bc <- list(c('eggd-1','eggd-2'),c('N2','eggd-1'),c('N2','eggd-2'),c('N2','eggd-1/eggd-2'),c('eggd-1','eggd-1/eggd-2'))
bc
ggplot(data = brood, aes(x = Strain, y = TotalBrood, fill = Strain)) + 
  geom_boxplot(width = 0.7) +    
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20),
        axis.text = element_text(size = 20)) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") +
  geom_jitter() + 
  stat_compare_means(method = 't.test', comparisons = bc)
