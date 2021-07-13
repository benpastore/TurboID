
# analyze mass spec data for ian 

library(BioVenn)
library(ggrepel)
source("/Users/benjaminpastore/Desktop/Projects/proteomics-analysis-source-code.R")

setwd("/Users/benjaminpastore/Desktop/Projects/turboid/")

data = read_delim("turboID_spectral_counts_cdl_ife_glh_deps_mut16.txt", delim = "\t")
data = data %>% select(Protien, IFE_1, IFE_2, IFE_3, DEPS_1, DEPS_2, DEPS_3, N2_1, N2_2, N2_3)
data = data %>% rename(gene = Protien)

ife_dat = compareGenotypesOneTailed(data, "N2_1", "N2_2", "N2_3", "IFE_1", "IFE_2", "IFE_3")
ife_plot = plotVolcano(ife_dat, "ife-1")

deps_dat = compareGenotypesOneTailed(data, "N2_1", "N2_2", "N2_3", "DEPS_1", "DEPS_2", "DEPS_3")
deps_plot =plotVolcano(deps_dat, "deps-1")

volcano_plots = ggarrange(ife_plot, deps_plot, ncol = 2, nrow = 1)
pdf("DEPS-IFE-volcano-plots.pdf", height = 10, width = 10)
volcano_plots 
dev.off()

#-- overlap of ife up and deps up 
ifeUp = ife_dat %>% filter(diff_expression == "up") %>% select(gene)
depsUp = deps_dat %>% filter(diff_expression == "up") %>% select(gene)
draw.venn(list_x = ifeUp$gene, list_y = depsUp$gene, list_z = "")

pdf("DEPS-IFE-overlap.pdf",  height = 5, width = 5)
draw.venn(list_x = ifeUp$gene, list_y = depsUp$gene, list_z = "")
dev.off()
both = ifeUp %>% filter(gene %in% depsUp$gene) %>% select(gene)

write.table(ifeUp$gene, "IFE-turboID-up.txt", quote = F, sep = '\t', row.names = F, col.names = F)
write.table(depsUp$gene, "DEPS-turboID-up.txt", quote = F, sep = '\t', row.names = F, col.names = F)
write.table(both$gene, "Both-turboID-up.txt", quote = F, sep = '\t', row.names = F, col.names = F)

# must convert gene id to wbgene id or alternatively change parse protein fasta script to keep seqid

# Iupred score 
iupred_data = read_delim("average_disorder_proteins.txt", delim = "\t", col_names = c("gene","disorder"))
deps1_only = depsUp %>% filter(!gene %in% both$gene)

enriched_genes = read_delim("up-genes.txt", delim = "\t", col_names = c("gene","stableID"))
enriched_genes = enriched_genes %>% 
  left_join(iupred_data, by = c("stableID" = "gene")) %>% 
  mutate(enriched_in = ifelse(gene %in% both$gene, "IFE/DEPS",
                ifelse(gene %in% deps1_only$gene, "DEPS", "IFE")))

write.table(enriched_genes, "Table-S1-for-ian.txt", sep = "\t", quote = F, col.names = T, row.names = F)

both_iupred = enriched_genes %>% filter(enriched_in == "IFE/DEPS") %>% select(-gene)
either_iupred = enriched_genes %>% mutate(enriched_in = "either") %>% select(-gene)

N_both = length(both_iupred$stableID)
N_either = length(either_iupred$stableID)

set.seed(100)
random_both = iupred_data %>% sample_n(N_both, replace = F)
random_both = random_both %>% mutate(enriched_in = "random") %>% rename(stableID = gene)
both_iupred = rbind(both_iupred, random_both)

set.seed(5)
random_either = iupred_data %>% sample_n(N_either, replace = F)
random_either = random_either %>% mutate(enriched_in = "random") %>% rename(stableID = gene)
either_iupred = rbind(either_iupred, random_either)

both_plot_iupred = plotIupredBoxplot(both_iupred)
either_plot_iuped = plotIupredBoxplot(either_iupred)

wilcox.test(either_iupred$disorder, random_either$disorder)
wilcox.test(both_iupred$disorder, random_both$disorder)

pdf("iupred-boxplot.pdf",  height = 10, width = 10)
ggarrange(both_plot_iupred, either_plot_iuped, ncol = 2, nrow = 1)
dev.off()

####### GO analysis 
data = read_delim("78.csv", delim = ",")
colnames(data)

go_plot = function(gost_result, gosource, color){
  
  data = subset(gost_result, source == paste0(gosource))
  print(head(data))
  data2 = data %>%  add_column(x = runif(nrow(.), min = .40, max = .60))
  labels = data2 %>% top_n(3, wt = negative_log10_of_adjusted_p_value)
  
  p = ggplot() + 
    geom_point(data = data2, aes(x = x, y = negative_log10_of_adjusted_p_value, size = negative_log10_of_adjusted_p_value),color = color) + 
    theme_classic() + 
    theme(aspect.ratio = 1,
          axis.ticks.length=unit(.25, "cm"),
          text = element_text(size=20),
          axis.text = element_text(size = 20),
          axis.text.x = element_blank()) + 
    theme(legend.position = "none") +
    coord_cartesian() + 
    coord_capped_cart(bottom = "both", left = "both") + 
    geom_text_repel(data = labels,aes(x = x, y = negative_log10_of_adjusted_p_value, label = term_name), box.padding = 4) + 
    scale_x_continuous(limits = c(0,1), breaks = c(.25,.75))+
    scale_y_continuous(limits = c(0,20), breaks = seq(0,20,4)) +
    labs(x = paste0(gosource), y = "-log10(p value)")
  
  return(p)  
  
  
}

mf_plot = go_plot(data, "GO:MF", "darkviolet")
bp_plot = go_plot(data, "GO:BP", "springgreen3")
cc_plot = go_plot(data, "GO:CC", "indianred2")

go_analysis = ggarrange(mf_plot, bp_plot, cc_plot, ncol = 3, nrow = 1)
pdf("go-analysis.pdf", height = 20, width = 20)
go_analysis
dev.off()

