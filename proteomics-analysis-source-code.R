
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readr)
library(lemon)
library(scales)
library(ggpubr)
library(BioVenn)
library(ggrepel)
library(ggpmisc)
library(GGally)
library(RColorBrewer)

compareGenotypesOneTailed = function(df, wildtype1, wildtype2, wildtype3, mutant1, mutant2, mutant3){
  
  my.table = df %>% select(gene, 
                           paste0(wildtype1), 
                           paste0(wildtype2),
                           paste0(wildtype3),
                           paste0(mutant1),
                           paste0(mutant2),
                           paste0(mutant3))
  
  my.table['wildtype'] = ( my.table[paste0(wildtype1)] + my.table[paste0(wildtype2)] + my.table[paste0(wildtype3)] ) / 3
  my.table['mutant'] = ( my.table[paste0(mutant1)] + my.table[paste0(mutant2)] + my.table[paste0(mutant3)] ) / 3
  my.table['fc'] = (my.table['mutant'] + 0.01) / (my.table['wildtype'] + 0.01)
  my.table['lfc'] = log2(my.table['fc'])
  
  my.table = my.table %>% filter(!is.na(gene))
  
  pvals = c()
  for ( i in 1:length(my.table$gene) ){
    
    x1 = my.table[[2]][i]
    x2 = my.table[[3]][i]
    x3 = my.table[[4]][i]
    
    y1 = my.table[[5]][i]
    y2 = my.table[[6]][i]
    y3 = my.table[[7]][i]
    
    pvals[i] = t.test( c(x1,x2,x3), c(y1,y2,y3), alternative = "less", var.equal = T)$p.value
    
  }
  
  my.table['pval'] = pvals
  my.table.filtered = subset(my.table, lfc != -Inf)
  
  my.table.filtered = my.table.filtered %>% mutate(diff_expression = ifelse(pval < 0.05 & lfc >= 3, "up",
                                                                            ifelse(pval < 0.05 & lfc <= -3, "down", "none")))
  
  my.table.filtered = my.table.filtered %>% select(gene, wildtype1, wildtype2, wildtype3, mutant1, mutant2, mutant3, wildtype, mutant, fc, lfc, pval, diff_expression)
  
  return(my.table.filtered)
  
}

plotVolcano = function(df, bait){
  
  up = df %>% filter(diff_expression == "up")
  down = df %>% filter(diff_expression == "down")
  none = df %>% filter(diff_expression == "none")
  bait = df %>% filter(gene == bait)
  nUp = length(up$gene)
  
  annot = data.frame(x = c("right"), y = c("top"), lab = c(paste0(nUp)))
  
  p = ggplot() + 
    geom_point(data = none, aes(x = lfc, y = -log10(pval)), color = "grey60",alpha = 0.5, size = 1.2) +
    geom_point(data = up, aes(x = lfc, y = -log10(pval)), color = "green2", size = 1.5) +
    geom_point(data = bait, aes(x = lfc, y = -log10(pval)), color = "red", size = 2) +
    ylim(0,10) +
    xlim(-10,10) +
    theme_classic() + 
    theme(aspect.ratio = 1,
          axis.ticks.length=unit(.25, "cm"),
          text = element_text(size=20),
          axis.text = element_text(size = 20)) + 
    coord_cartesian() + 
    coord_capped_cart(bottom = "both", left = "both") + 
    geom_vline(xintercept = 3, lty = 2, color = "grey50") + 
    geom_vline(xintercept = -3, lty = 2, color = "grey50") +
    geom_hline(yintercept = 1.3, lty = 2, color = "grey50") +
    labs(y = bquote(-log[10](p-value)),
         x = bquote(log[2]("Fold Change")),
         fill = "") + 
    geom_text_npc(data = annot, aes(npcx = x, npcy = y, label = lab), size = 6) 
  return(p)
  
}


plotIupredBoxplot = function(df){
  
  names(df) = c("gene", "disorder", "sample")
  
  ymax = max(df$disorder) * 1.3
  
  p = ggplot(df, aes(x = sample, y = disorder)) + 
    geom_boxplot(width = 0.5) + 
    theme_classic() + 
    theme(aspect.ratio = 1,
          axis.ticks.length=unit(.25, "cm"),
          text = element_text(size=20),
          axis.text = element_text(size = 20)) + 
    coord_cartesian() + 
    coord_capped_cart(bottom = "both", left = "both") + 
    scale_y_continuous(limits = c(0, 0.8)) + 
    labs(y = "Protein Disorder (%)", x = "") 
  
  return(p)
  
  
  
  
  
}




