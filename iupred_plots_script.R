library(tidyverse)
setwd('C:/Users/Ian/OneDrive - The Ohio State University/work/eggd_data/iupred')
e1iupred<- read_tsv('eggd-1_iupred.txt')
e2iupred<- read_tsv('eggd-2_iupred.txt')


a<- ggplot(e1iupred, aes(x = POS, y = `IUPRED SCORE`)) +
  geom_line() +
  lims( y = c(0,1))

a
a + theme_classic() + geom_hline(yintercept = 0.5) + scale_x_continuous(breaks = seq(0, 556, 100))

b<-ggplot(e2iupred, aes(x = POS, y = `IUPRED SCORE`)) +
  geom_line()+
  ylim(c(0,1)) +
  xlim(c(0,666))
b + theme_classic() + geom_hline(yintercept = 0.5)+ scale_x_continuous(breaks = seq(0, 666, 100))
