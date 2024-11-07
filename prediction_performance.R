library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('W:/deeplearn_pred/chrombpnet/test_prediction/outputs/')

list_of_files <- list.files(pattern='\\.txt$') %>%  
  gsub(pattern='_bwavg.txt', replacement='')

measured_adj <- fread('W:/katieData/ATACseq_peaks.filtered.paired.noDup.counts_batch.age.corrected.txt')
measured_count <- fread('W:/katieData/ATACseq_peaks.filtered.paired.noDup.counts')

for (id in list_of_files){
  print(id)
  prediction <- fread(id %&%'_bwavg.txt')
  measured_adj_red <- measured_adj %>% select(V1, all_of(id%&%'_NI'))
  measured_count_red <- measured_count %>% select(PeakID, all_of(id%&%'_NI'))

  full <- inner_join(prediction, measured_adj_red, by=c('V1'='V1')) %>% 
    inner_join(measured_count_red, by=c('V1'='PeakID')) %>% drop_na()
    
  colnames(full) <- c('name','size','covered','sum','mean0','mean','adj','count')
  
  t <- data.frame(c(cor(full$sum, full$adj, method='spearman'), 'adj'), 
  c(cor(full$sum, full$count, method='spearman'), 'count')) %>% t() %>%
    as.data.frame()
  
  if (exists('final.corr.df')){
    final.corr.df <- rbind(final.corr.df, t)
  } else { final.corr.df <- t}
  
}

final.corr.df$V1 <- as.numeric(final.corr.df$V1)
ggplot(final.corr.df) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1)


