library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction')

ids <- c('AF04','AF06','AF10','AF12','AF14','AF16','AF18','AF20','AF22',
         'AF24','AF26','AF28','AF30','AF34','AF36','AF38','EU03','EU05',
         'EU07','EU09','EU13','EU15','EU17','EU19','EU21','EU25','EU27',
         'EU29','EU33','EU37','EU39','EU41','EU43','EU47')
conditions <- c('NI', 'Flu')
folds <- c(0,1,2,3,4)
measured_adj <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_peaks.filtered.paired.noDup.counts_batch.age.corrected.txt')
measured_adj$V1 <- gsub('_', ':', measured_adj$V1)
measured_count <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_peaks.filtered.paired.noDup.counts')
measured_count$PeakID <- gsub('_', ':', measured_count$PeakID)

for (co in conditions){
  for (i in ids){
    print(i)
    measured_adj_red <- measured_adj %>% select(V1, all_of(i%&%'_'%&%co))
    measured_count_red <- measured_count %>% select(PeakID, all_of(i%&%'_'%&%co))
    
    for (f in folds){
      if (f==0){
        prediction <- fread(co%&%'_fold0_'%&%i%&%'_bwavg.txt') %>% select(V1,V4)
      } else {
        t <- fread(co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt') %>% select(V1,V4)
        prediction <- inner_join(prediction, t, by=c('V1'))
      }
    }
    averages <- rowMeans(prediction[,(ncol(prediction)-4):ncol(prediction)])
    prediction <- prediction %>% mutate(avg_sum=averages) %>% select(V1, avg_sum)
    full <- inner_join(prediction, measured_adj_red, by=c('V1'='V1')) %>%
      inner_join(measured_count_red, by=c('V1'='PeakID')) %>% drop_na()
    colnames(full) <- c('name','avg_sum','adj','count')
    
    t <- data.frame(c(cor(full$avg_sum, full$adj, method='spearman'), 'adj', co),
                    c(cor(full$avg_sum, full$count, method='spearman'), 'count', co)) %>%
      t() %>% as.data.frame()
    
    if (exists('final.corr.df')){
      final.corr.df <- rbind(final.corr.df, t)
    } else { final.corr.df <- t}
    
  }
}

final.corr.df$V1 <- as.numeric(final.corr.df$V1)
ggplot(final.corr.df) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1) + facet_wrap(~V3) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality')

ggsave('prediction_results_nonbiasmacromodel.pdf', height=5, width=5)