# load libraries
library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction')

# set vectors of indv IDs, the conditions, and different model folds
ids <- c('AF04','AF06','AF10','AF12','AF14','AF16','AF18','AF20','AF22',
         'AF24','AF26','AF28','AF30','AF34','AF36','AF38','EU03','EU05',
         'EU07','EU09','EU13','EU15','EU17','EU19','EU21','EU25','EU27',
         'EU29','EU33','EU37','EU39','EU41','EU43','EU47')
conditions <- c('NI', 'Flu')
folds <- c(0,1,2,3,4)

# read data frames containing the measured ATACseq peak values
measured_adj <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_peaks.filtered.paired.noDup.counts_batch.age.corrected.txt')
measured_adj$V1 <- gsub('_', ':', measured_adj$V1)
measured_count <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_peaks.filtered.paired.noDup.counts')
measured_count$PeakID <- gsub('_', ':', measured_count$PeakID)

# for every condition and indv
for (co in conditions){
  for (i in ids){
    print(i)
    
    # get the respective measured values for all peaks
    measured_adj_red <- measured_adj %>% select(V1, all_of(i%&%'_'%&%co))
    measured_count_red <- measured_count %>% select(PeakID, all_of(i%&%'_'%&%co))
    
    # read the result of all folds and combine them
    for (f in folds){
      t <- fread(co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt') %>% select(V1,V4) %>%
          rename(peak=V1, !!paste0('fold',f,'_pred'):=V4)
      
      if (exists('prediction')){
          prediction <- inner_join(prediction, t, by=c('peak'))
      } else { prediction <- t}
    }
    prediction <- inner_join(prediction, measured_adj_red, by=c('peak'='V1')) %>%
      inner_join(measured_count_red, by=c('peak'='PeakID')) %>% drop_na() %>%
      rename(adj:=!!paste0(i,'_',co,'.x'), count:=!!paste0(i,'_',co,'.y')) %>% 
      mutate(cond=co)
    
    # get spearman correlation coefficient per fold
    t <- data.frame(c(cor(prediction$fold0_pred, prediction$adj, method='spearman'), 'adj', 0, co),
                    c(cor(prediction$fold0_pred, prediction$count, method='spearman'), 'count', 0, co),
                    c(cor(prediction$fold1_pred, prediction$adj, method='spearman'), 'adj', 1, co),
                    c(cor(prediction$fold1_pred, prediction$count, method='spearman'), 'count', 1, co),
                    c(cor(prediction$fold2_pred, prediction$adj, method='spearman'), 'adj', 2, co),
                    c(cor(prediction$fold2_pred, prediction$count, method='spearman'), 'count', 2, co),
                    c(cor(prediction$fold3_pred, prediction$adj, method='spearman'), 'adj', 3, co),
                    c(cor(prediction$fold3_pred, prediction$count, method='spearman'), 'count', 3, co),
                    c(cor(prediction$fold4_pred, prediction$adj, method='spearman'), 'adj', 4, co),
                    c(cor(prediction$fold4_pred, prediction$count, method='spearman'), 'count', 4, co)) %>%
      t() %>% as.data.frame() 
    rm(prediction)
    
    # append to final data frame
    if (exists('final.corr.df')){
      final.corr.df <- rbind(final.corr.df, t)
    } else { final.corr.df <- t}
  }
}
final.corr.df$V1 <- as.numeric(final.corr.df$V1)

# plot data
ggplot(final.corr.df) + geom_violin(aes(x=V2,y=V1,fill=V3)) +
  geom_boxplot(aes(x=V2,y=V1,fill=V3), width=0.4, position=position_dodge(width=0.9)) + facet_wrap(~V4) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality')  +
  theme_bw()

# save plot figure
ggsave('prediction_results_nonbias_perfold.pdf', height=5, width=8)