# load libraries
library(data.table)
library(tidyverse)
library(rjson)
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

# initialize lists of different classes of chromosomes 
test_chr <- list()
train_chr <- list()
valid_chr <- list()

# assess which chromosomes, per fold, belong to the test, train, or validation sets
for (f in folds){
  t <- fromJSON(file='/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/splits/fold_'%&%f%&%'.json')
  test_chr[[f+1]] <- t$test
  train_chr[[f+1]] <- t$train
  valid_chr[[f+1]] <- t$valid
  }

# for every condition and indv
for (co in conditions){
  for (i in ids){
    print(i)
    
    # get the respective measured values for all peaks
    measured_adj_red <- measured_adj %>% select(V1, all_of(i%&%'_'%&%co))
    measured_count_red <- measured_count %>% select(PeakID, all_of(i%&%'_'%&%co))
    
    # read the result of all folds and combine them
    for (f in folds){
      t <- fread('chrombpnet_nobias_model/'%&%co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt') %>% select(V1,V4) %>%
          rename(peak=V1, !!paste0('fold',f,'_pred'):=V4)
      
      if (exists('prediction')){
          prediction <- inner_join(prediction, t, by=c('peak'))
      } else { prediction <- t}
    }
    prediction <- inner_join(prediction, measured_adj_red, by=c('peak'='V1')) %>%
      inner_join(measured_count_red, by=c('peak'='PeakID')) %>% drop_na() %>%
      rename(adj:=!!paste0(i,'_',co,'.x'), count:=!!paste0(i,'_',co,'.y')) %>% 
      mutate(cond=co) 
    
    # for each fold, get the respective spearman coefficient for each set of chromosomes (test, train, or validation)
    for (f in folds){
      train_pred <- prediction %>% separate(peak, into=c('chr','start','end'), sep=':') %>%
        filter(chr %in% train_chr[[f+1]]) %>% select(all_of('fold'%&%f%&%'_pred'),adj,count,cond)
      test_pred <- prediction %>% separate(peak, into=c('chr','start','end'), sep=':') %>% 
        filter(chr %in% test_chr[[f+1]]) %>% select(all_of('fold'%&%f%&%'_pred'),adj,count,cond)
      valid_pred <- prediction %>% separate(peak, into=c('chr','start','end'), sep=':') %>% 
        filter(chr %in% valid_chr[[f+1]]) %>% select(all_of('fold'%&%f%&%'_pred'),adj,count,cond)
      
      t <- data.frame(c(cor(train_pred[1], train_pred$adj, method='spearman'), 'adj', f, co, 'train'),
                      c(cor(train_pred[1], train_pred$count, method='spearman'), 'count', f, co, 'train'),
                      c(cor(test_pred[1], test_pred$adj, method='spearman'), 'adj', f, co, 'test'),
                      c(cor(test_pred[1], test_pred$count, method='spearman'), 'count', f, co, 'test'),
                      c(cor(valid_pred[1], valid_pred$adj, method='spearman'), 'adj', f, co, 'valid'),
                      c(cor(valid_pred[1], valid_pred$count, method='spearman'), 'count', f, co, 'valid')) %>%
        t() %>% as.data.frame()
      
      # append to final data frame
      if (exists('final.corr.df')){
        final.corr.df <- rbind(final.corr.df, t)
      } else { final.corr.df <- t}
    }
    rm(prediction)
  }
}
final.corr.df$V1 <- as.numeric(final.corr.df$V1)

# plot data
ggplot(final.corr.df) + geom_violin(aes(x=V2,y=V1,fill=V5)) +
  geom_boxplot(aes(x=V2,y=V1,fill=V5), width=0.4, position=position_dodge(width=0.9)) + 
  facet_grid(cols=vars(V4),rows=vars(V3)) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality') +
  theme_bw()

# save plot figure
ggsave('performance_byfold_bychrsets.pdf', height=10, width=7)