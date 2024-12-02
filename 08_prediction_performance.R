# load libraries
library(data.table)
library(tidyverse)
library(patchwork)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction')

# set vectors of indv IDs, the conditions, and different model folds
ids <- c('AF04','AF06','AF08','AF10','AF12','AF14','AF16','AF18','AF20','AF22',
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

# NO BIAS MODEL PREDICTION
# for every condition and indv
for (co in conditions){
  for (i in ids){
    print(i)
    
    # get the respective measured values for all peaks
    measured_adj_red <- measured_adj %>% select(V1, all_of(i%&%'_'%&%co))
    measured_count_red <- measured_count %>% select(PeakID, all_of(i%&%'_'%&%co))
    
    # read the result of all folds and combine them
    for (f in folds){
      if (f==0){
        prediction <- fread('chrombpnet_nobias_model/'%&%co%&%'_fold0_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
      } else {
        t <- fread('chrombpnet_nobias_model/'%&%co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
        prediction <- inner_join(prediction, t, by=c('V1'))
      }
    }
    
    # compute average predicted accessibility across folds
    averages <- rowMeans(prediction[,(ncol(prediction)-4):ncol(prediction)])
    prediction <- prediction %>% mutate(avg_sum=averages) %>% select(V1, avg_sum)
    full <- inner_join(prediction, measured_adj_red, by=c('V1'='V1')) %>%
      inner_join(measured_count_red, by=c('V1'='PeakID')) %>% drop_na()
    colnames(full) <- c('name','avg_sum','adj','count')
    
    # get spearman correlation coefficient
    t <- data.frame(c(cor(full$avg_sum, full$adj, method='spearman'), 'adj', co),
                    c(cor(full$avg_sum, full$count, method='spearman'), 'count', co)) %>%
      t() %>% as.data.frame()
    
    # append to final data frame
    if (exists('final.corr.df_nobias')){
      final.corr.df_nobias <- rbind(final.corr.df_nobias, t)
    } else { final.corr.df_nobias <- t}
    
  }
}
final.corr.df_nobias$V1 <- as.numeric(final.corr.df_nobias$V1)

# plot data
nobias <- ggplot(final.corr.df_nobias) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1) + facet_wrap(~V3) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality') +
  theme_bw()
nobias

# save plot figure
ggsave('prediction_results_nobias.pdf', height=5, width=5)

# BIAS MODEL PREDICTION
# for every condition and indv
for (co in conditions){
  for (i in ids){
    print(i)
    
    # get the respective measured values for all peaks
    measured_adj_red <- measured_adj %>% select(V1, all_of(i%&%'_'%&%co))
    measured_count_red <- measured_count %>% select(PeakID, all_of(i%&%'_'%&%co))
    
    # read the result of all folds and combine them
    for (f in folds){
      if (f==0){
        prediction <- fread('chrombpnetbias_model/'%&%co%&%'_fold0_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
      } else {
        t <- fread('chrombpnetbias_model/'%&%co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
        prediction <- inner_join(prediction, t, by=c('V1'))
      }
    }
    
    # compute average predicted accessibility across folds
    averages <- rowMeans(prediction[,(ncol(prediction)-4):ncol(prediction)])
    prediction <- prediction %>% mutate(avg_sum=averages) %>% select(V1, avg_sum)
    full <- inner_join(prediction, measured_adj_red, by=c('V1'='V1')) %>%
      inner_join(measured_count_red, by=c('V1'='PeakID')) %>% drop_na()
    colnames(full) <- c('name','avg_sum','adj','count')
    
    # get spearman correlation coefficient
    t <- data.frame(c(cor(full$avg_sum, full$adj, method='spearman'), 'adj', co),
                    c(cor(full$avg_sum, full$count, method='spearman'), 'count', co)) %>%
      t() %>% as.data.frame()
    
    # append to final data frame
    if (exists('final.corr.df_bias')){
      final.corr.df_bias <- rbind(final.corr.df_bias, t)
    } else { final.corr.df_bias <- t}
    
  }
}
final.corr.df_bias$V1 <- as.numeric(final.corr.df_bias$V1)

# plot data
withbias <- ggplot(final.corr.df_bias) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1) + facet_wrap(~V3) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality') +
  theme_bw()
withbias

# save plot figure
ggsave('prediction_results_withbias.pdf', height=5, width=5)

# merge both results
final.corr.df_bias <- final.corr.df_bias %>% mutate(model='bias')
final.corr.df_nobias <- final.corr.df_nobias %>% mutate(model='no.bias')
corr.merged <- rbind(final.corr.df_bias, final.corr.df_nobias)

# violin plots
ggplot(corr.merged) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1) + facet_grid(vars(V3),vars(model)) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality') +
  theme_bw()
ggsave('prediction_results_with.without.bias.pdf', height=5, width=5)

# scatter plot
corr.merged <- cbind(final.corr.df_bias, final.corr.df_nobias) 
colnames(corr.merged)[c(1,2,3,5)] <- c('bias_prediction', 'modality', 'condition','nobias_prediction')
corr.merged <- corr.merged %>% select(bias_prediction,modality,condition,nobias_prediction)
ggplot(corr.merged) + geom_point(aes(x=bias_prediction, y=nobias_prediction)) +
  facet_grid(vars(modality),vars(condition)) + 
  ylab('No-bias model Spearman correlation coefficient') + 
  xlab('Bias model Spearman correlation coefficient') +
  stat_smooth(aes(x=bias_prediction, y=nobias_prediction),
              method='lm', geom='smooth', formula=y~x) +
  geom_abline(slope=1) +
  theme_bw()
ggsave('scatterplot.prediction_results_with.without.bias.pdf', height=5, width=5)