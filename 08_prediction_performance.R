# load libraries
library(data.table)
library(tidyverse)
library(patchwork)
library(reshape2)
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

# get indv metadata
metadata <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_metadata.txt') %>% select(Sample, Admixture)
metadata$Sample <- str_remove(metadata$Sample, "_.*")
metadata <- metadata %>% unique() %>% mutate(Sample=case_when(
  Sample=='EU22'~'AF22', Sample=='EU36'~'AF36',
  Sample=='EU38'~'AF38', TRUE~Sample))
grm <- fread('/project/lbarreiro/USERS/daniel/katieData/WGS/35Samples.filtered.GRM.rel')
grm_ids <- fread('/project/lbarreiro/USERS/daniel/katieData/WGS/35Samples.filtered.GRM.rel.id', header=F)
grm_ids$V1 <- str_extract(grm_ids$V1, "[A-Z]{2}[0-9]+")
grm_ids$V2 <- str_extract(grm_ids$V2, "[A-Z]{2}[0-9]+")
rownames(grm) <- grm_ids$V1
colnames(grm) <- grm_ids$V2
grm <- grm %>% rownames_to_column('ID') %>% melt() %>% unique() %>% filter(variable=='AF04')
metadata <- inner_join(metadata, grm, by=c('Sample'='ID'))

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
    t <- data.frame(c(cor(full$avg_sum, full$adj, method='spearman'), 'adj', co, i),
                    c(cor(full$avg_sum, full$count, method='spearman'), 'count', co, i)) %>%
      t() %>% as.data.frame()
    
    # append to final data frame
    if (exists('final.corr.df_nobias')){
      final.corr.df_nobias <- rbind(final.corr.df_nobias, t)
    } else { final.corr.df_nobias <- t}
    
  }
}
final.corr.df_nobias$V1 <- as.numeric(final.corr.df_nobias$V1)
final.corr.df_nobias <- inner_join(final.corr.df_nobias, metadata, by=c('V4'='Sample'))

# plot data
nobias <- ggplot(final.corr.df_nobias) + geom_violin(aes(x=V2,y=V1)) +
  geom_boxplot(aes(x=V2,y=V1), width=0.1) + facet_wrap(~V3) +
  ylab('Spearman correlation coefficient') + xlab('ATAC data modality') +
  theme_bw()
nobias

# save plot figure
ggsave('prediction_results_nobias.pdf', height=5, width=5)

# see if prediction correlates with relatedness
# calculate Pearson correlation and p-value for each group
cor_data <- final.corr.df_nobias %>% filter(V4!='AF04') %>%
  group_by(V2, V3) %>% summarize(r_value=round(cor(V1, value), 2),
    p_value=format.pval(cor.test(V1, value)$p.value, digits = 3),
    .groups='drop') %>% mutate(label = paste0('r = ', r_value, '\nP = ', p_value))

# plot
grm_scatter <- final.corr.df_nobias %>% filter(V4!='AF04') %>% ggplot(.) + 
  geom_point(aes(y=V1,x=value)) + facet_grid(vars(V2), vars(V3)) + 
  stat_smooth(aes(x=value, y=V1), method='lm', geom='smooth', formula=y~x) + 
  geom_text(data=cor_data, aes(x=0.025, y=0.57, label=label), inherit.aes=FALSE, size=5, hjust=1) + 
  xlab('Genetic relatedness to training individual') + ylab('Spearman correlation') + 
  theme_bw()

# calculate Pearson correlation and p-value for each group
cor_data <- final.corr.df_nobias %>% filter(V4!='AF04') %>%
  group_by(V2, V3) %>% summarize(r_value=round(cor(V1, Admixture), 2),
    p_value=format.pval(cor.test(V1, Admixture)$p.value, digits = 3),
    .groups='drop') %>% mutate(label=paste0('r = ', r_value, '\nP = ', p_value))

# plot
adm_scatter <- final.corr.df_nobias %>% filter(V4!='AF04') %>% ggplot(.) + 
  geom_point(aes(y=V1,x=Admixture)) + facet_grid(vars(V2), vars(V3)) + 
  stat_smooth(aes(x=Admixture, y=V1), method='lm', geom='smooth', formula=y~x) + 
  geom_text(data=cor_data, aes(x=0.87, y=0.57, label=label), inherit.aes=FALSE, size=5, hjust=1) + 
  xlab('African admixture levels') + ylab('Spearman correlation') + 
  theme_bw()

grm_scatter + adm_scatter
ggsave('performance.vs.admixture_or_gr.scatter.pdf', height=10, width=13)

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