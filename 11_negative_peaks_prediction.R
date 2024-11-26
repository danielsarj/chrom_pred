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

# for every condition and indv
for (co in conditions){
  for (i in ids){
    print(i)
    
    # read the result of all folds and combine them
    for (f in folds){
      if (f==0){
        prediction <- fread('negative_preds/'%&%co%&%'_fold0_'%&%i%&%'_nobias_bwavg.txt.gz') %>% select(V1,V4)
      } else {
        t <- fread('negative_preds/'%&%co%&%'_fold'%&%f%&%'_'%&%i%&%'_nobias_bwavg.txt.gz') %>% select(V1,V4)
        prediction <- inner_join(prediction, t, by=c('V1'))
      }
    }
    
    # compute average predicted accessibility across folds
    averages <- rowMeans(prediction[,(ncol(prediction)-4):ncol(prediction)])
    prediction <- prediction %>% mutate(avg_sum=averages) %>% select(avg_sum) %>%
      mutate(model='neg_peaks',condition=co)

    # append to final data frame
    if (exists('final.df_nobias')){
      final.df_nobias <- rbind(final.df_nobias, prediction)
    } else { final.df_nobias <- prediction}
    
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
    prediction <- prediction %>% mutate(avg_sum=averages) %>% select(avg_sum) %>%
      mutate(model='called_peaks',condition=co)
    
    # append to final data frame
    if (exists('final.df_nobias')){
      final.df_nobias <- rbind(final.df_nobias, prediction)
    } else { final.df_nobias <- prediction}
  }
}

final.df_nobias <- final.df_nobias %>% filter(avg_sum!=1.149090e+287)

# plot data
ggplot(final.df_nobias) + geom_boxplot(aes(y=log10(avg_sum), x=model)) + 
  theme_bw() + facet_wrap(~condition)
ggsave('boxplot_peaksizes_called.vs.negative.pdf', height=5, width=5)