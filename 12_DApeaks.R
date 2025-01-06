# load libraries
library(data.table)
library(tidyverse)
library(limma)
library(edgeR)
library(UpSetR)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction')

# get Katie's count matrix
katie_count <- fread('/project/lbarreiro/USERS/daniel/katieData/ATAC/corrected_expression.txt') %>%
  as.data.frame()
katie_count$V1 <- gsub('_', ':', katie_count$V1)
row.names(katie_count) <- katie_count$V1
katie_count <- katie_count %>% select(-V1)
katie_weights <- fread('/project/lbarreiro/USERS/daniel/katieData/ATAC/weights.txt')  %>%
  as.data.frame()
katie_weights$V1 <- gsub('_', ':', katie_weights$V1)
row.names(katie_weights) <- katie_weights$V1
katie_weights <- katie_weights %>% select(-V1)

# make Katie's metadata
katie_meta <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_metadata.txt')
katie_meta <- katie_meta[which(!katie_meta$Condition=='Mock'),]
katie_meta <- katie_meta[!grepl('EU122', katie_meta$Genotyping_ID),]
## factorize certain columns from meta data
katie_meta$Batch <- as.factor(katie_meta$Batch)
katie_meta$Genotyping_ID <- as.factor(katie_meta$Genotyping_ID)
katie_meta$Condition <- factor(katie_meta$Condition, levels=c('NI','Flu'))

# run differential analysis
katie_design <- model.matrix(~0+Genotyping_ID+Condition, data=katie_meta)
katie_design <- katie_design[,colSums(katie_design!=0)>0]
katie_fit <- lmFit(katie_count, design=katie_design, weights=katie_weights) %>% eBayes()

# get results
katie_results <- topTable(katie_fit, coef=ncol(katie_fit), number=Inf)

###

# get predicted count matrix
# set vectors of indv IDs, the conditions, and different model folds
ids <- c('AF04','AF06','AF08','AF10','AF12','AF14','AF16','AF18','AF20','AF22',
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
        prediction <- fread('chrombpnet_nobias_model/'%&%co%&%'_fold0_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
      } else {
        t <- fread('chrombpnet_nobias_model/'%&%co%&%'_fold'%&%f%&%'_'%&%i%&%'_bwavg.txt.gz') %>% select(V1,V4)
        prediction <- inner_join(prediction, t, by=c('V1'))
      }
    }
    
    # compute average predicted accessibility across folds
    averages <- rowMeans(prediction[,(ncol(prediction)-4):ncol(prediction)])
    prediction <- prediction %>% mutate(!!c(i %&%'_'%&% co):=as.integer(round(averages))) %>% select(V1, c(i %&%'_'%&% co))
    
    # append to final data frame
    if (exists('predicted_count')){
      predicted_count <- inner_join(predicted_count, prediction, by=c('V1'))
    } else { predicted_count <- prediction}
  }
}
predicted_count <- predicted_count %>% as.data.frame()
row.names(predicted_count) <- predicted_count$V1
predicted_count <- predicted_count %>% select(-V1) %>% drop_na() %>% DGEList() %>% 
  calcNormFactors() 

# predicted voom
predicted_groups <- c(rep('NI', length(ids)), rep('Flu', length(ids))) %>% 
  factor(levels=c('NI','Flu'))
predicted_ids <- rep(ids, length(conditions)) %>% as.factor()
predicted_design <- model.matrix(~0+predicted_ids+predicted_groups)
predicted_voom <- voom(predicted_count, predicted_design, plot=T)
predicted_fit <-  lmFit(predicted_voom) %>% eBayes()

# get results
predicted_results <- topTable(predicted_fit, coef=ncol(predicted_fit), number=Inf)

###

# compare katie's results to predicted results

# get shared peaks
katie_results_df <- katie_results %>% as.data.frame() %>% rownames_to_column('peakID')
predicted_results_df <- predicted_results %>% as.data.frame() %>% rownames_to_column('peakID')
combined_dfs <- inner_join(katie_results_df, predicted_results_df, by=c('peakID'))
common_peaks <- combined_dfs %>% select(peakID) %>% pull()

ggplot(combined_dfs) + geom_point(aes(x=logFC.x, y=logFC.y)) + 
  xlab('Measured DA logFC') + ylab('Predicted DA logFC') + geom_abline(slope=1, color='red') +
  stat_smooth(aes(x=logFC.x, y=logFC.y), method='lm', geom='smooth', formula=y~x) + theme_bw()
ggsave('predicted.measured.DAlogFCcorrelation.pdf', height=4, width=4)

# filter peaks by pval and combine dfs
katie_results_df <- katie_results_df %>% mutate(data='real', sig=adj.P.Val<0.05 & abs(logFC)>1) %>% 
  filter(peakID %in% common_peaks)
predicted_results_df <- predicted_results_df %>% mutate(data='pred', sig=adj.P.Val<0.05 & abs(logFC)>1) %>% 
  filter(peakID %in% common_peaks)
combined_DA <- rbind(katie_results_df, predicted_results_df)
fwrite(combined_DA, 'DApeaks_dataframe.txt', sep=' ')

# make volcano plots of DA peaks
ggplot(combined_DA) + geom_point(aes(x=logFC, y=-log10(adj.P.Val), color=sig)) +
  facet_wrap(~data) + theme_bw()
ggsave('DApeaks_volcanoplots.pdf', height=5, width=8)

# get list of peaks per data type (predicted and measured)
katie_sig_results <- katie_results_df %>% filter(adj.P.Val<0.05 & abs(logFC)>1) %>% 
  select(peakID) %>% pull()
predicted_sig_results <- predicted_results_df %>% filter(adj.P.Val<0.05 & abs(logFC)>1) %>% 
  select(peakID) %>% pull()
sig_list <- list(pred=predicted_sig_results, real=katie_sig_results)

# UpSet plot
pdf(file='DApeaks_upsetplot.pdf', height=5, width=6, onefile=F)
upset(fromList(sig_list), point.size=3.5, line.size=2, text.scale=c(1.3, 1.3, 1, 1, 2, 1.5))
dev.off()