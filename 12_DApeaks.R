# load libraries
library(data.table)
library(tidyverse)
library(DESeq2)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('W:/deeplearn_pred/DIY/prediction')

# get Katie's count matrix
katie_count <- fread('W:/katieData/ATACseq_peaks.filtered.paired.noDup.counts') %>%
  as.data.frame()
katie_count$PeakID <- gsub('_', ':', katie_count$PeakID)
row.names(katie_count) <- katie_count$PeakID
katie_count <- katie_count %>% select(-PeakID)

# make Katie's metadata
katie_meta <- colnames(katie_count) %>% as.data.frame()
colnames(katie_meta) <- 'ids'
katie_meta$condition <- substr(katie_meta$ids, 6, nchar(katie_meta$ids))
row.names(katie_meta) <- katie_meta$ids
katie_meta <- katie_meta %>% select(-ids)

# create Katie's DESeq2 dataset
katie_dds <- DESeqDataSetFromMatrix(countData=katie_count, colData=katie_meta, design=~condition)

# run differential analysis
katie_dds <- DESeq(katie_dds)

# get results
katie_results <- results(katie_dds, contrast=c('condition', 'NI', 'Flu'))
katie_sig_results <- katie_results[which(katie_results$padj<0.05),]

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
predicted_count <- predicted_count %>% select(-V1) %>% drop_na()

# make predicted metadata
predicted_meta <- colnames(predicted_count) %>% as.data.frame()
colnames(predicted_meta) <- 'ids'
predicted_meta$condition <- substr(predicted_meta$ids, 6, nchar(predicted_meta$ids))
row.names(predicted_meta) <- predicted_meta$ids
predicted_meta <- predicted_meta %>% select(-ids)

# create predicted DESeq2 dataset
predicted_dds <- DESeqDataSetFromMatrix(countData=predicted_count, colData=predicted_meta, design=~condition)

# run differential analysis
predicted_dds <- DESeq(predicted_dds)

# get results
predcited_results <- results(predicted_dds, contrast=c('condition', 'NI', 'Flu'))
predicted_sig_results <- predcited_results[which(predcited_results$padj<0.05),]

###









