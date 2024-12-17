# load libraries and databases
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/prediction')
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
DAfile <- fread('DApeaks_dataframe.txt')

# make model-specific GRange objects
pred_DA <- DAfile %>% filter(data=='pred', sig==TRUE) %>% dplyr::select(peakID) %>%
  separate(peakID, into=c('chr', 'start', 'end'), sep=':') %>% makeGRangesFromDataFrame()
real_DA <- DAfile %>% filter(data=='real', sig==TRUE) %>% dplyr::select(peakID) %>%
  separate(peakID, into=c('chr', 'start', 'end'), sep=':') %>% makeGRangesFromDataFrame()

# annotate peaks
pred_annot = annotatePeak(pred_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end)
real_annot = annotatePeak(real_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end)

# get info about peak subsets (real/pred only; intersection)
pred_peaks <- DAfile %>% filter(data=='pred', sig==TRUE) %>% dplyr::select(peakID) %>% pull()
real_peaks <- DAfile %>% filter(data=='real', sig==TRUE) %>% dplyr::select(peakID) %>% pull()

pred.only_annot <- pred_annot %>% filter(peakID %in% pred_peaks & !(peakID %in% real_peaks))
real.only_annot <- real_annot %>% filter(peakID %in% real_peaks & !(peakID %in% pred_peaks))
pred.real_annot <- pred_annot %>% filter(peakID %in% pred_peaks & peakID %in% real_peaks)

# quick summary
pred.only_annot$annotation <- gsub("\\(.*?\\)", "", pred.only_annot$annotation)
pred.only_annot$annotation <- trimws(pred.only_annot$annotation)
pred.only_summary <- pred.only_annot %>% group_by(annotation) %>% 
  summarise(count=n(), perc=(n()/nrow(pred.only_annot))*100) %>%
  mutate(DApeak='pred')

real.only_annot$annotation <- gsub("\\(.*?\\)", "", real.only_annot$annotation)
real.only_annot$annotation <- trimws(real.only_annot$annotation)
real.only_summary <- real.only_annot %>% group_by(annotation) %>% 
  summarise(count=n(), perc=(n()/nrow(real.only_annot))*100) %>%
  mutate(DApeak='real')

pred.real_annot$annotation <- gsub("\\(.*?\\)", "", pred.real_annot$annotation)
pred.real_annot$annotation <- trimws(pred.real_annot$annotation)
pred.real_summary <- pred.real_annot %>% group_by(annotation) %>% 
  summarise(count=n(), perc=(n()/nrow(pred.real_annot))*100)  %>%
  mutate(DApeak='pred&real')

# barplot idk
merged_summary <- rbind(pred.only_summary, real.only_summary, pred.real_summary)
ggplot(merged_summary) + geom_col(aes(x=annotation, y=perc, fill=DApeak), position='dodge') +
  theme_bw() + coord_flip()
ggsave('peaks.annotation.barplot.pdf', height=6, width=5)

# enriched GO terms
pred.GO <- enrichGO(pred.only_annot$geneId, org.Hs.eg.db, ont='MF') %>% 
  as.data.frame() %>% mutate(DApeak='pred')
real.GO <- enrichGO(real.only_annot$geneId, org.Hs.eg.db, ont='MF') %>% 
  as.data.frame() %>% mutate(DApeak='real')
pred.real.GO <- enrichGO(pred.real_annot$geneId, org.Hs.eg.db, ont='MF') %>% 
  as.data.frame() %>% mutate(DApeak='pred&real')

# join dfs
jointGOs <- rbind(pred.GO, real.GO, pred.real.GO) %>% select(Description, FoldEnrichment, p.adjust, DApeak)
ggplot(jointGOs) + geom_point(aes(x=DApeak, y=Description, color=FoldEnrichment, size=p.adjust)) + 
  theme_bw() + scale_size(transform='reverse')
ggsave('GOenrichment.heatmap.pdf', height=20, width=7)
fwrite(jointGOs, 'GOenrichment_joint.df.txt', sep=' ')

# enriched pathways terms
pred.path <- enrichPathway(pred.only_annot$geneId) %>% as.data.frame() %>% mutate(DApeak='pred')
real.path <- enrichPathway(real.only_annot$geneId) %>% as.data.frame() %>% mutate(DApeak='real')
pred.real.path <- enrichPathway(pred.real_annot$geneId) %>% as.data.frame() %>% mutate(DApeak='pred&real')

# join dfs
jointpaths <- rbind(pred.path, real.path, pred.real.path) %>% select(Description, FoldEnrichment, p.adjust, DApeak)
ggplot(jointpaths) + geom_point(aes(x=DApeak, y=Description, color=FoldEnrichment, size=p.adjust)) + 
  theme_bw() + scale_size(transform='reverse')
ggsave('Pathwayenrichment.heatmap.pdf', height=20, width=10)
fwrite(jointpaths, 'Pathwayenrichment_joint.df.txt', sep=' ')