# load libraries and databases
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)
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
pred_annot <- annotatePeak(pred_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end) %>% 
  group_by(SYMBOL) %>% slice(which.min(abs(distanceToTSS)))
real_annot <- annotatePeak(real_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end) %>% 
  group_by(SYMBOL) %>% slice(which.min(abs(distanceToTSS)))

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
ggsave('GOenrichment.heatmap.pdf', height=20, width=10)
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

# GSEA
# make model-specific GRange objects for all peaks
pred_DA <- DAfile %>% filter(data=='pred') %>% dplyr::select(peakID) %>%
  separate(peakID, into=c('chr', 'start', 'end'), sep=':') %>% makeGRangesFromDataFrame()
real_DA <- DAfile %>% filter(data=='real') %>% dplyr::select(peakID) %>%
  separate(peakID, into=c('chr', 'start', 'end'), sep=':') %>% makeGRangesFromDataFrame()

# annotate peaks
pred_annot <- annotatePeak(pred_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end) %>% 
  group_by(SYMBOL) %>% slice(which.min(abs(distanceToTSS)))
real_annot <- annotatePeak(real_DA, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db') %>%
  as.data.frame() %>% mutate(peakID=seqnames%&%':'%&%start%&%':'%&%end) %>% 
  group_by(SYMBOL) %>% slice(which.min(abs(distanceToTSS)))

# perform gsea
pred_gsea <- DAfile %>% filter(data=='pred') %>% 
  right_join(pred_annot) %>% arrange(logFC) 
pred_ranks <- pred_gsea %>% select(geneId, logFC) 
pred_pathways <- reactomePathways(pred_ranks$geneId)
pred_ranks <- setNames(pred_ranks$logFC, pred_ranks$geneId)
pred_gsea <- fgsea(pathways=pred_pathways, stats=pred_ranks, minSize=15, maxSize=500)

pred.topPathwaysUp <- pred_gsea[ES > 0][head(order(pval), n=10), pathway]
pred.topPathwaysDown <- pred_gsea[ES < 0][head(order(pval), n=10), pathway]
pred.topPathways <- c(pred.topPathwaysUp, rev(pred.topPathwaysDown))
plotGseaTable(pred_pathways[pred.topPathways], pred_ranks, pred_gsea, gseaParam=0.5)

real_gsea <- DAfile %>% filter(data=='real') %>% 
  right_join(real_annot) %>% arrange(logFC) 
real_ranks <- real_gsea %>% select(geneId, logFC) 
real_pathways <- reactomePathways(real_ranks$geneId)
real_ranks <- setNames(real_ranks$logFC, real_ranks$geneId)
real_gsea <- fgsea(pathways=real_pathways, stats=real_ranks, minSize=15, maxSize=500)

real.topPathwaysUp <- real_gsea[ES > 0][head(order(pval), n=10), pathway]
real.topPathwaysDown <- real_gsea[ES < 0][head(order(pval), n=10), pathway]
real.topPathways <- c(real.topPathwaysUp, rev(real.topPathwaysDown))
plotGseaTable(real_pathways[real.topPathways], real_ranks, real_gsea, gseaParam=0.5)

joint_gsea <- rbind(mutate(pred_gsea, data='pred'), mutate(real_gsea, data='real'))
fwrite(joint_gsea, 'GSEA_joint.df.txt', sep=' ')