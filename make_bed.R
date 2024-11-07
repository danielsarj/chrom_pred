library(tidyverse)
library(data.table)
setwd('/project/lbarreiro/USERS/daniel/katieData')

atac_peaks <- fread('ATACseq_peaks.filtered.paired.noDup.counts') %>%
    select('PeakID') %>% separate(PeakID, c('chr', 'start', 'end')) %>%
    filter(chr %in% c('chrX', 'chrY', 'chrM')==F)

atac_peaks$start <- as.integer(atac_peaks$start)
atac_peaks$end <- as.integer(atac_peaks$end)
atac_peaks <- atac_peaks %>% drop_na()

empty_cols <- c('4', '5', '6','7','8','9')
atac_peaks[ , empty_cols] <- '.'
atac_peaks <- atac_peaks %>% mutate(summit=round((end-start)/2))

atac_peaks

fwrite(atac_peaks, 'ATACseq_peaks.bed', sep='\t', col.names=F)
