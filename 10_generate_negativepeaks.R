# load libraries
library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
# function to see if a given region overlaps with any region contained within a data frame
is_overlapping <- function(start, end, regions) {
  any((start >= regions$start & start <= regions$end) | 
        (end >= regions$start & end <= regions$end) | 
        (start <= regions$start & end >= regions$end))
}

# get the distribution of sizes of called peaks
peaks <- fread('/project/lbarreiro/USERS/daniel/katieData/ATACseq_peaks.bed') %>% select(V1,V2,V3) %>%
  mutate(len=V3-V2)
ggplot(peaks) + geom_histogram(aes(x=log10(len))) + theme_bw() +
  xlab('Peaks size log(10)')

# simulate peaks whose sizes follow the distribution of called peaks
non_peaks_length <- runif(nrow(peaks), min=max(300, mean(peaks$len)-sqrt(var(peaks$len))), 
                          max=mean(peaks$len)+sqrt(var(peaks$len))) %>% round() %>%
  as.data.frame() %>% mutate(chr=NA,start=NA,end=NA,V4='.',V5='.',V6='.',V7='.',V8='.',V9='.',summit=NA)
colnames(non_peaks_length)[1] <- 'len'

# also get background peaks (used in training)
for (f in seq(0,4)){
  t <- fread('/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/nonpeaks_fold'%&%f%&%'_negatives.bed')
  if (exists('backgroundpeaks')){
    backgroundpeaks <- rbind(backgroundpeaks, t)
  } else { backgroundpeaks <- t}
  rm(t)
}

# merge all peak files. those regions shouldn't be in the final output of this script
backgroundpeaks <- backgroundpeaks %>% select(V1,V2,V3) %>% 
  unique() %>% arrange(V1,V2,V3)
banned.peaks <- rbind(peaks[,1:3], backgroundpeaks) %>% arrange(V1,V2,V3)
rm(peaks, backgroundpeaks,f)

# get chr length and format data frame 
chr_sizes <- fread('/project/lbarreiro/USERS/daniel/katieData/WGS/chr_length.txt')
non_peaks_length$chr <- sample(unique(banned.peaks$V1), size=nrow(non_peaks_length), replace=T)
non_peaks_length <- left_join(non_peaks_length, chr_sizes, by=c('chr'='V1'))
rm(chr_sizes)
colnames(banned.peaks) <- c('chr','start','end')
colnames(non_peaks_length)[12] <- c('size')

# get regions for each chromosome at a time
# (find a way to run this in parallel?)
for (chrom in unique(non_peaks_length$chr)){
  print(chrom)
  
  # subset dfs for the specific chromosome
  banned_sub <- banned.peaks %>% filter(chr==chrom)
  nonpeaks_sub <- non_peaks_length %>% filter(chr==chrom)
  
  # for each future region
  for (i in 1:nrow(nonpeaks_sub)){
    # while the algorithm hasnt found a region that fits the criteria yet
    while (is.na(nonpeaks_sub$start[i])){
      # randomly select region start position and get its end position
      temptive_start <-  sample(1:(nonpeaks_sub$size[i]-nonpeaks_sub$len[i]), 1)
      temptive_end <- temptive_start+nonpeaks_sub$len[i]-1 
      
      # if the region doesn't overlap with any of the banned regions,
      # save the start and end positions; else, while loop runs again
      if(!is_overlapping(temptive_start, temptive_end, banned_sub)){
        nonpeaks_sub$start[i] <- temptive_start
        nonpeaks_sub$end[i] <- temptive_end
      }
    }
  }
  # replace regions in the empty data frame
  non_peaks_length[non_peaks_length$chr==chrom,] <- nonpeaks_sub
}
rm(banned_sub, nonpeaks_sub)

# remove unnecessary columns and prepare for chrombpnet
non_peaks_length <- non_peaks_length %>% select(-len,-size) %>% arrange(chr,start,end)  
non_peaks_length$summit <- round((non_peaks_length$end-non_peaks_length$start)/2)

# save output file
fwrite(non_peaks_length, '/project/lbarreiro/USERS/daniel/deeplearn_pred/DIY/negative_peaks.bed', col.names=F, sep='\t')