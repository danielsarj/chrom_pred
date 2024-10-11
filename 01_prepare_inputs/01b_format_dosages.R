library(data.table)
library(tidyverse)
'%&%' = function(a,b) paste (a,b,sep='')
setwd('/project/lbarreiro/USERS/daniel/katieData/WGS/dosages')

# get indv. IDs
ids <- fread('35Samples.noIPSC.filtered.chr.traw.gz') %>% select(contains('WGS')) %>%
  colnames() %>% str_split_i('_1_', 1)

for (chr in c(1:22)){
  # open unformatted dosage file for chromossome 'chr'
  dosage <- fread('35Samples.noIPSC.filtered.chr'%&% chr %&%'.txt.gz')

  # rename columns
  colnames(dosage) <- c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'MAF', ids)

  # edit SNP_ID column
  dosage$SNP_ID <- dosage$CHR %&%'_'%&% dosage$POS %&%'_'%&% dosage$REF %&%'_'%&% dosage$ALT

  # save formatted dosage file
  fwrite(dosage, '35Samples.noIPSC.filtered.formatted.chr'%&% chr %&%'.txt', col.names=T, sep=' ')
  system('gzip 35Samples.noIPSC.filtered.formatted.chr'%&% chr %&%'.txt')
}