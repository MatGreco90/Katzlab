#This function converts any fasta file into table format (csv)
#If your input is an amino acid fasta file you can comment out line 10 and use line 11 instead

library(DECIPHER)
library(tidyverse)
library(magrittr)


fasta2csv<-function(x){
  dna <- readDNAStringSet(x)
  #dna <-readAAStringSet(x)
  dnadf<-as.data.frame(dna)
  dnadf$seq_name<-rownames(dnadf)
  final_df<-dnadf %>% 
    select(2,1) %>%
    set_colnames(c('seq_id','seq'))
    write.csv(final_df, file='my_final_file.csv',row.names = FALSE)
}
