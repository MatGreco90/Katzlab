#uses a fasta file to subset as input 
#modify path at line 10 to txt file with sequence ID list

library(Biostrings)
library(tidyverse)
library(DECIPHER)

subset_fasta<- function(x) {
  
  list_seqs_file<-read_table('list_seq_to_grepl.txt',col_names = FALSE) %>% #change path to file with list of sequences ID to subset
    mutate(X1=ifelse(grepl('>',X1),gsub('>','',X1),X1))
  
  dna <- readDNAStringSet(x)
  dnadf<-as.data.frame(dna)
  dnadf$seq_name<-rownames(dnadf)
  colnames(dnadf)<-c('seq','seq_ID')
  
  final_df<-dnadf %>%
    filter(seq_ID %in% list_seqs_file$X1) 
    
    
    
  seq = final_df$seq
  names(seq) = final_df$seq_ID
  dna_final = DNAStringSet(seq)
  writeXStringSet(dna_final, "subsetted.fasta")
  
}



 

