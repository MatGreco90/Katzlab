#This function converts any fasta file into table format (csv)
#If your input is an amino acid fasta file you can comment out line 10 and use line 11 instead
#to deploy the functin use the deploy_fasta2csv.r script

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
   file.name<-paste(outdir,gsub('fasta_to_convert/|.fasta','',x),'.csv',sep='') #change fasta_to_convert to your input folder  
    write.csv(final_df, file='my_final_file.csv',row.names = FALSE)
}
