# Written by Mattia Greco (mattia_greco@outlook.com)

#this function takes a fasta file of rDNA sequences and uses the assignTaxonomy within the DADA2 package (more here https://benjjneb.github.io/dada2/assign.html)
#the function outputs a table with sequence identifier, sequence and taxon assignment
#the function requires the user to download the Silva database silva_nr99_v138_train_set.fa.gz from here https://zenodo.org/record/3986799#.YFDC8ZNKjgE
#After downloading you will need to modify the path to the database in line 20
#to deploy function use the script 'deploy_silva_assign'

library(dada2)
library(DECIPHER)
library(tidyverse)


Silva_16S_assign<-function(x){
  dna <- readDNAStringSet(x)
  dnadf<-as.data.frame(dna) 
  dnadf$seq_name<-rownames(dnadf)
  colnames(dnadf)<-c('seq','seq_name')
  myseq=tolower(dnadf$seq)
  taxa <- assignTaxonomy(myseq, "Databases/silva_nr99_v138_train_set.fa.gz",
                         multithread=FALSE,tryRC = TRUE)
  colnames(taxa) <- paste("Silva", colnames(taxa), sep = "_")
  finaldf<-cbind(dnadf,taxa)
  file.name<-paste(outdir,gsub('Arcellinida_rRNA_bin/|.fasta','',x),'.csv',sep='') #change 'Arcellinida_rRNA_bin/' with name of the folder where your fasta files are located
  write.csv(finaldf, file=file.name)
}
