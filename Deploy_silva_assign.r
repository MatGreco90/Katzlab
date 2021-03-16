#Written by Mattia Greco (mattia_greco@outlook.com)

#this script is needed to deploy the function Silva 16S assign on multiple fasta files
#you will need to create a folder where you want your output 
#then write the path to the folder in line 7

outdir<-'Arcellinida_16S_table_Silva/' #change me to output folder

arcell_fasta<-list.files('Arcellinida_rRNA_bin',pattern = '.fasta',full.names = TRUE) # path to fasta file (input) folder


lapply(arcell_fasta,Silva_16S_assign)
