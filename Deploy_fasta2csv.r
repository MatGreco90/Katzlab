#This script is needed to deploy the fasta2csv function on multiple files
#it requires the user to create an output folder where the csv files will be written


my_files<-list.files('input folder',pattern = '.fas',full.names = TRUE) # change to input folder

outdir<-'output_folder/'# change me to output folder

lapply(my_files,fasta2csv)
