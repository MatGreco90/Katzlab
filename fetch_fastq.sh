#Written by Mattia Greco

#requires sratoolkit on your machine https://www.ncbi.nlm.nih.gov/books/NBK242621/
#go to folder on your machine where you want to download the data

#!/bin/bash

#change path to statoolkit on your machine

export PATH=$PATH:/Users/katzlab/Desktop/sratoolkit.2.11.0-mac64/bin

#change to SRA Run (format= SRRnnnnnnn) identifier of interest, each between quotes and separated by space

mySRAs=("SRR8447774" "SRR8447785")

for i in "${mySRAs[@]}"
do
	export sraid=$i

	prefetch $sraid

	fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/"$sraid".sra

done

