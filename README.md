# Katzlab

### CODE:

* Silva_16S_assing.R: this function takes a fasta file of rDNA sequences and uses the assignTaxonomy within the DADA2 package 
* Deploy_silva_assign.R: this script is needed to deploy the function Silva 16S assign on multiple fasta files allowing to control file names
* Tutorial_explore_silva_assign_data.r: this script takes the output of silva_assign.r and performs conversion to phyloseq object for further analyses
* fasta2csv.R: this function takes a fasta file as input and returns a table in csv format
* Deploy_fasta2csv.R:this script is needed to deploy the function fasta2csv on multiple fasta files allowing to control file names
* get_mono_clades.r: this function takes a newick tree as input and extracts the main monophyletic clades producing a summary table
* fetch_fastq.sh: this script automatically downloads fastq files from GenBank using SRA Run ID as input

## Author
Mattia Greco (mattia_greco@outlook.com)
