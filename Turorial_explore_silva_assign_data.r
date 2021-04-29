library(tidyverse)
library(janitor)
library(readxl)
library(phyloseq)

################################################################################
#Silva assignment results

#1)Kingdom
#2)Phylum
#3)Class
#4)Order
#5)Family
#6)Genus

#change to working directory where table were produced for the Silva_16S_assing.R function 
setwd(dir = "path_to_folder/")

#these lines will merge all the different tables and parse the sequence id to retrieve info such as clade, length and coverage
Silva_16S_tot_df<-rio::import_list(dir(pattern = ".csv"), rbind = TRUE, rbind_label = "source") %>%  
  select(2:10) %>% 
  row_to_names(row_number = 1) %>%   
  rename('cell_id'=tail(colnames(.),n=1)) %>%  
  mutate(cell_id=gsub('_rRNAseqs.csv','',cell_id)) %>%   
  filter(!seq=='seq') %>%   
  separate(seq_name,c('clade','bin','specimen_id','n_contigs','len','cov'), sep='_',extra = 'merge') %>% 
  mutate(len=ifelse(is.na(len),specimen_id,len),
         cov=ifelse(is.na(cov),n_contigs,cov),
         n_contigs=ifelse(clade=='Contig',bin,n_contigs)) %>% 
  select(-c(specimen_id,clade,bin)) %>% 
  mutate(len=as.numeric(gsub('Len','',len)),
         cov=as.numeric(gsub('Cov','',cov)),
         n_contigs=as.numeric(n_contigs)) %>%  
  filter(! Silva_Family=='Mitochondria') %>% #here we get rid of the sequences assigned to mitochondria 
  filter(cov>median(cov)) #we only keep sequences with a coverage higher than the median (excluded the mitochondria which are more abundant in the dataset)

#now proceed w phyloseq

#we are going to use the coverage as read abundance and use the package phyloseq to create an object that can be easily manipulated like in the amplicon studies
#to this end we need to create an otu table, a taxa table and a sample table (more info here https://joey711.github.io/phyloseq/ )


#we will now sum coverage belogning to sequences assigned to the same Family


otu_matrix<-Silva_16S_tot_df %>% 
  select(11,4,5:9) %>% 
  mutate(tax_lab=paste(Silva_Kingdom,Silva_Phylum,Silva_Class,Silva_Order,Silva_Family, sep="_")) %>%
  mutate(cell_id=gsub(' ','',cell_id)) %>% 
  group_by(cell_id,tax_lab) %>% 
  summarise(abund=sum(cov)) %>%
  ungroup() %>% 
  pivot_wider(names_from='cell_id' ,values_from ='abund' ) %>%  
  mutate_all(funs(replace_na(.,0))) %>%  
  mutate(otu=paste('OTU',1:length(tax_lab),sep='_'))

#taxa table

final_tax<-otu_matrix %>% 
  separate(tax_lab,into=c('Silva_Kingdom','Silva_Phylum','Silva_Class','Silva_Order','Silva_Family')) %>%   
  select(length(.),1:5) %>%  
  column_to_rownames('otu')


TAX_tab<-as.matrix(final_tax)

#otu table

otu_mat<-otu_matrix %>% 
  select(length(.),2:(length(.)-1))%>% 
  column_to_rownames('otu')


OTU_tab<-as.matrix(otu_mat)



#sample taxonomy file (needs to be in your working directory)

taxonomy_file<-read_xlsx('taxon_sheet_blast_annotated_taxonomy.xlsx') %>% 
  clean_names() %>% 
  select(3,4) %>% 
  rename('cell_id'='code')


sample_list<-Silva_16S_tot_df %>% 
  select(cell_id) %>% 
  separate(cell_id, into=c('cell_id','lkh_number'), sep='-') %>%  
  distinct() %>%  
  mutate(cell_id=gsub(' ','',cell_id)) %>%  
  left_join(taxonomy_file,by='cell_id') %>% 
  distinct() %>% 
  mutate(sample=paste(cell_id,lkh_number, sep='-')) %>% 
  select(sample,cell_id,lkh_number, taxon) %>% 
  column_to_rownames('sample')

OTU = otu_table(OTU_tab, taxa_are_rows = TRUE)
TAX = tax_table(TAX_tab)
samples = sample_data(sample_list)


#now that we have the three tables we can create a phuloseq object

phylo_dat <- phyloseq(OTU, TAX, samples)

#from here we can start analyze the data and play with different visualization techniques https://bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html#:~:text=The%20phyloseq%20package%20is%20a,taxonomic%20assignment%20of%20the%20OTUs.

