library(tidytree)
library(tidyverse)
library(purrr)
library(phytools)

get_monophyletic_clades<- function (x) { #read tree
  
  my_tree<-phytools::read.newick(x)
  tree_df <- tidytree::as_tibble(my_tree) %>% 
    mutate(major_clade=gsub("_.*", "",label)) %>% #define major clades
    dplyr::mutate(isTip = !.data$node %in% .data$parent) #get tips
  
  my_off<- function (y) {
    out_res<- offspring(tree_df,y)
    return(out_res)
  }
  
  list_clades<-lapply(tree_df$parent,my_off)  #calculates offspring branches from each internal node and outputs a table
  
  tot_df <- map_df(list_clades, ~as.data.frame(.x), .id="id")  #merge tables for each combination node-offspring
  
  #print(tot_df)
  final_df_mono<-tot_df %>%
    na.omit() %>% 
    group_by(id) %>% #for each combination calculates the number of the major taxon present in each clade
    mutate(n_clades=length(unique(major_clade))) %>% 
    ungroup() %>% 
    filter(n_clades<3) %>% #remove poliphyletic clades (allowing for monophyletic w singletons)
    group_by(id) %>% 
    add_count(name='clade_size') %>%  #calculates size of the clade
    ungroup() %>% 
    group_by(id, major_clade) %>% #counts the occurrence of each major taxon in the clade
    add_tally(name='taxa_occurrence') %>% 
    ungroup() %>%
    mutate(freq=taxa_occurrence/clade_size) %>% 
    group_by(id) %>% #for allowing singletons i create treshold which is only relative to the size of each clade
    mutate(max_tres=(max(taxa_occurrence)-1)/clade_size, min_tres=1/clade_size) %>%
    ungroup() %>% 
    mutate(to_keep= ifelse(freq >=max_tres|freq==min_tres,'keep','discard')) %>% 
    group_by(id) %>% 
    mutate(final_decision= length(unique(to_keep))) %>%
    ungroup() %>% 
    filter(final_decision==1) %>%  #only keep the ones where there is not ambiguous nodes (all keep)
    filter(!clade_size==n_clades) %>%   #removes clades with size two and two taxa in it
    group_by(id) %>% 
    mutate(main_taxon= unique(major_clade[max(taxa_occurrence)]))%>% #define main id of the monophyletic clade
    ungroup() %>%  
    select(-c(11:15)) %>%  
    mutate(node_ran=paste(parent, node,sep='_')) %>% #removes duplicate clades only keeping the one with the bigger size based on parent-to-node combination
    group_by(node_ran) %>% 
    mutate(max_clade_comb=max(clade_size)) %>%
    ungroup() %>% 
    filter(!clade_size!=max_clade_comb) %>% 
    select(2:13) %>% 
    distinct() %>% 
    select(-c(1,2,3,6)) #write table
    write.csv(final_df_mono, file=gsub('.tre','.csv',x),row.names = FALSE) 
}
