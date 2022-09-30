#inputs are the BLAST alignments of the unique gene models (cds) against the NR database
#libraries ---- 
library(tidyverse)
library(data.table)
library(ggplot2)
names_species <- c('cajan', 'cicar', 'glyma', 'lotja', 'lupan','medtr', 'phavu', 'pissa', 'tripr', 'vigan', 'vigra', 'vigun')   
names_species[names_species == "cajan"] <- "cajca"
names_species <- tibble(names_species)
colnames(names_species) <- c('species')

all_unique_gene models <- tibble()
for (query in names_species$species){
  unique_gene models <- fread(paste0(query,'_300cds_no_other_legume_list.txt'), header = F)
  unique_gene models <- unique_gene models %>% mutate(query_spec = query)
#  unique_gene models <- unique_gene models %>% mutate(ref_spec=sapply(strsplit(V2,'[.]'),'[',1))  
  assign(paste0(query,'_300bp_unique_genes'),unique_genes)
  all_unique_gene models <- bind_rows(all_unique_genes,unique_genes)
}

totals <- all_unique_gene models %>% group_by(query_spec) %>% summarise(n=n())
 
#All legume Species vs NR BLAST p-value histogram
pvalue_matrix= tibble()
for (query in names_species$species){
  df <- fread(paste0(query,'_300cds__nolegume_butnrhits.txt'), header = F)  
  df <- df %>% 
    mutate(query_spec=sapply(strsplit(V1,'[.]'),'[',1))   %>% 
    arrange(V1,-V14,V13) %>% 
    distinct(V1,.keep_all = T)
  pvalue_matrix= bind_rows(pvalue_matrix, df)
}

translatetibble <- tibble(query_spec = unique(pvalue_matrix$query_spec),
                          query=c('Cajanus cajan',
                                  'Cicer arietinum',
                                  'Glycine max',
                                  'Lotus japonicus',
                                  'Lupinus angustifolius',
                                  'Medicago truncatula',
                                  'Pisum sativum',
                                  'Phaseolus vulgaris',
                                  'Trifolium pratense',
                                  'Vigna angularis',
                                  'Vigna radiata',
                                  'Vigna unguiculata')) %>% 
  mutate(query_spec=gsub('Medtr','medtr',query_spec))

all_unique_gene models <- all_unique_gene models %>% left_join(translatetibble) %>% 
  group_by(query_spec) %>% tally()

pvalue_matrix$pvalue_high <- as.factor((pvalue_matrix[,13]>0.01))
#nr hits anything----
nr_hits <- pvalue_matrix %>% 
  mutate(query_spec=gsub('Medtr','medtr',query_spec))%>% 
  group_by(query_spec, pvalue_high) %>% 
  summarise(n=n()) %>% 
  spread(pvalue_high,n) %>% 
  left_join(all_unique_genes) %>% 
  rename(bad=`TRUE`, good=`FALSE`) %>%
  mutate(bad=ifelse(is.na(bad),0,bad),
         good_rel=round(good/n*100,digits = 3),
         bad_rel=round(bad/n*100,digits = 3)) %>% 
  arrange(-good_rel) %>% 
  mutate(type='legume gene models that aligned with nr')

#hit_medians <- nr_hits %>% group_by(query_spec,good_rel) %>% summarise()
#colnames(hit_medians)<- c('query_spec','medians')

#plot
ggplot(nr_hits, aes(x='type',y=good_rel))+
  geom_boxplot(width = .4)+
  geom_jitter(aes(col = query_spec),width = 0.1, size = 4,alpha=0.5)

#nr query species does not hit itself ----

nr_noselfmatch= tibble()
for (query in names_species$species){
  df <- fread(paste0(query,'_300cds_nolegume_noselfnrhits.txt'), header = F) %>% 
    mutate(query_spec = query)
  nr_noselfmatch= bind_rows(nr_noselfmatch, df)
}
  
#nr query species hits itself ----
nr_list= tibble()
for (query in names_species$species){
  df <- fread(paste0(query,'_300cds__nolegume_butnrhits.txt'), header = F)  %>% 
    mutate(query_spec = query) %>% 
    distinct(V1,.keep_all = T)
  nr_list= bind_rows(nr_list, df)
}
nr_list_names <- nr_list %>% 
  select(V1)
nr_selfmatch <- nr_list_names %>% filter(!(V1 %in% nr_noselfmatch$V1)) 
nr_selfmatch_blasts <- nr_list %>% filter(V1 %in% nr_selfmatch$V1) 
nr_selfmatch_blasts$pvalue_high <- as.factor((nr_selfmatch_blasts[,13]>0.01))
nr_selfmatch_blasts <- nr_selfmatch_blasts %>%
  group_by(query_spec, pvalue_high) %>% 
  summarise(n=n()) %>% 
  spread(pvalue_high,n)

nr_selfmatch_blasts <- nr_selfmatch_blasts %>% 
  left_join(all_unique_genes) 

nr_selfmatch_blasts <- nr_selfmatch_blasts %>%
  rename(good=`FALSE`) %>% 
  mutate(bad=(n-good)) %>% 
  select(query_spec,good,bad,n) %>% 
  mutate(good_rel=round(good/n*100,digits = 3),
         bad_rel=round(bad/n*100,digits = 3)) %>% 
  arrange(-good_rel) %>% 
  mutate(type='nr hits aligning with the query')

# NR _ mismatch
#number of hits in NR which do not match their query species
hits <- nr_hits %>% 
  select(good,query_spec) %>% 
  rename(hits=good)


nr_noselfmatch_blasts <- nr_list_names %>% filter(V1 %in% nr_noselfmatch$V1)
nr_noselfmatch_blasts <- nr_list %>% filter(V1 %in% nr_noselfmatch_blasts$V1) 
nr_noselfmatch_blasts$pvalue_high <- as.factor((nr_noselfmatch_blasts[,13]>0.01)) 
taxids<- fread("D:/Rscripts/chap2/taxid")
colnames(taxids) <- c("query_spec","taxid")
nr_noselfmatch_blasts <- left_join(nr_noselfmatch_blasts,taxids, by=c("query_spec"))
nr_noselfmatch_blasts <- nr_noselfmatch_blasts %>%  filter(str_detect(as.character(taxid),as.character(V3), negate = T))

nr_noselfmatch_blasts <- nr_noselfmatch_blasts %>% 
  group_by(query_spec, pvalue_high) %>% 
  summarise(n=n()) %>% 
  spread(pvalue_high,n) 

nr_noselfmatch_blasts[is.na(nr_noselfmatch_blasts)] <- 0
nr_noselfmatch_blasts <- nr_noselfmatch_blasts %>% 
  left_join(all_unique_genes) %>% 
  rename(bad=`TRUE`, good=`FALSE`) %>% 
  mutate(good_rel=round(good/n*100,digits = 3),
         bad_rel=round(bad/n*100,digits = 3)) %>% 
  arrange(-good_rel) %>% 
  mutate(type='nr hits aligning with a different species to the query')

ggplot(nr_noselfmatch_blasts, aes(x='type',y=good_rel))+
  geom_boxplot(width = .4)+
  geom_jitter(aes(col = query_spec),width = 0.1, size = 4)+
  labs(x='', y = '% hit NR (e < 0.01)', title = 'Unique gene models hitting nr') +
  ylim(c(0, 100))
 
#combined plot
all_nr <- bind_rows(nr_hits, nr_selfmatch_blasts, nr_noselfmatch_blasts)
 
name_shenans <- all_nr
colnames(name_shenans)<- c('Species name','good','bad','n','good_rel','bad_rel','type')

ggplot(name_shenans, aes(x=type,y=good_rel))+
  geom_boxplot(width = .4)+
  geom_jitter(aes(col = `Species name` ),  width = 0.1, size = 3)+
  labs(x='', y = '% hit NR (e < 0.01)', title = 'Unique legume gene models aligning with the NR database')
