#Alignment of the CDS of all legume species against the CDS of every other species, where d.tsv is the name of the species alignment: eg: cajca_cds_vs_other_legumes.tsv
library(tidyverse) 
library(data.table) 
library(UpSetR) 
names_species <- c('cajca', 'cicar', 'glyma', 'lotja', 'lupan', 'medtr', 'phavu', 'pissa', 'tripr', 'vigan', 'vigra', 'vigun')

all_matrix <- tibble()
colnames <- c('qseqid','sseqid','sscinames','scomnames','pident','length','mismatch','gapopen','qstart','qend','sstart','send','ppos','evalue','bitscore')
for (d in names_species){
  query_species = unlist(strsplit(d,'_'))[1]
  gene_matrix= tibble()
  for (f in list.files(d, '.tsv')){
    df <- fread(paste0(d,'/',f), header= FALSE) %>% 
      as_tibble() 
    colnames(df) = colnames
    df = df %>% 
      distinct(qseqid,sseqid,.keep_all = T) %>% 
      arrange(qseqid,bitscore,-evalue) %>% 
      distinct(qseqid,.keep_all = T) %>% 
      select(qseqid,evalue,bitscore)  %>% 
      mutate(species=sapply(strsplit(f,'_'),'[',5))
    }
        gene_matrix= bind_rows(gene_matrix, df)
  }
  gene_matrix = mutate(gene_matrix, query = query_species)
  all_matrix <- bind_rows(all_matrix, gene_matrix)
}

#stepwiser blast count per species (for plotting)
all_matrix_count <- all_matrix %>% 
  group_by(query,gene) %>% 
  summarize(count=n()) %>% 
  group_by(query,count) %>%
  summarize(n=n())

# gained genes per species added
all_extra_genes <- all_matrix %>%
  group_by(query,gene, species) %>% 
  select(query,gene,species) %>% 
  mutate(n=n()) %>% 
  spread(species,n)
  
#define start point
#data prep
species <- unique(all_extra_genes$query)
bootstrep=tibble()
for (i in 1:10){
  #species random order
  species <- species[sample(1:13,13)]
  #reorder by species order
  all_extra_genes <- all_extra_genes %>% 
    select(c('query','gene',species))
  #steps loop
  steps= tibble(x= 1, spec = species[1] , new_genes = filter(all_extra_genes, query==species[1]) %>% nrow())
  for(st in 2:13){
    x=st
    step = tibble(x= st , spec = species[st] ,new_genes = all_extra_genes %>% filter(query==species[st]) %>%  filter_at(seq(3,2+st,1), all_vars(is.na(.))) %>% nrow())
    steps = bind_rows(steps,step)
  }
  steps= steps %>% mutate(core=0)
  for (r in 2:nrow(steps)){
  steps[r,4] <- steps[r-1,3]+steps[r-1,4]
  }
  steps = steps %>% mutate(round=i)
  bootstrep =  bind_rows(bootstrep,steps)
}
 
# all steps bootstrapped

bootstrep <- bootstrep %>% mutate(all_genes= core+ new_genes) 
#points
ggplot(bootstrep, aes(x=x,y=all_genes))+
  geom_point()+
  scale_x_continuous(breaks=1:13, labels=1:13)+
  scale_y_continuous(limits = c(0,110000))+
  geom_smooth()+
  labs(title = 'Genes in legumes', x = 'steps', y='genes')
 
#boxplot

ggplot(bootstrep, aes(x=x,y=all_genes))+
  geom_boxplot(aes(group=x))+
  geom_smooth(se=F)+
  scale_x_continuous(breaks=1:13, labels=1:13)+
  scale_y_continuous(limits = c(0,110000))+
  labs(title = 'Genes in legumes', x = 'steps', y='genes')

# individual steps facetted

all_steps <- bootstrep %>% select(-all_genes) %>% 
gather(type,count,new_genes:core)
all_steps <- all_steps %>% mutate(type=fct_relevel(type,c('new_genes','core')),
                                  round= paste0('round ',round),
                                  spec=ifelse(type=='core','',spec),
                                  round=fct_relevel(round,paste('round',1:10)))


ggplot(all_steps, aes(x=x,y = count, fill= type))+
  geom_bar(stat='identity')+
  facet_wrap(~round, scales='free_x')+
  geom_text(aes(y = 2000,label=spec), angle =90, size=3, position=position_dodge(0.5), hjust=0)+
  scale_x_continuous(breaks=1:13)+
  labs(title = 'Genes in Legumes', x = 'Number of Species', y= 'Gene count')
