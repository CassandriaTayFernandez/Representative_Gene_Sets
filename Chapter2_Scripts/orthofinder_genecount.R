#input file is Orthogroups.GeneCount.tsv which is provided by R using default settings
library(tidyverse)
 
orthoall <- fread("Orthogroups.GeneCount.tsv", header=T)
colnames(orthoall) <- c('Orthogroup', 'Ca', 'Ci' ,'Gl' ,'Lo', 'Lu', 'Me' ,'Ph' ,'Pi', 'Tr', 'Va', 'Vr' ,'Vu','total')
.cols <- setdiff(colnames(orthoall), "Orthogroup")
orthoall$No.ofspecies <- rowSums(orthoall != "0") - 2

species_combination <- orthoall %>%
 gather(species,count,2:13) %>% filter(count !=0) %>% 
 group_by(Orthogroup) %>% summarize(present_species=paste(species,collapse = '_'))

orthoall2 <- orthoall %>% left_join(species_combination)
a <- orthoall2 %>% select(No.ofspecies,present_species,total)
colnames(a) <- c("V1","V2","V3")

ggplot(a, aes(x=V1,y=V3))+
 geom_boxplot(aes(group=V1))+
# geom_smooth(se=F)+
 scale_x_continuous(breaks=1:13, labels=1:13)+
 scale_y_continuous(limits = c(0,50))+
 labs(title = 'No. of genes in orthogroups on based on the number of legumes', x = 'No. of legume species', y='No. genes per orthogroup')
