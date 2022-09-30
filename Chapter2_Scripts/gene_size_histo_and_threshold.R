#Input is the tblastn (output 6) of peas aligned with all other legumes. 
library(tidyverse)
library(data.table) 

blasts <- fread('renamed.pissa.tblastn.chapt2.tsv',header = T,nThread=8)
colnames(blasts) <- c('qseqid', 'sseqid','pident' ,'length', 'mismatch', 'gapopen' ,'qstart' ,'qend', 'sstart', 'send', 'ppos' ,'evalue' ,'bitscore')
blasts <- blasts %>% arrange(qseqid,-pident,evalue,-bitscore) %>% 
 distinct(qseqid, .keep_all = T) %>% 
 mutate(evalue=round(evalue,digits = 2)) %>% 
 rename(subject_length=length)

# Gene size ----
#download and read pea cds (needs package 'BIOSTRINGS' installed)
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("Biostrings")

peacds <- Biostrings::readDNAStringSet('pissa.formatted.renamed.fna') %>% as.data.frame() %>% rownames_to_column(var='gene') %>% 
 rename(seq=x) %>% as_tibble()

#pissa.list.txt is the list of all pea genes in the cds file.
pea_list <- fread('pissa.list.txt', header = F)
pea_list_seq <- peacds %>% filter(gene %in% pea_list$V1) %>% 
 #count sequence length for gene size
 mutate(length=nchar(seq)) %>% 
 #add blast data with left_join()
 left_join(blasts, by=c('gene'='qseqid')) %>% 
 #classify hits by e-value
 mutate(type=ifelse(is.na(evalue),'Did not align with legumes',
      ifelse(evalue <= 0.01,'Aligns with legumes evalue < 0.01','Aligns with legumes evalue > 0.01')),
   type=fct_relevel(type,c('Aligns with legumes evalue < 0.01','Aligns with legumes evalue > 0.01','Did not align with legumes')))

#length summary table
summary <- pea_list_seq %>% 
 group_by(type) %>% 
 summarize(average=mean(length),
   median=median(length),
   sd=sd(length),
   n=n(),
   se=sd/sqrt(n))

cbPalette <- c("#56B4E9", "#009E73", "#CC79A7")
#plot histogram
ggplot(pea_list_seq) + 
 geom_histogram(aes(length,fill=type), binwidth = 50, alpha=0.8) +
 facet_wrap(~type)+
 scale_x_continuous(limits = c(0,3000))+
 labs(title = 'Pea gene length compared to the quality of alignment with other legumes', x='Gene Size (bp)',y='Number of hits with the pea CDS')+
 theme(legend.position = c(0.84, 0.8)) + scale_fill_manual(values=cbPalette)

#threshold estimation----
#make 100bp bins for sequences between 0 and 3000 base pairs
bins <- seq(0,3000,100)
# categorise data into a bin
threshold <- pea_list_seq %>% 
 select(length,type) %>%
 #adjust everything bigger than 3000bp to be shown as 3000bp
 mutate(length=ifelse(length >= 3000, 3000,length),
   #bin data
   bin=cut(length,bins, right=F, include.lowest = T)) %>% 
 #summarize
 group_by(type,bin) %>% 
 summarize(n_bin=n())

#generate n per bin
bin_size <- threshold %>%
 group_by(bin) %>% 
 summarize(n=sum(n_bin))
#add bins to data and calculate relative amounts of hit types per bin 
threshold <- threshold %>% 
 #add bin data
 left_join(bin_size) %>% 
 #calculate relative %
 mutate(perc=n_bin/n*100,
   #adjust bin text
   bin_mod=sapply(strsplit(as.character(bin),','),'[',2),
   bin_mod=gsub(')','',bin_mod,fixed = T)) 
#Plot estimates

ggplot(threshold, aes(as.numeric(bin_mod),perc,col=type))+
 geom_smooth(se = F)+
 scale_x_continuous(breaks = bins)+
 theme(axis.text.x = element_text(size=6))+labs(x='Gene size bin (bp)', y = 'percent (%)', title='Threshold of the length of a true gene')+
 theme(legend.position = c(0.8, 0.5)) + scale_colour_manual(values=cbPalette)
