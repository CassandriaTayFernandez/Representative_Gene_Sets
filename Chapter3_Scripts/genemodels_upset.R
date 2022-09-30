#Input sequences are all of the gene names of sequences larger than 300bp and the BLAST results of every legume against other legumes (cds vs genome)
library(tidyverse)
library(data.table)
library(UpSetR)

names_species <- c('cajca', 'cicar', 'glyma', 'lotja', 'lupan', 'medtr', 'phavu', 'pissa', 'tripr', 'vigan', 'vigra', 'vigun')
names_species <- tibble(names_species)
colnames(names_species) <- c('species')
for (ref in names_species$species){
  ref_name <- fread(paste0(ref,'.cds.300.names.list'), header = F)  
  assign(paste0(ref,'_300bp_gene_list'),ref_name)
}

all_gene models <- tibble()
for (ref in names_species$species){
  all_refs <- tibble()
  for (query in names_species$species){
    if (query == ref) {next}
    ref_blast <- fread(paste0(query,'_query_',ref,'.cds.300_ref.tsv_hit_names_list.txt'), header = F)  
    ref_blast <- ref_blast %>% mutate(query_spec = query, ref_spec = ref)
    print(paste0(ref,'_ref_',query,'_cds_list'))
    all_refs <- bind_rows(all_refs,ref_blast)
  }
  all_gene models <- bind_rows(all_genes,all_refs)
  all_refs <- all_refs %>% distinct(V1, .keep_all=T)
  assign(paste0(ref,'_genes_shared'),all_refs)
} 
for (query in names_species$species){
  here <- all_genes[all_genes[, "query_spec"] == query,]
  assign(paste0(query,'_cds_genes'),here)
}
listInput <- list(cajca = cajca_cds_genes$V1,cicar = cicar_cds_genes$V1, glyma = glyma_cds_genes$V1, lotja = lotja_cds_genes$V1, lupan = lupan_cds_genes$V1, medtr = medtr_cds_genes$V1, phavu = phavu_cds_genes$V1, pissa = pissa_cds_genes$V1, tripr = tripr_cds_genes$V1, vigan = vigan_cds_genes$V1, vigra = vigra_cds_genes$V1, vigun = vigun_cds_genes$V1)

upset(fromList(listInput), nsets = 12, order.by = "freq")
