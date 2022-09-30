#input is Orthogroups.GeneCount.tsv produced from OrthoFinder
library(UpSetR)
library(tidyverse)
library(data.table)
Genecount <- fread('Orthogroups.GeneCount.tsv', header = T)
colnames(Genecount) <- names_species <- c('Orthogroup','cajan', 'cicar', 'glyma', 'lotja', 'lupan','medtr', 'phavu', 'pissa', 'tripr', 'vigan', 'vigra', 'vigun','Total') 

allgenes <- Genecount
Genecount2 <- Genecount
allgenes[allgenes == 0] <- NA
alllist <- na.omit(allgenes)

Genecount <- Genecount %>% filter(!(Genecount$Orthogroup %in% alllist$Orthogroup))

cajcaorthogenes <- Genecount %>% filter(!(cajan == 0))
cicarorthogenes <- Genecount %>% filter(!(cicar == 0))
glymaorthogenes <- Genecount %>% filter(!(glyma == 0))
lotjaorthogenes <- Genecount %>% filter(!(lotja == 0))
lupanorthogenes <- Genecount %>% filter(!(lupan == 0))
medtrorthogenes <- Genecount %>% filter(!(medtr == 0))
pissaorthogenes <- Genecount %>% filter(!(pissa == 0))
phavuorthogenes <- Genecount %>% filter(!(phavu == 0))
triprorthogenes <- Genecount %>% filter(!(tripr == 0))
viganorthogenes <- Genecount %>% filter(!(vigan == 0))
vigraorthogenes <- Genecount %>% filter(!(vigra == 0))
vigunorthogenes <- Genecount %>% filter(!(vigun == 0))

listInputX <- list(cajca = cajcaorthogenes$Orthogroup, cicar = cicarorthogenes$Orthogroup, glyma = glymaorthogenes$Orthogroup, lotja = lotjaorthogenes$Orthogroup, lupan = lupanorthogenes$Orthogroup, medtr = medtrorthogenes$Orthogroup, phavu = phavuorthogenes$Orthogroup, pissa = pissaorthogenes$Orthogroup, tripr = triprorthogenes$Orthogroup, vigan = viganorthogenes$Orthogroup, vigra = vigraorthogenes$Orthogroup, vigun = vigunorthogenes$Orthogroup) 
upset(fromList(listInputX), nsets = 12, order.by = "freq")
