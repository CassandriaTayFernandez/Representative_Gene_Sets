#usage: take output file from orthogroup_counter_1.py and put into R
library(tidyverse)
legume_combos <- read_tsv('output.file')
legume_combos <- legume_combos %>% mutate(total_genes= new_total_genes - extra_genes)
legume_combos <- legume_combos %>% select(-new_total_genes)

new <- legume_combos %>% pivot_longer(-Number_of_species, names_to = 'type', values_to ='count')
new <- new %>% ggplot(aes(x=Number_of_species, y=count, fill=type)) + geom_bar(stat="identity")
new + labs( x = "Number of Legume Species", y = "Number of Orthogroups with genes") + labs(fill = "Type of Genes") + scale_x_continuous(breaks=1:12, labels=1:12)
