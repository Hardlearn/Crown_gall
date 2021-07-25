library(dplyr)
library(ggplot2)
## 1. loading data
Am_hog <- read.table("./1_fliter/Am_hog_final.txt")
colnames(Am_hog) <- c("gene",'phylostrata',"hog")
Ej_hog <- read.table("./1_fliter/Ej_hog_final.txt")
colnames(Ej_hog) <- c("gene","phylostrata","hog")
Go_hog <- read.table("./1_fliter/Go_hog_final.txt")
colnames(Go_hog) <- c("gene","phylostrata","hog")

## 2. 去除冗余基因
Am_gene_hog <- Am_hog %>% distinct(gene, .keep_all = TRUE) %>% 
     select(gene, hog, phylostrata)
Ej_gene_hog <- Ej_hog %>% distinct(gene, .keep_all = TRUE) %>%
     select(gene, hog, phylostrata)
Go_gene_hog <- Go_hog %>% distinct(gene, .keep_all = TRUE) %>%
     select(gene, hog, phylostrata)

## 3. 年龄划分
Am_gene_hog$age <- "init"
Am_gene_hog[Am_gene_hog$phylostrata == "LUCA",]$age <- 'UC'
Am_gene_hog[Am_gene_hog$phylostrata == "Eukaryota",]$age <- 'UC'
Am_gene_hog[Am_gene_hog$phylostrata == "Viridiplantae",]$age <- 'UC'
Am_gene_hog[Am_gene_hog$age == "init",]$age <- "MC"

Ej_gene_hog$age <- "init"
Ej_gene_hog[Ej_gene_hog$phylostrata == "LUCA",]$age <- 'UC'
Ej_gene_hog[Ej_gene_hog$phylostrata == "Eukaryota",]$age <- 'UC'
Ej_gene_hog[Ej_gene_hog$phylostrata == "Viridiplantae",]$age <- 'UC'
Ej_gene_hog[Ej_gene_hog$age == "init",]$age <- "MC"

Go_gene_hog$age <- "init"
Go_gene_hog[Go_gene_hog$phylostrata == "LUCA",]$age <- 'UC'
Go_gene_hog[Go_gene_hog$phylostrata == "Eukaryota",]$age <- 'UC'
Go_gene_hog[Go_gene_hog$phylostrata == "Viridiplantae",]$age <- 'UC'
Go_gene_hog[Go_gene_hog$age == "init",]$age <- "MC"

# 4.画图
bar_data <- data.frame(species = rep(c("Artemisia_montana","Eurya_japonica","Glochidion_obovatum"),2),
                       age = c("UC","UC","UC","MC","MC","MC"),
                       count = c(8913,5075,5100,10662,7363,6997))
ggplot(bar_data, aes(x = species, y = count, fill = age)) + 
     geom_bar(stat = "identity", position = "dodge")

write.table(Am_gene_hog, "Am_hog_age.txt")
write.table(Ej_gene_hog, "Ej_hog_age.txt")
write.table(Go_gene_hog, "Go_hog_age.txt")
