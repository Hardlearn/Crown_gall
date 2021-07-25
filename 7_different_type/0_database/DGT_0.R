setwd("F:/project/Grown_gall/7_diff_gall_type/0_database")

library(dplyr)
library(ggplot2)
# !.数据读取
Am_hog <- read.table("Am_oma_result.1.txt", sep = " ")
colnames(Am_hog) <- c("transcript", "hog", "id", "phylostrata")
Ej_hog <- read.table("Ej_oma_result.1.txt", sep = " ")
colnames(Ej_hog) <- c("transcript", "hog", "id", "phylostrata")
Go_hog <- read.table("Go_oma_result.1.txt", sep = " ")
colnames(Go_hog) <- c("transcript", "hog", "id", "phylostrata")

Am.trans2gene <- read.table("Am.gene_trans_map.txt")
colnames(Am.trans2gene) <- c('gene', 'transcript')
Ej.trans2gene <- read.table("Ej.gene_trans_map.txt")
colnames(Ej.trans2gene) <- c('gene', 'transcript')
Go.trans2gene <- read.table("Go.gene_trans_map.txt")
colnames(Go.trans2gene) <- c('gene', 'transcript')

# 2.将转录本转化为基因
Am_hog <- merge(Am_hog, Am.trans2gene, by = 'transcript')
Ej_hog <- merge(Ej_hog, Ej.trans2gene, by = 'transcript')
Go_hog <- merge(Go_hog, Go.trans2gene, by = 'transcript')

# 3.去除冗余基因
Am_gene_hog <- Am_hog %>% distinct(gene, .keep_all = TRUE) %>% 
     select(gene, hog, phylostrata)
Ej_gene_hog <- Ej_hog %>% distinct(gene, .keep_all = TRUE) %>%
     select(gene, hog, phylostrata)
Go_gene_hog <- Go_hog %>% distinct(gene, .keep_all = TRUE) %>%
     select(gene, hog, phylostrata)

# 4.年龄划分
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

# 5.画图
bar_data <- data.frame(species = rep(c("Artemisia_montana","Eurya_japonica","Glochidion_obovatum"),2),
                       age = c("UC","UC","UC","MC","MC","MC"),
                       count = c(8912,8124,7610,10663,11812,10161))
ggplot(bar_data, aes(x = species, y = count, fill = age)) + 
     geom_bar(stat = "identity", position = "dodge")
