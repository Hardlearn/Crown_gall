#loading 6.GO analysis in workflow
X35d_anno_gene = read.csv("35d.anno.gene.csv",row.names = 1)
X6d_anno_gene = read.csv("6d.anno.gene.csv", row.names = 1)

# phylostrata level 
ara.oma = read.csv("arabidopsis.oma.csv",row.names = 1)
colnames(ara.oma)[5] = "AGI"
X35d_anno_gene = merge(X35d_anno_gene, ara.oma, by="AGI")
X6d_anno_gene = merge(X6d_anno_gene, ara.oma, by="AGI")

# filter the rows which didn't annotate.
library(dplyr)
X35d_anno_gene = filter(X35d_anno_gene, !is.na(Gene.symbol))
X6d_anno_gene = filter(X6d_anno_gene, !is.na(Gene.symbol))
X35d_anno_gene = X35d_anno_gene %>% distinct(AGI, .keep_all = TRUE)
X6d_anno_gene = X6d_anno_gene %>% distinct(AGI, .keep_all = TRUE)

# adj.p.value < 0.05
gene.35d = filter(X35d_anno_gene, adj.P.Val < 0.05)
gene.6d = filter(X6d_anno_gene, adj.P.Val < 0.05)

# UC_MC gene
gene.35d = gene.35d[,c(1,2,3,4,7,13,14)]
gene.6d = gene.6d[,c(1,2,3,4,7,13,14)]

gene.35d$age = NA
gene.35d[gene.35d$level == 1 | gene.35d$level == 2 | gene.35d$level == 3,]$age = "UC"
gene.35d[is.na(gene.35d$age), ]$age = "MC"

gene.6d$age = NA
gene.6d[gene.6d$level ==1 | gene.6d$level ==2 | gene.6d$level == 3,]$age = 'UC'
gene.6d[is.na(gene.6d$age), ]$age = "MC"

write.csv(gene.35d,file = "35d.diff.gene.csv")
write.csv(gene.6d,file = "6d.diff.gene.csv")
