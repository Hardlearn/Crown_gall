ara.oma <- read.csv("../1_diff_gene/arabidopsis.oma.csv", row.names = 1)
probe_id <- read.csv("../1_diff_gene/35d.anno.gene.csv", row.names = 1)

colnames(probe_id)[9] = "gene_id"
ara.oma$age = NA
ara.oma[ara.oma$level == 1 | ara.oma$level == 2 | ara.oma$level == 3,]$age = "UC"
ara.oma[is.na(ara.oma$age), ]$age = "MC"

UCMC_anno <- merge(ara.oma[,c(5,6)], probe_id[,c(1,9)], by="gene_id")
write.csv(UCMC_anno, "anno_UCMC.csv")


