library(clusterProfiler)
library(org.At.tair.db)
library(GO.db)

DEG <- read.csv("../1_diff_gene/35d.diff.gene.csv", row.names = 1)

ego_BP <- enrichGO(DEG$AGI, OrgDb = org.At.tair.db, ont = "BP", keyType = "TAIR")
ego_BP <- simplify(ego_BP,cutoff=0.6, by="p.adjust", select_fun=min)
result <- ego_BP@result

write.csv(result, "./result/ego_result.csv")
