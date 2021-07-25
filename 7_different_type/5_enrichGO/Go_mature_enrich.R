library(dplyr)
library(GO.db)
library(clusterProfiler)

Go_GOTERM = read.table("../3_GOannotation/Go_GOTERM.txt")
Go_pro2gene = read.table("../3_GOannotation/Go_Pro2Gene.txt")
colnames(Go_pro2gene) = c('gene','protein')

Go_GOTERM = merge(Go_GOTERM, Go_pro2gene, by = "protein")

Go_GOTERM = Go_GOTERM %>% distinct(GOID,Ont,gene,.keep_all = TRUE) ##move same row

#BP
Go_GOTERM_BP = Go_GOTERM[Go_GOTERM$Ont == "biological_process ",]
Go_gomap_BP = data.frame(GOID = Go_GOTERM_BP$GOID, gene = Go_GOTERM_BP$gene)
Go_gomap_BP = buildGOmap(Go_gomap_BP)
Go_goname_BP <- data.frame(GO=Go_GOTERM_BP$GOID, name=Go_GOTERM_BP$description)

#DE gene

DE_result <- read.table("../4_DE/Go.mature.voom.3310496.dir/Go_matrue.salmon.gene.counts.matrix.Go_mature_gall_vs_Go_mature_leaf.voom.DE_results")


#FDR < 0.01, abs(logFC) > 1, TPM > 3 per sample
DE_result_fliter = DE_result[abs(DE_result$logFC)>1 & DE_result$FDR < 0.01,]
DE_gene_name = rownames(DE_result_fliter)


#enrich
ego_BP <- enricher(DE_gene_name,TERM2GENE = Go_gomap_BP, TERM2NAME=Go_goname_BP)
ego_BP@ontology <- "BP"
ego_BP_slim = simplify(ego_BP,cutoff=0.8, by="p.adjust", select_fun=min)
result <- ego_BP_slim@result

write.csv(result, "Go_mature_clusterprofiler_result.csv")
write.table(Go_gomap_BP, "Go_gomap_BP.txt")
