library(dplyr)
library(GO.db)
library(clusterProfiler)


Ej_GOTERM = read.table("../3_GOannotation/Ej_GOTERM.txt")
Ej_pro2gene = read.table("../3_GOannotation/Ej_Pro2Gene.txt")
colnames(Ej_pro2gene) = c('gene','protein')
colnames(Ej_GOTERM)[2] = "protein"

Ej_GOTERM = merge(Ej_GOTERM, Ej_pro2gene, by = "protein")

Ej_GOTERM = Ej_GOTERM %>% distinct(GOID,Ont,gene,.keep_all = TRUE) ##move same row

#BP
Ej_GOTERM_BP = Ej_GOTERM[Ej_GOTERM$Ont == "biological_process ",]
Ej_gomap_BP = data.frame(GOID = Ej_GOTERM_BP$GOID, gene = Ej_GOTERM_BP$gene)
Ej_gomap_BP = buildGOmap(Ej_gomap_BP)
Ej_goname_BP <- data.frame(GO=Ej_GOTERM_BP$GOID, name=Ej_GOTERM_BP$description)

#DE gene

DE_result <- read.table("../4_DE/Ej.voom.3310017.dir/salmon.gene.counts.matrix.Ej_gall_vs_Ej_leaf.voom.DE_results")
Ej_TPM_matrix = read.table("../4_DE/salmon.gene.TPM.not_cross_norm")

#FDR < 0.1, abs(logFC) > 1, TPM > 3 per sample
DE_result_fliter = DE_result[abs(DE_result$logFC)>1 & DE_result$FDR < 0.1,]
DE_gene_name = rownames(DE_result_fliter)


#enrich
ego_BP <- enricher(DE_gene_name,TERM2GENE = Ej_gomap_BP, TERM2NAME=Ej_goname_BP)
ego_BP@ontology <- "BP"
ego_BP_slim = simplify(ego_BP,cutoff=0.7, by="p.adjust", select_fun=min)
result <- ego_BP_slim@result

write.csv(result, "Ej_clusterprofiler_result.csv")
write.table(Ej_gomap_BP, "Ej_gomap_BP.txt")
