library(dplyr)
library(GO.db)
library(clusterProfiler)


Am_GOTERM = read.table("../3_GOannotation/Am_GOTERM.txt")
AM_pro2gene = read.table("../3_GOannotation/Am_Pro2Gene.txt")
colnames(AM_pro2gene) = c('gene','protein')

Am_GOTERM = merge(Am_GOTERM, AM_pro2gene, by = "protein")

Am_GOTERM = Am_GOTERM %>% distinct(GOID,Ont,gene,.keep_all = TRUE) ##move same row

#BP
Am_GOTERM_BP = Am_GOTERM[Am_GOTERM$Ont == "biological_process ",]
Am_gomap_BP = data.frame(GOID = Am_GOTERM_BP$GOID, gene = Am_GOTERM_BP$gene)
Am_gomap_BP = buildGOmap(Am_gomap_BP)
Am_goname_BP <- data.frame(GO=Am_GOTERM_BP$GOID, name=Am_GOTERM_BP$description)

#DE gene

DE_result <- read.table("../4_DE/AM.voom.3213776.dir//salmon.gene.counts.matrix.Am_gall_vs_Am_leaf.voom.DE_results")
Am_TPM_matrix = read.table("../4_DE/salmon.gene.TPM.not_cross_norm")

#FDR < 0.01, abs(logFC) > 1, TPM > 3 per sample
DE_result_fliter = DE_result[abs(DE_result$logFC)>1 & DE_result$FDR < 0.01,]
DE_gene_name = rownames(DE_result_fliter)


#enrich
ego_BP <- enricher(DE_gene_name,TERM2GENE = Am_gomap_BP, TERM2NAME=Am_goname_BP)
ego_BP@ontology <- "BP"
ego_BP_slim = simplify(ego_BP,cutoff=0.6, by="p.adjust", select_fun=min)
result <- ego_BP_slim@result

write.csv(result, "Am_clusterprofiler_result.csv")
write.table(Am_gomap_BP, "Am_gomap_BP.txt")
