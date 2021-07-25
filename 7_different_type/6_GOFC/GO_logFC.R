library(dplyr)
library(ggplot2)

ego_result = read.csv("../5_enrichGO/Am_clusterprofiler_result.csv")
tpm_matrix = read.table("../4_DE/Am.salmon.gene.TPM.not_cross_norm")
DE_result = read.table("../4_DE/AM.voom.3213776.dir/salmon.gene.counts.matrix.Am_gall_vs_Am_leaf.voom.DE_results")
GO_age = read.csv("../5_enrichGO/Am_GO_age_permutation.csv", row.names = 1)
gene_age = read.table("../2_GeneCount/Am_hog_age.txt")
colnames(gene_age)[1] = "GeneID"

#remove the NA age row
GO_age_filter = GO_age[complete.cases(GO_age),]
colnames(GO_age_filter)[1] = "ID"
ego_result = merge(ego_result, GO_age_filter[,c(1,9)], by = "ID")

#FDR < 0.01, abs(logFC) > 1
DE_result_fliter = DE_result[abs(DE_result$logFC)>1 & DE_result$FDR < 0.01,]
DE_gene_name = rownames(DE_result_fliter)

DE_tpm_matrix = tpm_matrix[DE_gene_name,]

#mean value matrix
averafe_count_matrix = data.frame(GeneID = rownames(DE_tpm_matrix),
                                  AveGall = (DE_tpm_matrix$Am_gall_rep1_TranscriptAbundance +
                                                  DE_tpm_matrix$Am_gall_rep2_TranscriptAbundance +
                                                  DE_tpm_matrix$Am_gall_rep3_TranscriptAbundance)/3,
                                  AveLeaf = (DE_tpm_matrix$Am_leaf_rep1_TranscriptAbundance +
                                                  DE_tpm_matrix$Am_leaf_rep2_TranscriptAbundance +
                                                  DE_tpm_matrix$Am_leaf_rep3_TranscriptAbundance)/3)

#enrich gene 2 GO
strsplit(ego_result$geneID[1],split = '/')[[1]][2]

item_1 = c()
item_2 = c()
for(i in 1:dim(ego_result)[1]){
     for(j in 1:ego_result$Count[i]){
          item_1 = append(item_1, ego_result$ID[i])
          item_2 = append(item_2, strsplit(ego_result$geneID[1],split = '/')[[1]][j])
     }
     
}

GO_gene_enrich = data.frame(GOID = item_1, GeneID = item_2)
GO_gene_enrich = merge(GO_gene_enrich, gene_age[,c(1,4)], by = "GeneID" )

GO_gene_enrich = merge(GO_gene_enrich, averafe_count_matrix, by = "GeneID")

#total count per GO
GO_leaf_sum <- GO_gene_enrich %>% group_by(GOID, age) %>% summarise(AveLeafSum = sum(AveLeaf))
GO_gall_sum <- GO_gene_enrich %>% group_by(GOID, age) %>% summarise(AveGallSum = sum(AveGall))

GO_leaf_gall_sum <- cbind(GO_gall_sum, GO_leaf_sum[,3])

#log2FC
GO_leaf_gall_sum$logFC <- log2(GO_leaf_gall_sum$AveGallSum + 1) - log2(GO_leaf_gall_sum$AveLeafSum + 1)

colnames(GO_leaf_gall_sum)[1] = "ID" 
GO_leaf_gall_sum = merge(GO_leaf_gall_sum, ego_result[,c(1,3)], by = "ID")
write.csv(GO_leaf_gall_sum, "Am_GO_leaf_gall_sum.csv")
