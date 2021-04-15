library(dplyr)
library(ggplot2)
## 文件读取
ego_result_indi <- read.csv("../3_Goslim/result/ego_result_indi.csv",header = FALSE)
colnames(ego_result_indi) = c("GOID","GOTERM","GeneID")
matrix_value_35d = read.table("../1_diff_gene/GSE13927_series_matrix.txt/GSE13927_series_matrix.txt", 
                              comment.char = "!", header = T)
colnames(matrix_value_35d) <- c("probe_id","Ref1","Ref2","Ref3","Ref4",
                                "Tumor1","Tumor2","Tumor3","Tumor4")
probe_tair <- read.csv("../1_diff_gene/35d.anno.gene.csv",row.names = 1)
DEG_35d <- read.csv("../1_diff_gene/35d.diff.gene.csv", row.names = 1)
EnrichGO_age <- read.csv("../3_Goslim/result/EnrichGO_age.csv", row.names = 1)
anno_UCMC <- read.csv("../1.1_MCUC_anno/anno_UCMC.csv",row.names = 1)

## 平均表达量
Average_value_35d <- data.frame(GeneID = matrix_value_35d$AGI, 
                                AveRef = (matrix_value_35d$Ref1 + matrix_value_35d$Ref2 + matrix_value_35d$Ref3 + 
                    matrix_value_35d$Ref4)/4, 
                    AveTumor = (matrix_value_35d$Tumor1 + matrix_value_35d$Tumor2 + matrix_value_35d$Tumor3 + matrix_value_35d$Tumor4)/4)



## 富集基因表达矩阵
ego_gene_value <- merge(ego_result_indi, Average_value_35d, by = "GeneID")

## 计算GOTERM的基因表达总量
ego_Ref_sum <- 
     ego_gene_value %>% group_by(GOTERM) %>% summarise(AveRefSum = sum(AveRef))

ego_Tumor_sum <-
     ego_gene_value %>% group_by(GOTERM) %>% summarise(AveTumorSum = sum(AveTumor))

ego_value_sum <- merge(ego_Ref_sum, ego_Tumor_sum, by = "GOTERM")

## 计算logFC
ego_value_sum$logFC <- log2(ego_value_sum$AveTumorSum) - log2(ego_value_sum$AveRefSum)

## 查看UC、MC GOTERM的基因表达量上下调情况
colnames(ego_value_sum)[1] = "Description"
EnrichGO_age <- merge(EnrichGO_age, ego_value_sum, by = "Description")

ggplot(data = EnrichGO_age, aes(x = logFC, y = -log10(p.adjust), color = age)) + 
     geom_point()

## 查看每个TERM种UC、MC基因的变化情况
colnames(anno_UCMC)[1] = "GeneID"
ego_gene_value <- merge(ego_gene_value, anno_UCMC, by = "GeneID")

ego_Ref_sum_UCMC <- 
     ego_gene_value %>% group_by(GOTERM, age) %>% summarise(AveRefSum = sum(AveRef))

ego_Tumor_sum_UCMC <-
     ego_gene_value %>% group_by(GOTERM, age) %>% summarise(AveTumorSum = sum(AveTumor))

ego_UCMC_logFC <- data.frame(GOTERM = ego_Ref_sum_UCMC$GOTERM, 
                             age = ego_Ref_sum_UCMC$age,
                             logFC = log2(ego_Tumor_sum_UCMC$AveTumorSum) - log2(ego_Ref_sum_UCMC$AveRefSum))

EnrichGO_UCMC <- na.omit(EnrichGO_age)
EnrichGO_UCMC <- EnrichGO_UCMC[,c(1,17)]
colnames(EnrichGO_UCMC) <- c("GOTERM","TERMAGE")
ego_UCMC_logFC <- merge(ego_UCMC_logFC, EnrichGO_UCMC, by = "GOTERM")
ego_UCMC_logFC <- ego_UCMC_logFC[order(ego_UCMC_logFC$TERMAGE,decreasing = TRUE),]
ego_UCMC_logFC$GOTERM <- factor(ego_UCMC_logFC$GOTERM, 
                                levels = unique(ego_UCMC_logFC$GOTERM) )

ggplot(ego_UCMC_logFC, aes(x=logFC, y=GOTERM, color = age)) + 
     geom_point() + 
        annotate("segment",x = 0.5, xend = 0.5, y = "anion transport", yend = "tetraterpenoid metabolic process") + 
        annotate("segment",x = 0.49, xend = 0.5, y = "anion transport", yend = "anion transport") + 
        annotate("segment",x = 0.49, xend = 0.5, y = "tetraterpenoid metabolic process", yend = "tetraterpenoid metabolic process") +
        annotate("text", x = 0.54, y = "nuclear DNA replication",label = "UC") + 
        annotate("segment", x = 0.49, xend = 0.5, y = "anatomical structure formation involved in morphogenesis", yend = "anatomical structure formation involved in morphogenesis") + 
        annotate("segment", x = 0.5, xend = 0.5, y = "anatomical structure formation involved in morphogenesis", yend = "xylem development") + 
        annotate("segment", x = 0.49, xend = 0.5, y = "xylem development", yend = "xylem development") + 
        annotate("text", x = 0.54, y = "immune system process", label = "MC")
