library(clusterProfiler)
library(dplyr)
library(GO.db)

# input data
ATH_GOTERM = read.table("./ATH_GO_GOSLIM.txt/ATH_GO_TERM.txt", sep="\t")
colnames(ATH_GOTERM) = c("geneID","GOTerm","Ont", "GOSlim","Source")
GOSLIMS = read.table("./ATH_GO_GOSLIM.txt/GOslim.txt", sep = "\t", header = TRUE)
colnames(GOSLIMS) = c("Ont","GO","GOslim_ID","GOSlim")

# GO to GOslim
ATH_GOTERM <- merge(GOSLIMS[,c(3,4)], ATH_GOTERM, by="GOSlim")

# delete the rebundant terms
ND_Term = which(ATH_GOTERM$Source == "ND")
ATH_GOTERM = ATH_GOTERM[-ND_Term,]
ATH_GOTERM = ATH_GOTERM %>% distinct(GOTerm,geneID,Ont,GOSlim,.keep_all = TRUE)


# BP
ATH_GOTERM_BP <- ATH_GOTERM[ATH_GOTERM$Ont == "P",] 
gomap_BP <- data.frame(GOslim=ATH_GOTERM_BP$GOslim_ID, gene=ATH_GOTERM_BP$geneID)
gomap_BP <- buildGOmap(gomap_BP) #间接注释
colnames(gomap_BP)[1] = "GOID"
# MF
ATH_GOTERM_MF <- ATH_GOTERM[ATH_GOTERM$Ont == "F",]
gomap_MF <- data.frame(GOslim=ATH_GOTERM_MF$GOslim_ID, gene=ATH_GOTERM_MF$geneID)
gomap_MF <- buildGOmap(gomap_MF) #间接注释
colnames(gomap_MF)[1] = "GOID"
# CC
ATH_GOTERM_CC <- ATH_GOTERM[ATH_GOTERM$Ont == "C",]
gomap_CC <- data.frame(GOslim=ATH_GOTERM_CC$GOslim_ID, gene=ATH_GOTERM_CC$geneID)
gomap_CC <- buildGOmap(gomap_CC) #间接注释
colnames(gomap_MF)[1] = "GOID"

'''
goname_BP <- AnnotationDbi::select(x=GO.db, keys = gomap_BP$GOslim,  keytype = "GOID",columns = "TERM" )
goname_MF <- AnnotationDbi::select(x=GO.db, keys = gomap_MF$GOslim,  keytype = "GOID",columns = "TERM" )
goname_CC <- AnnotationDbi::select(x=GO.db, keys = gomap_CC$GOslim,  keytype = "GOID",columns = "TERM" )
'''

goname_BP <- data.frame(GOID=ATH_GOTERM_BP$GOslim_ID, name=ATH_GOTERM_BP$GOSlim)
goname_MF <- data.frame(GOID=ATH_GOTERM_MF$GOslim_ID, name=ATH_GOTERM_MF$GOSlim)
goname_CC <- data.frame(GOID=ATH_GOTERM_CC$GOslim_ID, name=ATH_GOTERM_CC$GOSlim)

# enricher delete the unmaping gene 
DEG_GENE_ID = read.csv("../1_diff_gene/35d.diff.gene.csv",row.names = 1) #diff gene
colnames(DEG_GENE_ID)[1] = "geneID"

x = ATH_GOTERM %>% distinct(geneID, .keep_all = TRUE)
DEG_GENE_ID = merge(x, DEG_GENE_ID, by="geneID")
DEG_GENE_ID <- DEG_GENE_ID[DEG_GENE_ID$logFC > 0,] #上调
DEG_GENE_ID <- DEG_GENE_ID[DEG_GENE_ID$logFC < 0,] #下调
ego_BP <- enricher(DEG_GENE_ID$geneID,TERM2GENE = gomap_BP, TERM2NAME=goname_BP)
ego_MF <- enricher(DEG_GENE_ID$geneID,TERM2GENE = gomap_MF, TERM2NAME=goname_MF)
ego_CC <- enricher(DEG_GENE_ID$geneID,TERM2GENE = gomap_CC, TERM2NAME=goname_CC)

# remove redundant
ego_BP@ontology <- "BP" 
simplify(ego_BP)
ego_MF@ontology <- "MF"
simplify(ego_MF)
ego_CC@ontology <- "CC"
simplify(ego_CC)

# plot
dotplot(ego_BP)
dotplot(ego_MF)
dotplot(ego_CC)

result_BP = ego_BP@result[ego_BP@result$p.adjust<0.05,]
result_MF = ego_MF@result[ego_MF@result$p.adjust<0.05,]
result_CC = ego_CC@result[ego_CC@result$p.adjust<0.05,]
result_BP$Ont <- "BP"
result_MF$Ont <- "MF"
result_CC$Ont <- "CC" 
result = rbind(result_BP,result_CC,result_MF)
write.csv(result, "./result/down_GOslim.csv")

# MC_UC
result_GOslim = read.csv("./result/result_GOslim_indi.csv",header = F)
colnames(result_GOslim) <- c("ID","GOslim","gene_id","Ont")
result_GOslim_MCUC = merge(result_GOslim, ara.oma, by="gene_id")
write.csv(result_GOslim_MCUC, "./result/result_GOslim_MCUC.csv")
