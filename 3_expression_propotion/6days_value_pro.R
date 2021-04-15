library(affy)
library(ggplot2)

# 1.reading the .cel files
setwd("F:/project/Grown_gall/6d/GSE13930_RAW/GSE13930_RAW")

data.raw = ReadAffy()
sampleNames(data.raw) <- c("GV_6d_3","GV_6d_2","GV_6d_1","6d_3","6d_2","6d_1",
                           "6d_ref_3","6d_ref_2","6d_ref_1")

sampleNames(data.raw)
head(geneNames(data.raw))

# 2.normalization
eset.rma = rma(data.raw)

# 3.calculate gene expression
emat.rma.log2 = exprs(eset.rma)
head(emat.rma.log2,1)

# 4.calculate the mean value
result.rma = data.frame((emat.rma.log2[,c(1,4,7)] + emat.rma.log2[,c(2,5,8)] + 
                              emat.rma.log2[,c(3,6,9)])/3)

colnames(result.rma)[1] = "GV_6d_exp_value"
colnames(result.rma)[2] = "6d_exp_value"
colnames(result.rma)[3] = "6d_ref_exp_value"
write.csv(result.rma, "rawdata_6d.csv")

# 5. the total expression value of UC and MC genes
gene.anno.6d = read.csv("x6d_anno_gene.csv", row.names = 1)
rawdata_6d = read.csv("rawdata_6d.csv")
colnames(rawdata_6d)[1] = "probe_id"

gene.anno.6d$age = NA
gene.anno.6d[gene.anno.6d$level == 1 | gene.anno.6d$level == 2 | gene.anno.6d$level == 3,]$age = "UC"
gene.anno.6d[is.na(gene.anno.6d$age), ]$age = "MC"

gene.exp.6d = merge(gene.anno.6d, rawdata_6d, by="probe_id")
write.csv(gene.exp.6d, "gene.exp.6d.csv")

# 6. plot the distribution of the gene expression value
rawdata_6d = read.csv("rawdata_6d.csv")
GV = as.data.frame(rawdata_6d$GV_6d_exp_value)
colnames(GV)[1] = "value"
GV$Class = "GV"
Tumor = as.data.frame(rawdata_6d$X6d_exp_value)
colnames(Tumor)[1] = "value"
Tumor$Class = "Tumor"
Normal = as.data.frame(rawdata_6d$X6d_ref_exp_value)
colnames(Normal)[1] = "value"
Normal$Class = "Normal"

data = rbind(GV, Tumor, Normal)

ggplot(data, aes(x=value, fill=Class)) + 
     geom_density(alpha=0.4) +
     labs(title="The ditribution of 6days gene expression value",
          x="Expression Value", y = "Density") +
     theme_classic()
