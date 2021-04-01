library(affy)
library(ggplot2)
# 1.reading the .cel files
setwd("F:/project/Grown_gall/35d/GSE13927_RAW/GSE13927_RAW")

data.raw = ReadAffy()
sampleNames(data.raw) <- c("35d_4","35d_3","35d_2","35d_1","35d_ref_4",
                           "35d_ref_3","35d_ref_2","35d_ref_1")
sampleNames(data.raw)
head(geneNames(data.raw))

# 2.normalization
eset.rma = rma(data.raw)

# 3.calculate gene expression
emat.rma.log2 = exprs(eset.rma)
head(emat.rma.log2,1)

# 4.calculate the mean value
result.rma = data.frame((emat.rma.log2[,c(1,5)] + emat.rma.log2[,c(2,6)]
                         + emat.rma.log2[,c(3,7)] + emat.rma.log2[,c(4,8)])/4)
colnames(result.rma)[1] = "d35_gene_exp_value"
colnames(result.rma)[2] = "d35_ref_gene_exp_value"
write.csv(result.rma,"rawdata_35d.csv")

# 5.plot of total expression distribution
Tumor = as.data.frame(result.rma$d35_gene_exp_value)
Tumor$Class = "Tumor"
colnames(Tumor)[1] = "value"
Normal = as.data.frame(result.rma$d35_ref_gene_exp_value)
Normal$Class = "Normal"
colnames(Normal)[1] = "value"
data = rbind(Tumor,Normal)

ggplot(data, aes(x=value, fill=Class)) + 
     geom_density(alpha=0.7) +
     labs(title="The ditribution of 35days gene expression value",
          x="Expression Value", y = "Density") +
     scale_color_brewer(palette="Paired") +
     scale_fill_brewer(palette="Blues") +
     theme_classic()

# 6. the total expression value of UC and MC genes
gene.anno.35d = read.csv("x35d_anno_gene.csv", row.names = 1)
rawdata_35d = read.csv("rawdata_35d.csv")
colnames(rawdata_35d)[1] = "probe_id"

gene.anno.35d$age = NA
gene.anno.35d[gene.anno.35d$level == 1 | gene.anno.35d$level == 2 | gene.anno.35d$level == 3,]$age = "UC"
gene.anno.35d[is.na(gene.anno.35d$age), ]$age = "MC"

gene.exp.35d = merge(gene.anno.35d, rawdata_35d, by="probe_id")

