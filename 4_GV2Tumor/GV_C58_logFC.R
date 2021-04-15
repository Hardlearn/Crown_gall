library(ggplot2)

Tumor_diff_6d <- read.csv("F:/project/Grown_gall/1_diff_gene/6d.anno.gene.csv",row.names = 1)
GV_diff_6d <- read.csv("GV_diff_6d.csv")
MCUC <- read.csv("../1.1_MCUC_anno/anno_UCMC.csv",row.names = 1)

# GV gene
colnames(GV_diff_6d)[1] = "probe_id"
GV_diff_6d <- merge(GV_diff_6d, MCUC, by="probe_id")
Tumor_diff_6d <- merge(Tumor_diff_6d, MCUC, by="probe_id")
Tumor_diff_6d$class <- "Tumor"
GV_diff_6d$class <- "GV"

# UC gene fold change
Tumor_UC <- Tumor_diff_6d[Tumor_diff_6d$age == "UC",]
Tumor_UC <- Tumor_UC[,c(9,6,2,13)]
colnames(Tumor_UC)[1] <- "gene_id"
colnames(Tumor_UC)[2] <- "logFC.Tumor" 
colnames(Tumor_UC)[3] <- "adj.P.Val.Tumor"

GV_UC <- GV_diff_6d[GV_diff_6d$age == "UC", ]
GV_UC <- GV_UC[,c(9,6,2,10)]
colnames(GV_UC)[2] <- "logFC.GV"
colnames(GV_UC)[3] <- "adj.P.Val.GV"

UC_plot <- merge(Tumor_UC, GV_UC, by = "gene_id")
UC_plot <- UC_plot[UC_plot$adj.P.Val.Tumor < 0.05 | UC_plot$adj.P.Val.GV < 0.05,]
UC_plot$class <- NA
UC_plot[UC_plot$adj.P.Val.Tumor < 0.05 & UC_plot$adj.P.Val.GV < 0.05,]$class <- "Both"
UC_plot[UC_plot$adj.P.Val.Tumor < 0.05 & UC_plot$adj.P.Val.GV >= 0.05,]$class <- "Tumor"
UC_plot[UC_plot$adj.P.Val.Tumor >= 0.05 & UC_plot$adj.P.Val.GV < 0.05,]$class <- "GV" 

ggplot(data=UC_plot, aes(x=logFC.Tumor, y = logFC.GV, color = class)) + geom_point()

# MC gene fold change 
Tumor_MC <- Tumor_diff_6d[Tumor_diff_6d$age == "MC",]
Tumor_MC <- Tumor_MC[,c(9,6,2,13)]
colnames(Tumor_MC)[1] <- "gene_id"
colnames(Tumor_MC)[2] <- "logFC.Tumor" 
colnames(Tumor_MC)[3] <- "adj.P.Val.Tumor"
GV_MC <- GV_diff_6d[GV_diff_6d$age == "MC",]
GV_MC <- GV_MC[,c(9,6,2,10)]
colnames(GV_MC)[2] <- "logFC.GV"
colnames(GV_MC)[3] <- "adj.P.Val.GV"

MC_plot <- merge(Tumor_MC, GV_MC, by = "gene_id")
MC_plot <- MC_plot[MC_plot$adj.P.Val.Tumor < 0.05 | MC_plot$adj.P.Val.GV < 0.05,]
MC_plot$class <- NA
MC_plot[MC_plot$adj.P.Val.Tumor < 0.05 & MC_plot$adj.P.Val.GV < 0.05,]$class <- "Both"
MC_plot[MC_plot$adj.P.Val.Tumor < 0.05 & MC_plot$adj.P.Val.GV >= 0.05,]$class <- "Tumor"
MC_plot[MC_plot$adj.P.Val.Tumor >= 0.05 & MC_plot$adj.P.Val.GV < 0.05,]$class <- "GV"

ggplot(data=MC_plot, aes(x=logFC.Tumor, y=logFC.GV, color = class)) + geom_point()

