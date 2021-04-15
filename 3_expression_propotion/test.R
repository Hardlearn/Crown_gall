library(ggplot2)
gene.exp.35d = read.csv("gene.exp.35d.csv",row.names = 1)
gene.exp.6d = read.csv("gene.exp.6d.csv", row.names = 1)

# UC and MC genes annotation
gene.exp.35d$age = NA
gene.exp.35d[gene.exp.35d$level <= 3, ]$age = "UC"
gene.exp.35d[gene.exp.35d$level > 3, ]$age = "MC"

#35days
total_exp_35d = sum(gene.exp.35d$d35_gene_exp_value)
UC_total_exp_35d = sum(gene.exp.35d[gene.exp.35d$age == "UC",]$d35_gene_exp_value)
UC_propotion_35d = UC_total_exp_35d/total_exp_35d
MC_total_exp_35d = sum(gene.exp.35d[gene.exp.35d$age == "MC",]$d35_gene_exp_value)
MC_propotion_35d = MC_total_exp_35d/total_exp_35d

#6days
total_exp_6d = sum(gene.exp.6d$X6d_exp_value)
UC_total_exp_6d = sum(gene.exp.6d[gene.exp.6d$age=="UC",]$X6d_exp_value)
UC_propotion_6d = UC_total_exp_6d/total_exp_6d
MC_total_exp_6d = sum(gene.exp.6d[gene.exp.6d$age == "MC",]$X6d_exp_value)
MC_propotion_6d = MC_total_exp_6d/total_exp_6d

GV_total_6d = sum(gene.exp.6d$GV_6d_exp_value)
UC_GV = sum(gene.exp.6d[gene.exp.6d$age == "UC",]$GV_6d_exp_value)
MC_GV = sum(gene.exp.6d[gene.exp.6d$age == "MC",]$GV_6d_exp_value)


# total_exp_ref = sum(gene.exp.35d$d35_ref_gene_exp_value)
# UC_total_exp_ref = sum(gene.exp.35d[gene.exp.35d$age == "UC",]$d35_ref_gene_exp_value)
# UC_ref_propotion = UC_total_exp_ref/total_exp_ref
# MC_total_exp_ref = sum(gene.exp.35d[gene.exp.35d$age == "MC",]$d35_ref_gene_exp_value)
# MC_ref_propotion = MC_total_exp_ref/total_exp_ref

# plot
df <- data.frame(value = c(UC_total_exp_35d, MC_total_exp_35d, 
                           UC_total_exp_6d, MC_total_exp_6d),
                 age = rep(c("UC","MC"),2),
                 group = c("35days","35days","6days","6days"))

ggplot(df, aes(x=group,y=value)) +
  geom_bar(stat = 'identity',aes(fill = age)) +
  geom_text(mapping = aes(label = round(value,2),vjust=2.4,)) +
  theme_bw()



            
