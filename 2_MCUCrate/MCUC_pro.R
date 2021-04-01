library(ggplot2)
load("F:/project/Grown_gall/1_diff_gene/.RData")

## differential expression genes propotion in each phylostrata
diff_35d = table(gene.35d$level)
all_35d = table(X35d_anno_gene$level)
ratio_35d = diff_35d/all_35d
ratio_35d = as.data.frame(ratio_35d)
ratio_35d$days = "35days"
colnames(ratio_35d)[1] = "phylostrata"
colnames(ratio_35d)[2] = "propotion"

diff_6d = table(gene.6d$level)
all_6d = table(X6d_anno_gene$level)
ratio_6d = diff_6d/all_6d
ratio_6d = as.data.frame(ratio_6d)
ratio_6d$days = "6days"
colnames(ratio_6d)[1] = "phylostrata"
colnames(ratio_6d)[2] = "propotion"

df = rbind(ratio_35d,ratio_6d)
ggplot(df,aes(x=phylostrata, y=propotion, group=days)) + 
     geom_line(aes(linetype=days, color=days)) +
     geom_point(aes(shape=days, color=days)) +
     theme_bw() +
     annotate("segment", x = 1, xend = 3, y = 0.45, yend = 0.45) +
     annotate("segment", x = 1 ,xend = 1, y = 0.443, yend = 0.45) +
     annotate("segment", x = 3 ,xend = 3, y = 0.443, yend = 0.45) +
     annotate("text",x = 2, y = 0.46, label = "UC genes") + 
     annotate("segment", x = 3, xend = 14, y = 0.45, yend = 0.45) +
     annotate("segment", x = 14, xend = 14, y = 0.443, yend = 0.45) +
     annotate("text",x = 8, y = 0.46, label = "MC genes")

## propotion of up regulate genes and down regulate genes in each phylostrata
# 35days
up_35d = gene.35d[gene.35d$logFC>0,]
up_35d = table(up_35d$level)
up_35d = as.data.frame(up_35d)
colnames(up_35d)[1] = "phylostrata"
colnames(up_35d)[2] = "count"
up_35d$days = "35days"
up_35d$regulate = "up"

down_35d = gene.35d[gene.35d$logFC<0,]
down_35d = table(down_35d$level)
down_35d = as.data.frame(down_35d)
colnames(down_35d)[1] = "phylostrata"
colnames(down_35d)[2] = "count"
down_35d$days = "35days"
down_35d$regulate = "down"

df_35d = rbind(up_35d,down_35d)

# 6days
up_6d = table(gene.6d[gene.6d$logFC>0,]$level)
up_6d = as.data.frame(up_6d)
colnames(up_6d)[1] = "phylostrata"
colnames(up_6d)[2] = "count"
up_6d$days = "6days"
up_6d$regulate = "up"

down_6d = table(gene.6d[gene.6d$logFC<0,]$level)
down_6d = as.data.frame(down_6d)
colnames(down_6d)[1] = "phylostrata"
colnames(down_6d)[2] = "count"
down_6d$days = "6days"
down_6d$regulate = "down"

df_6d = rbind(up_6d,down_6d)

# plot 
df = rbind(df_6d,df_35d)
ggplot(df_35d, aes(x=phylostrata, y = count)) + 
     geom_bar(stat = 'identity', aes(fill = regulate)) +
     geom_text(stat = 'identity', aes(label=count), position=position_stack(0.5)) +
     theme_bw()

ggplot(df_6d, aes(x=phylostrata, y = count)) + 
     geom_bar(stat = 'identity', aes(fill = regulate)) +
     geom_text(stat = 'identity', aes(label=count), position=position_stack(0.5)) +
     theme_bw()




