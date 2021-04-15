# 1.基因年龄及划分

通过两种方法（R包`phylostratr`，OMA数据库）我们可以追溯到每个基因最早出现的时期。如下图所示我们可以看到使用`phylostratr`和OMA数据库得到拟南芥的基因年龄。 

> 文件放置在./Crown_gall/0_database_result/路径下

![Image](https://github.com/Hardlearn/Crown_gall/tree/main/IMG-Folder/image-20210415103758960.png)

![image-20210415103821912](C:\Users\13284\AppData\Roaming\Typora\typora-user-images\image-20210415103821912.png)



# 2.基因表达谱及差异基因

> 文件放置在./Crown_gall/2_diff_gene/下

## 2.1 数据来源

- 3 hours(T-DNA未整合) ： https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13929
- 6 days(T-DNA已在宿主细胞中表达):https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13930
- 35 days(tumors 已经形成):https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13927



## 2.2 数据下载

主要探究的是tumor形成过程中的基因表达情况，所以主要分析6 、35 days的数据。下载6、35天的基因表达矩阵，使用NCBI的GEO2R在线交互网站进行差异表达基因的筛选，或者直接使用`limma`包进行筛选（adj.p.value < 0.05)。



## 2.3 数据注释

根据TAIR数据库的基因编号，确定每个差异基因的年龄，如下图所示

![image-20210415112124542](C:\Users\13284\AppData\Roaming\Typora\typora-user-images\image-20210415112124542.png)

最终得到有年龄标记的差异表达基因结果为：

- 6 days: 3944个基因
- 35days: 5788个基因



## 2.4 表达量分布

> 文件放置在./Crown_gall/3_expression_propotion/下

如下图所示展示了6天和35天的基因表达量的分布情况

![6d_value_distribution](F:\project\Grown_gall\3_expression_propotion\6d_value_distribution.png)



![35d_distribution_value](F:\project\Grown_gall\3_expression_propotion\35d_distribution_value.png)





## 2.5 MC、UC基因占比

同时也统计了UC、MC基因的总表达量的在两组数据中的占比，以及基因个数的占比情况

![UC_MC_sum_value](F:\project\Grown_gall\3_expression_propotion\UC_MC_sum_value.png)

![propotion of diff genes](F:\project\Grown_gall\3_MCUCrate\propotion of diff genes.png)





上下调基因的个数也已条形图展示出来。



![6days_regulate](F:\project\Grown_gall\3_MCUCrate\6days_regulate.png)

![35days_regulate](F:\project\Grown_gall\3_MCUCrate\35days_regulate.png)



# 3.  GV菌和C58菌的差异

> 文件放置在./Crown_gall/4_GV2Tumor下

为了区分由外源菌入侵植物所造成的紊乱现象，我们比较了6天时GV（不会诱导Tumor产生的菌）和C58（诱导Tumor产生的菌）入侵时拟南芥的差异表达基因，并绘制了UC、MC基因的散点图。

![MC_plot](F:\project\Grown_gall\4_GV2Tumor\MC_plot.png)![UC_plot](F:\project\Grown_gall\4_GV2Tumor\UC_plot.png)

> Tumor: logFC Tumor P-Value < 0.05, logFC GV P-Value >= 0.05.
>
> GV : logFC Tumor P-Value >= 0.05, logFC GV P-Value < 0.05.
>
> Both: logFC Tumor P- Value < 0.05, logFC GV P-Value < 0.05



# 4. Gene ontology analysis

> 文件放置在./Crown_gall/5_GO/下

使用R包`clusterProfiler`对差异基因进行富集，并使用其中的`simplify`函数对结果去冗余，`cutoff = 0.6`，最后得到85个Biological Process Term。根据每个GOTERM的UC、MC所占比例进行Permutation test来判断Term的年龄（UCTerm or MCTerm）。结果如下图所示

![image-20210415135311822](C:\Users\13284\AppData\Roaming\Typora\typora-user-images\image-20210415135311822.png)





# 5. GOTerm的上下调

> 文件放置在./Crown_gall/6_GOTERM_FC/下

确定好GOTERM的年龄之后可以在生物学功能的基础上对UC、MC的改变模式进行观察，如下图所示

![GOTERM_logFC](F:\project\Grown_gall\6_GOTERM_FC\GOTERM_logFC.png)

> logFC: 在Tumor组中每个TERM中所包含的基因的表达总量与对照组的差异



更进一步的对每个TERM中UC基因、MC基因各自的基因总表大量改变进行了统计，如下图所示



![UCMC_gene_logFC](F:\project\Grown_gall\6_GOTERM_FC\UCMC_gene_logFC.png)
