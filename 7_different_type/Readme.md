# 1. 基本统计

通过OMAmer对三个植物得基因进行年龄注释，同样以绿色植物为分界线来区分UC和MC基因。结果如下图所示。

![Image](https://github.com/Hardlearn/Crown_gall/blob/main/7_different_type/0_database/bar_plot.png)

# 2. 对OMA结果进行过滤

因为同一个基因可能出现几个不同的年龄注释，对于这种情况，我们比较了每个基因转录本的表达量，选取表达量最高的isoform对应的年龄注释为最终基因的年龄注释。
>  具体看./1_fliter/basic_stat.R

# 3. 差异表达分析

使用Trinity中的脚本进行分析。具体参考：https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression

# 4. GO富集分析

对不同物种肿瘤的差异表达基因进行GO富集分析，并计算GOTERM的年龄。

> 结果及处理脚本于./5_enrichGO中

# 5. GOTERM logFC

比较能成功注释年龄的GOTERM中UC基因表达量与MC基因的表达量的logFC。

Am：

![Image](https://github.com/Hardlearn/Crown_gall/blob/main/7_different_type/6_GOFC/6_GOFC/Am_logFC.png)


Ej：差异表达基因过少，富集通路过少，无法计算GOTERM的年龄，所以没有进行此分析


Go_mature(Go_young情况与Ej一致):

![Image](https://github.com/Hardlearn/Crown_gall/blob/main/7_different_type/6_GOFC/6_GOFC/Go_mature_logFC.png)
