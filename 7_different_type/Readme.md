# 1. 基本统计

通过OMAmer对三个植物得基因进行年龄注释，同样以绿色植物为分界线来区分UC和MC基因。结果如下图所示。

![Image](https://github.com/Hardlearn/Crown_gall/blob/main/7_different_type/0_database/bar_plot.png)

# 2. 对OMA结果进行过滤

因为同一个基因可能出现几个不同的年龄注释，对于这种情况，我们比较了每个基因转录本的表达量，选取表达量最高的isoform对应的年龄注释为最终基因的年龄注释。
>  具体看./1_fliter/basic_stat.R

# 3. 差异表达分析

使用Trinity中的脚本进行分析。具体参考：https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression
