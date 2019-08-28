# Introduction

- Single-cell sequencing provide more information about single-cell transcriptome and chromatin
accessibility features.

```
map = c(
    "Megakaryocyte-Erythroid_precursor_cell" = "MEP",
    "Erythroid_proge`nitor" = "EPC",
    "Erythroid_progenitor_like" = "EPC-like",
    "Erythroid_Bpgm_high" = "ERY1",
    "Erythroid_Slfn14_high" = "ERY2",
    "Macrophage" = "MAC",
    "Hepatoblast" = "HEP",
    "Endothelial_cell" = "END"
)
``` 

尊敬的各位老师大家好

我是来自生物科学系的赵洵，我的科研实践II的内容是胚胎器官发育中单细胞测序分析的应用，指导老师是北京基因组研究所的刘江老师。

我将从研究背景、方法、结果以及未来计划四个方面进行汇报。

首先是研究背景，核酸测序技术被广泛应用在各种研究领域中，但是过去对细胞群体测序并不能很好的体现出细胞与细胞之间的异质性，同时例如胚胎器官发育的研究中，生物样本过于缺乏并不适合细胞群体测序，所以出现了基于高通量测序技术，又具有高分辨率低样本量特点的单细胞测序技术，它能获得各个细胞的各个特征的值，例如RNA测序中的各个基因表达量、ATAC测序中染色质区域的开放程度等等，最终以矩阵的形式表达。

我所做的主要工作是对测序数据进行分析，利用的工具是一个名为Seurat的R语言工具包。它可以读取测序后的矩阵数据，进行过滤、降维、聚类等统计分析。在聚类后，一方面可以直观的看到聚类的可视化结果，另一方面也可以通过热图的形式了解不同类别所具有的特征，经过查阅文献或查询其他数据库，将特征和细胞类型对应起来，综合后得到每个类别的细胞类型。

在这里我以小鼠着床13.5天的胚胎肝脏数据为例，已经有文献证明在13.5天时，肝脏是胚胎中一个主要的造血器官。在这幅图当中，大部分细胞被注释为血细胞类群，而与成熟肝脏相关的肝细胞则只占很少一部分。另一方面，这样的RNA测序数据表明可能存在的细胞类型变化过程