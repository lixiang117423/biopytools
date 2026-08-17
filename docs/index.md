<div class="hero" markdown>
# 每个生信模块，一份零基础也能读懂的文档

BioPyTools 的 210 个工具模块，用统一模板写成使用文档——干什么、输入什么、输出什么、结果怎么看、参数怎么选，全部讲清楚。

<div class="badge-row" markdown>
<span class="chip">210 个模块</span>
<span class="chip">中英双语</span>
<span class="chip">零基础友好</span>
</div>
</div>

<div class="grid cards" markdown>

-   :material-dna: **CIM 复合区间作图**

    ---

    BSA 群体 QTL 定位：VCF + 表型 → LOD 曲线 + 峰值表，含零基础概念速览

    [:octicons-arrow-right-24: 查看文档](cim.md)

-   :material-book-open-variant: **文档怎么看**

    ---

    每页固定结构：概念速览 → 参数说明 → 结果解读 → 参数建议 → FAQ

    [:octicons-arrow-right-24: 了解结构](#structure)

-   :material-folder-multiple-outline: **全部模块去哪找**

    ---

    210 个模块按功能域分组，收纳在左侧导航栏，顶部搜索框可直达（模块名/命令名最准）

    [:octicons-arrow-right-24: 查看分组](#categories)

</div>

## 模块分类 | Categories { #categories }

全量模块不放在首页，而是按**功能域分组**收纳在左侧导航栏。分组规划如下（边界可调整，后续按此分批上线）：

| 分组<br>Group | 覆盖模块示例<br>Example modules |
|------|--------------|
| 群体遗传 | cim、admixture、fst、pixy、treemix、vcf2pca |
| 基因组组装 | hifiasm、purge_dups、yahs、ragtag、telocomp |
| 基因组注释 | braker、braker4phyto、annorefine、busco、gffcompare |
| 转录组 | rnaseq、longrnaseq、dual_rnaseq、sra2fastq |
| 变异检测 | bwa_gatk、gatk_joint、vcf_filter、annovar |
| 泛基因组 | cactus、pggb、ngenomesyn、panman |
| 三维基因组 | allhic、cphasing、haphic、hic_qc |
| 系统发育 | iqtree、raxml、mafft_fasttree |
| 蛋白质分析 | tmhmm、signalp、deeploc、interproscan |
| 重复序列 | repeatmask、edta、lai |
| 可视化 | samplot、plotsr、jcvi、hic_heatmap |
| 工具类 | fastp、seq_extract、vcf_sampler、kmertools |

## 文档结构 | What each page covers { #structure }

| 章节<br>Section | 解决什么问题<br>What it answers |
|------|--------------|
| **功能概述** | 这工具干什么、解决什么问题 |
| **零基础概念速览** | 术语的生活化解释（面向无生信基础的读者） |
| **参数说明** | 每个参数管什么、调了会怎样、什么时候需要动 |
| **结果解读** | 输出文件怎么看、好坏怎么判 |
| **参数选择建议** | 按场景给出推荐值 |
| **常见问题** | 真实踩坑点与断点续传等使用陷阱 |
