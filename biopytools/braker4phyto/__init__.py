"""
BRAKER4Phyto疫霉菌基因组注释工具包|BRAKER4Phyto Oomycete Genome Annotation Toolkit

功能: 基于 annorefine 的疫霉/卵菌端到端基因注释
      (BRAKER 注释 + 同源查漏补缺 → 整合 GFF3),
      与 annorefine 流程完全一致,仅默认关闭重复序列屏蔽
      (疫霉效应子多位于 repeat 区,不 mask 才能保住)
|Features: annorefine-based end-to-end annotation for Phytophthora/oomycetes
     (BRAKER + homology gap-filling → integrated GFF3). Same pipeline as
     annorefine, but repeat masking is off by default (effectors in repeats).
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-08-17

使用示例|Usage Examples:
    biopytools braker4phyto -g genome.fa -s my_phyto -p proteins.fa -o out/
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"
