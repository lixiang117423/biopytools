"""
func_anno: 蛋白功能注释流水线|Protein functional annotation pipeline.

端到端: interproscan(结构域) + eggnog-mapper(GO/KEGG 源) → 标准 GO/KEGG 表(衔接下游 R).
|End-to-end: interproscan (domains) + eggnog-mapper (GO/KEGG source) → standard
GO/KEGG tables for downstream R enrichment.

约束|Constraint: 不改 interproscan/eggnog_mapper 源码, 仅 import 调用(braker4ps 模式).
|import-only (braker4ps pattern).
"""

from .table_builder import build_tables, load_go_dict
from .kegg_db import KEGGDatabase

__version__ = "1.0.0"
__all__ = ["build_tables", "load_go_dict", "KEGGDatabase"]
