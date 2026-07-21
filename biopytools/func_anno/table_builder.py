"""
解析 eggnog .emapper.annotations → 标准 GO/KEGG 表(衔接下游 R 富集).
|Parse eggnog .emapper.annotations → standard GO/KEGG tables for downstream R enrichment.

输出严格 4 列 TSV(纯表头, 无中文, 直接喂 R 函数)|Output strict 4-column TSV
(bare headers, no Chinese, ready for R functions):
    GO 表|GO table:   gene  go_id  go_term  go_ontology
    KEGG 表|KEGG table: gene  kegg_id  kegg_term  kegg_category

数据来源|Data sources:
    GO: eggnog GOs 列 + go_data.py(GO id → name + BP/MF/CC 缩写, 复用 interproscan 内置表)
    KEGG: eggnog KEGG_Pathway 列 + kegg_data.py(ko id → name + category)
"""

import csv
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .kegg_db import KEGGDatabase

GO_HEADER = ["gene", "go_id", "go_term", "go_ontology"]
KEGG_HEADER = ["gene", "kegg_id", "kegg_term", "kegg_category"]


def load_go_dict() -> Dict[str, Tuple[str, str]]:
    """
    直接按路径加载 interproscan/go_data.py 的 GO_DATABASE(绕过 interproscan/__init__,
    避免触发 pandas 等重依赖)|Load GO_DATABASE from interproscan/go_data.py by path,
    bypassing interproscan/__init__ to avoid heavy deps (pandas, etc.).

    Returns:
        {go_id: (name, ontology_缩写|BP/MF/CC)}
    """
    go_data_path = Path(__file__).parent.parent / "interproscan" / "go_data.py"
    spec = importlib.util.spec_from_file_location("_func_anno_go_data", go_data_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.GO_DATABASE


def parse_annotations(path: str) -> Tuple[List[str], List[List[str]]]:
    """
    解析 .emapper.annotations|Parse annotations file.

    Returns:
        (headers, rows): headers 来自 #query 行; rows 为数据行(按列拆分).
    """
    headers: List[str] = []
    rows: List[List[str]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#query"):
                headers = line.lstrip("#").split("\t")
            elif line.startswith("#"):
                continue
            else:
                rows.append(line.split("\t"))
    return headers, rows


def _col_index(headers: List[str], *names: str) -> int:
    """按列名找索引(支持多候选名, 抗 eggnog 版本变化)|Find column index by name."""
    name_set = set(names)
    for i, h in enumerate(headers):
        if h in name_set:
            return i
    return -1


def _cell(row: List[str], idx: int) -> str:
    """安全取单元格|Safe cell access."""
    return row[idx] if 0 <= idx < len(row) else ""


def build_go_table(headers: List[str], rows: List[List[str]],
                   go_dict: Dict[str, Tuple[str, str]]) -> Tuple[List[tuple], int]:
    """
    构建 GO 长表|Build GO long table.

    Returns:
        (records, miss_term): records 每行 (gene, go_id, go_term, go_ontology);
        miss_term 为查不到 term 的 GO id 数|count of GO ids with no term.
    """
    idx_query = _col_index(headers, "query")
    idx_go = _col_index(headers, "GOs", "GO_terms")
    records: List[tuple] = []
    miss_term = 0

    for r in rows:
        gene = _cell(r, idx_query)
        go_str = _cell(r, idx_go)
        if not go_str or go_str == "-":
            continue
        for go_id in go_str.split(","):
            go_id = go_id.strip()
            # 去掉 eggnog 可能附加的来源标记, 如 "GO:0003700(PANTHER)"
            go_id = go_id.split("(")[0].strip()
            if not go_id.startswith("GO:"):
                continue
            name, onto = go_dict.get(go_id, ("", ""))
            if not name:
                miss_term += 1
            records.append((gene, go_id, name, onto))
    return records, miss_term


def build_kegg_table(headers: List[str], rows: List[List[str]],
                     kegg_db: KEGGDatabase) -> Tuple[List[tuple], int]:
    """
    构建 KEGG 长表|Build KEGG long table.

    Returns:
        (records, miss_term): records 每行 (gene, kegg_id, kegg_term, kegg_category);
        miss_term 为查不到 term 的 pathway 数.
    """
    idx_query = _col_index(headers, "query")
    idx_pw = _col_index(headers, "KEGG_Pathway")
    records: List[tuple] = []
    miss_term = 0

    for r in rows:
        gene = _cell(r, idx_query)
        pw_str = _cell(r, idx_pw)
        if not pw_str or pw_str == "-":
            continue
        for kid in pw_str.split(","):
            kid = kid.strip()
            if not kid or kid == "-":
                continue
            info = kegg_db.get_pathway_info(kid)
            if not info["name"]:
                miss_term += 1
            records.append((gene, kid, info["name"], info["category"]))
    return records, miss_term


def write_tsv(path: str, header: List[str], rows: List[tuple]):
    """写严格 4 列 TSV(无 BOM, \\n 换行)|Write strict 4-column TSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build_tables(
    annotations_file: str,
    output_dir: str,
    sample: str,
    go_dict: Optional[Dict[str, Tuple[str, str]]] = None,
    kegg_db: Optional[KEGGDatabase] = None,
    logger=None,
) -> dict:
    """
    主入口: 解析 eggnog annotations → 写 GO.tsv + KEGG.tsv|Main entry.

    Args:
        annotations_file: eggnog .emapper.annotations 路径.
        output_dir: 03_tables 输出目录(已存在).
        sample: 样本名(输出文件前缀).
        go_dict: 预加载的 GO 映射( None 则自动 load_go_dict).
        kegg_db: 预加载的 KEGGDatabase(None 则自动构造).

    Returns:
        统计 dict|stats dict.
    """
    if go_dict is None:
        go_dict = load_go_dict()
    if kegg_db is None:
        kegg_db = KEGGDatabase(logger=logger)

    headers, rows = parse_annotations(annotations_file)
    if not headers:
        raise ValueError(f"未找到注释表头|No #query header in: {annotations_file}")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # GO 表|GO table
    go_recs, go_miss = build_go_table(headers, rows, go_dict)
    go_path = out / f"{sample}.go.tsv"
    write_tsv(str(go_path), GO_HEADER, go_recs)

    # KEGG 表|KEGG table
    kegg_recs, kegg_miss = build_kegg_table(headers, rows, kegg_db)
    kegg_path = out / f"{sample}.kegg.tsv"
    write_tsv(str(kegg_path), KEGG_HEADER, kegg_recs)

    stats = {
        "go_rows": len(go_recs),
        "go_miss_term": go_miss,
        "kegg_rows": len(kegg_recs),
        "kegg_miss_term": kegg_miss,
        "go_table": str(go_path),
        "kegg_table": str(kegg_path),
    }
    if logger:
        logger.info(
            f"GO 表|GO table: {stats['go_rows']} 行|rows "
            f"(term 缺失|missing term: {stats['go_miss_term']}) → {go_path}"
        )
        logger.info(
            f"KEGG 表|KEGG table: {stats['kegg_rows']} 行|rows "
            f"(term 缺失|missing term: {stats['kegg_miss_term']}) → {kegg_path}"
        )
    return stats
