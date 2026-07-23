"""
解析 eggnog .emapper.annotations → 标准 GO/KEGG 表(衔接下游 R 富集).
|Parse eggnog .emapper.annotations → standard GO/KEGG tables for downstream R enrichment.

输出严格 4 列 TSV|Output strict 4-column TSV:
    GO 表|GO table:   gene  go_id  go_term  go_ontology
    KEGG 表|KEGG table: gene  kegg_id  kegg_term  kegg_category

KEGG 过滤(去除 eggnog KO 注释误挂的人类/动物通路)|KEGG filtering:
    1. name 黑名单(默认 DEFAULT_PLANT_EXCLUDE, 植物无关词 cancer/estrogen 等)
    2. category 黑名单(可选, 如 Human Diseases 大类; 默认空, 因整块排会误伤植物同源如 P450)
"""

import csv
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .kegg_db import KEGGDatabase

GO_HEADER = ["gene", "go_id", "go_term", "go_ontology"]
KEGG_HEADER = ["gene", "kegg_id", "kegg_term", "kegg_category"]

# 植物无关的纯人类/动物通路 name 关键词(子串匹配, 小写)|plant-irrelevant keywords.
# ⚠️ 不含 chemical carcinogenesis/drug metabolism/xenobiotics —— 植物有 P450/GST/UGT 同源.
# |Excludes P450/drug-metabolism terms (plants have homologs).
DEFAULT_PLANT_EXCLUDE = (
    "cancer,carcinoma,leukemia,lymphoma,melanoma,sarcoma,hepatoma,myeloma,glioma,"
    "cholangiocarcinoma,gastric acid secretion,"
    "atherosclerosis,shear stress,myocardial,cardiomyopathy,hypertrophic,"
    "viral myocarditis,vascular smooth muscle,"
    "estrogen,androgen,progesterone,prolactin,oxytocin,cortisol,glucagon,aldosterone,"
    "relaxin,gnrh,thyroid hormone,insulin resistance,diabetes mellitus,maturity onset,"
    "long-term potentiation,long-term depression,dopaminergic,serotonergic,axon guidance,"
    "synaptic vesicle,endocannabinoid,neurotrophin,circadian entrainment,"
    "t cell,b cell,antigen processing,graft-versus-host,natural killer,intestinal immune,"
    "hematopoietic,fc gamma,complement and coagulation,"
    "hepatitis,influenza,epstein-barr,herpes,tuberculosis,malaria,leishmaniasis,chagas,"
    "amoebiasis,staphylococcal,vibrio,pertussis,legionellosis,leptospirosis,"
    "platinum drug,drug resistance,antifolate resistance,egfr,chemokine,"
    "alzheimer,parkinson,huntington,prion,amyotrophic,adderall,cocaine,morphine,"
    "renal cell carcinoma,prostate,bladder cancer,pancreatic,thyroid cancer"
)


def load_go_dict() -> Dict[str, Tuple[str, str]]:
    """
    直接按路径加载 interproscan/go_data.py 的 GO_DATABASE(绕过 interproscan/__init__,
    避免触发 pandas 等重依赖)|Load GO_DATABASE by path, bypassing interproscan/__init__.

    Returns:
        {go_id: (name, ontology_缩写|BP/MF/CC)}
    """
    go_data_path = Path(__file__).parent.parent / "interproscan" / "go_data.py"
    spec = importlib.util.spec_from_file_location("_func_anno_go_data", go_data_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.GO_DATABASE


def parse_annotations(path: str) -> Tuple[List[str], List[List[str]]]:
    """解析 .emapper.annotations|Parse annotations file. Returns (headers, rows)."""
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
    """按列名找索引(支持多候选名)|Find column index by name."""
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
    """构建 GO 长表|Build GO long table. Returns (records, miss_term)."""
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
            go_id = go_id.split("(")[0].strip()
            if not go_id.startswith("GO:"):
                continue
            name, onto = go_dict.get(go_id, ("", ""))
            if not name:
                miss_term += 1
            records.append((gene, go_id, name, onto))
    return records, miss_term


def _name_excluded(name: str, keyword_list) -> bool:
    """通路 name 含任一关键词(小写子串)→ 排除|Exclude if name contains keyword."""
    if not keyword_list or not name:
        return False
    name_lower = name.lower()
    return any(kw and kw.lower() in name_lower for kw in keyword_list)


def _category_excluded(cat_a: str, cat_b: str, exclude_list) -> bool:
    """category_A/B 任一包含排除关键词 → True(子串匹配)|Exclude if A/B matches keyword."""
    if not exclude_list:
        return False
    combined = f"{cat_a}; {cat_b}"
    return any(ex and ex in combined for ex in exclude_list)


def build_kegg_table(headers: List[str], rows: List[List[str]],
                     kegg_db: KEGGDatabase,
                     exclude_keywords=None, exclude_categories=None
                     ) -> Tuple[List[tuple], int, int]:
    """
    构建 KEGG 长表|Build KEGG long table.

    Args:
        exclude_keywords: name 黑名单关键词列表(子串匹配, 小写).
        exclude_categories: category A/B 黑名单关键词列表(子串匹配).

    Returns:
        (records, miss_term, excluded): records 每行 (gene, kegg_id, kegg_term, kegg_category_A).
    """
    idx_query = _col_index(headers, "query")
    idx_pw = _col_index(headers, "KEGG_Pathway")
    records: List[tuple] = []
    miss_term = 0
    excluded = 0

    for r in rows:
        gene = _cell(r, idx_query)
        pw_str = _cell(r, idx_pw)
        if not pw_str or pw_str == "-":
            continue
        for kid in pw_str.split(","):
            kid = kid.strip()
            if not kid or kid == "-":
                continue
            # eggnog KEGG_Pathway 同一通路冗余给 ko/map; 用户选 ko → 只留 ko
            if not kid.startswith("ko"):
                continue
            info = kegg_db.get_pathway_info(kid)
            if not info["name"]:
                miss_term += 1
            # name 黑名单(人类/动物通路)|name blacklist (human/animal pathways)
            if _name_excluded(info["name"], exclude_keywords):
                excluded += 1
                continue
            # category 黑名单(可选, 如 Human Diseases 大类)|category blacklist
            if _category_excluded(info["category_a"], info["category_b"], exclude_categories):
                excluded += 1
                continue
            records.append((gene, kid, info["name"], info["category_a"]))
    return records, miss_term, excluded


def write_tsv(path: str, header: List[str], rows: List[tuple]):
    """写严格 4 列 TSV(无 BOM, \\n 换行)|Write strict 4-column TSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _parse_csv(opt: Optional[str]) -> List[str]:
    """逗号分隔串 → 列表(空/None → [])|Parse comma-separated string."""
    if not opt:
        return []
    return [x.strip() for x in opt.split(",") if x.strip()]


def build_tables(
    annotations_file: str,
    output_dir: str,
    sample: str,
    go_dict: Optional[Dict[str, Tuple[str, str]]] = None,
    kegg_db: Optional[KEGGDatabase] = None,
    kegg_exclude_keywords: Optional[str] = None,
    kegg_exclude_categories: Optional[str] = None,
    logger=None,
) -> dict:
    """
    主入口: 解析 eggnog annotations → 写 GO.tsv + KEGG.tsv|Main entry.

    Args:
        kegg_exclude_keywords: 通路 name 黑名单(逗号分隔). None=用内置 DEFAULT_PLANT_EXCLUDE
            (植物无关的癌症/激素/免疫等). ""=不过滤.
        kegg_exclude_categories: category A/B 黑名单(逗号分隔, 如 "Human Diseases"). 默认空.

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

    # KEGG name 黑名单: None → 内置植物无关词|None → built-in plant-irrelevant
    if kegg_exclude_keywords is None:
        kw_list = _parse_csv(DEFAULT_PLANT_EXCLUDE)
    else:
        kw_list = _parse_csv(kegg_exclude_keywords)
    cat_list = _parse_csv(kegg_exclude_categories)

    kegg_recs, kegg_miss, kegg_excluded = build_kegg_table(
        headers, rows, kegg_db,
        exclude_keywords=kw_list, exclude_categories=cat_list,
    )
    kegg_path = out / f"{sample}.kegg.tsv"
    write_tsv(str(kegg_path), KEGG_HEADER, kegg_recs)

    stats = {
        "go_rows": len(go_recs),
        "go_miss_term": go_miss,
        "kegg_rows": len(kegg_recs),
        "kegg_miss_term": kegg_miss,
        "kegg_excluded": kegg_excluded,
        "go_table": str(go_path),
        "kegg_table": str(kegg_path),
    }
    if logger:
        logger.info(
            f"GO 表|GO table: {stats['go_rows']} 行|rows "
            f"(term 缺失|missing: {stats['go_miss_term']}) → {go_path}"
        )
        logger.info(
            f"KEGG 表|KEGG table: {stats['kegg_rows']} 行|rows "
            f"(term 缺失|missing: {stats['kegg_miss_term']}, "
            f"按黑名单排除|excluded: {stats['kegg_excluded']}) → {kegg_path}"
        )
    return stats
