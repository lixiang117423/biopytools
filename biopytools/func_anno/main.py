"""
func_anno 主程序|func_anno Main Entry.

端到端: interproscan(结构域) + eggnog-mapper(GO/KEGG 源) → 标准 GO/KEGG 表(衔接下游 R).
|End-to-end: interproscan (domains) + eggnog-mapper (GO/KEGG source) → standard
GO/KEGG tables for downstream R enrichment.

约束|Constraint: 不改 interproscan/eggnog_mapper 源码, 仅 import 调用(braker4ps 模式).
|import-only (braker4ps pattern).

阶段|Phases:
    1. IPS(可选): 结构域注释, 不影响 GO/KEGG 表. 支持 --ips-result 复用已跑结果.
    2. eggnog(必需): GO/KEGG 唯一数据源. 支持 --eggnog-result 复用.
    3. 建表: 解析 eggnog annotations → 标准 GO.tsv + KEGG.tsv.
"""

import argparse
import logging
import os
import sys
from pathlib import Path

from .config import FuncAnnoConfig


def _setup_logger(logs_dir: str) -> logging.Logger:
    """统一日志(stdout INFO + stderr WARNING + file, 遵循 §2.3)|Unified logger."""
    logger = logging.getLogger("func_anno")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter(
        "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setLevel(logging.INFO)
    stdout_h.setFormatter(fmt)
    logger.addHandler(stdout_h)

    stderr_h = logging.StreamHandler(sys.stderr)
    stderr_h.setLevel(logging.WARNING)
    stderr_h.setFormatter(fmt)
    logger.addHandler(stderr_h)

    os.makedirs(logs_dir, exist_ok=True)
    file_h = logging.FileHandler(os.path.join(logs_dir, "func_anno.log"), encoding="utf-8")
    file_h.setLevel(logging.DEBUG)
    file_h.setFormatter(fmt)
    logger.addHandler(file_h)
    return logger


def parse_arguments():
    """解析命令行参数|Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="func_anno: 蛋白功能注释(IPS+eggnog→GO/KEGG标准表)"
        "|Protein functional annotation (IPS+eggnog→GO/KEGG tables)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Example: biopytools func-anno -i proteins.fa -o out/ -t 24",
    )
    parser.add_argument("-i", "--input", required=True, help="蛋白序列 FASTA|Protein FASTA")
    parser.add_argument("-o", "--output-dir", required=True, help="输出目录|Output dir")
    parser.add_argument("-t", "--threads", type=int, default=12, help="线程数|Threads (default 12)")
    parser.add_argument("-s", "--sample-name", default=None,
                        help="样本名/输出前缀(默认输入文件名)|Sample name (default: input stem)")
    # 复用|reuse
    parser.add_argument("--ips-result", default=None,
                        help="复用已有 IPS 结果目录(跳过 IPS)|Reuse existing IPS dir")
    parser.add_argument("--eggnog-result", default=None,
                        help="复用已有 .emapper.annotations(跳过 eggnog)|Reuse existing annotations")
    parser.add_argument("--skip-ips", action="store_true",
                        help="跳过 IPS(只要 GO/KEGG)|Skip IPS")
    parser.add_argument("--skip-eggnog", action="store_true",
                        help="跳过 eggnog(仅整理已有结果)|Skip eggnog")
    # KEGG|KEGG
    parser.add_argument("--kegg-map", default=None,
                        help="外部 KEGG 映射 TSV(ko_id\\tname\\tcategory, 补 category)"
                        "|External KEGG map to fill category")
    # eggnog 透传|eggnog passthrough
    parser.add_argument("-m", "--mode", default="mmseqs",
                        choices=["mmseqs", "diamond", "hmmer"], help="搜索模式|Search mode")
    parser.add_argument("--data-dir", default=None, help="eggnog DB 目录|eggnog DB dir")
    parser.add_argument("--emapper-path", default=None,
                        help="emapper.py 路径覆盖|emapper.py path override")
    return parser.parse_args()


def _find_ips_tsv(ips_dir: str, sample: str, logger) -> str:
    """在已有 IPS 目录找 TSV(优先 {sample}.tsv, 否则任意非 cn.tsv)|Find IPS TSV."""
    d = Path(ips_dir)
    # 优先样本名匹配|prefer sample-name match
    candidates = [d / f"{sample}.tsv", d / f"{sample}.proteins.tsv"]
    for c in candidates:
        if c.exists():
            return str(c)
    # 任意 TSV(排除中文重排版 cn.tsv)|any TSV (exclude cn.tsv)
    tsvs = sorted(d.rglob("*.tsv"))
    tsvs = [t for t in tsvs if not t.name.endswith(".cn.tsv")]
    if tsvs:
        return str(tsvs[0])
    return ""


def run_ips_phase(cfg: FuncAnnoConfig, logger: logging.Logger) -> str:
    """阶段1: interproscan(结构域, 可选, 不影响 GO/KEGG 表)|Phase 1 IPS (optional)."""
    if cfg.skip_ips:
        logger.info("跳过 IPS|Skip IPS (domains not needed)")
        return ""

    ips_tsv = cfg.ips_dir / f"{cfg.sample_name}.tsv"

    # 复用已有结果|reuse existing
    if cfg.ips_result:
        found = _find_ips_tsv(cfg.ips_result, cfg.sample_name, logger)
        if found:
            logger.info(f"复用 IPS 结果(跳过 IPS)|Reuse IPS (skip): {found}")
            return found
        logger.warning(f"未在 {cfg.ips_result} 找到 IPS TSV, 将重跑|"
                       f"No IPS TSV in {cfg.ips_result}, will rerun")

    # 断点续传|resume
    if ips_tsv.exists():
        logger.info(f"IPS 已完成(断点续传)|IPS done (resume): {ips_tsv}")
        return str(ips_tsv)

    logger.info("-" * 70)
    logger.info("阶段1: interproscan 结构域注释|Phase 1: InterProScan domains")
    logger.info("-" * 70)
    cfg.ips_dir.mkdir(parents=True, exist_ok=True)

    # 延迟 import(避免 help 时加载重依赖)|lazy import
    from ..interproscan import InterProScanAnnotator

    annotator = InterProScanAnnotator(
        input_file=cfg.input_file,
        output_prefix=str(cfg.ips_dir / cfg.sample_name),
        output_format="tsv,xml",
        threads=cfg.threads,
    )
    ok = annotator.run_analysis()
    if ok and ips_tsv.exists():
        logger.info(f"IPS 完成|IPS done: {ips_tsv}")
        return str(ips_tsv)
    logger.warning("IPS 未成功或无 TSV(不影响 GO/KEGG 表)|"
                   "IPS failed/no TSV (does not affect GO/KEGG tables)")
    return ""


def run_eggnog_phase(cfg: FuncAnnoConfig, logger: logging.Logger) -> str:
    """阶段2: eggnog-mapper(必需, GO/KEGG 源)|Phase 2 eggnog (required)."""
    annotations = cfg.eggnog_dir / f"{cfg.sample_name}.emapper.annotations"

    # 复用已有结果|reuse existing
    if cfg.eggnog_result:
        if os.path.exists(cfg.eggnog_result):
            logger.info(f"复用 eggnog 结果(跳过 eggnog)|Reuse eggnog (skip): {cfg.eggnog_result}")
            return cfg.eggnog_result
        raise RuntimeError(f"指定的 eggnog 结果不存在|eggnog result not found: {cfg.eggnog_result}")

    # 断点续传|resume
    if annotations.exists():
        logger.info(f"eggnog 已完成(断点续传)|eggnog done (resume): {annotations}")
        return str(annotations)

    logger.info("-" * 70)
    logger.info("阶段2: eggnog-mapper 功能注释(GO/KEGG 源)|Phase 2: eggNOG-mapper")
    logger.info("-" * 70)
    cfg.eggnog_dir.mkdir(parents=True, exist_ok=True)

    from ..eggnog_mapper.config import EggnogMapperConfig
    from ..eggnog_mapper.utils import EggnogMapperRunner

    ekwargs = dict(
        input_file=cfg.input_file,
        output_dir=str(cfg.eggnog_dir),
        itype="proteins",
        mode=cfg.mode,
        cpu=cfg.threads,
        prefix=cfg.sample_name,
    )
    if cfg.data_dir:
        ekwargs["data_dir"] = cfg.data_dir
    if cfg.emapper_path:
        ekwargs["emapper_path"] = cfg.emapper_path

    ecfg = EggnogMapperConfig(**ekwargs)
    ecfg.validate()
    ok = EggnogMapperRunner(ecfg, logger).run()
    if not ok or not annotations.exists():
        raise RuntimeError("eggnog-mapper 失败或无 annotations|eggnog-mapper failed/no annotations")
    logger.info(f"eggnog 完成|eggnog done: {annotations}")
    return str(annotations)


def run_table_phase(cfg: FuncAnnoConfig, annotations: str, logger: logging.Logger):
    """阶段3: 建 GO/KEGG 标准表|Phase 3: build standard tables."""
    go_tsv = cfg.tables_dir / f"{cfg.sample_name}.go.tsv"
    kegg_tsv = cfg.tables_dir / f"{cfg.sample_name}.kegg.tsv"

    # 断点续传|resume
    if go_tsv.exists() and kegg_tsv.exists():
        logger.info(f"表已完成(断点续传)|Tables done (resume): {go_tsv}, {kegg_tsv}")
        return

    logger.info("-" * 70)
    logger.info("阶段3: 建 GO/KEGG 标准表(衔接 R)|Phase 3: Build GO/KEGG tables")
    logger.info("-" * 70)
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)

    from .table_builder import build_tables, load_go_dict
    from .kegg_db import KEGGDatabase

    logger.info("加载 GO 映射(内置 go_data)|Loading GO map (built-in go_data)...")
    go_dict = load_go_dict()
    logger.info(f"GO 映射|GO map: {len(go_dict)} 条|entries")
    kegg_db = KEGGDatabase(kegg_map_file=cfg.kegg_map, logger=logger)

    stats = build_tables(annotations, str(cfg.tables_dir), cfg.sample_name,
                         go_dict=go_dict, kegg_db=kegg_db, logger=logger)
    logger.info("=" * 70)
    logger.info(f"GO 表|GO table: {stats['go_rows']} 行|rows "
                f"(term 缺失|missing: {stats['go_miss_term']})")
    logger.info(f"KEGG 表|KEGG table: {stats['kegg_rows']} 行|rows "
                f"(term 缺失|missing: {stats['kegg_miss_term']})")


def main():
    """主入口: IPS → eggnog → 建表|Main entry."""
    args = parse_arguments()
    cfg = FuncAnnoConfig(
        input_file=args.input,
        output_dir=args.output_dir,
        threads=args.threads,
        sample_name=args.sample_name,
        ips_result=args.ips_result,
        eggnog_result=args.eggnog_result,
        skip_ips=args.skip_ips,
        skip_eggnog=args.skip_eggnog,
        kegg_map=args.kegg_map,
        data_dir=args.data_dir,
        mode=args.mode,
        emapper_path=args.emapper_path,
    )
    cfg.validate()

    logger = _setup_logger(str(cfg.logs_dir))

    try:
        logger.info("=" * 70)
        logger.info(f"func_anno: IPS + eggnog → GO/KEGG 表|End-to-end | "
                    f"sample={cfg.sample_name}")
        logger.info("=" * 70)

        # 阶段1 IPS(可选)|Phase 1 IPS (optional)
        run_ips_phase(cfg, logger)

        # 阶段2 eggnog(必需)|Phase 2 eggnog (required)
        annotations = run_eggnog_phase(cfg, logger)

        # 阶段3 建表|Phase 3 tables
        run_table_phase(cfg, annotations, logger)

        logger.info("=" * 70)
        logger.info("func_anno 完成|func_anno done")
        logger.info(f"输出|Output: {cfg.sample_dir}")
        logger.info("=" * 70)
        sys.exit(0)

    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
