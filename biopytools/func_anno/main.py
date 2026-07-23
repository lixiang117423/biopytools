"""
func_anno 主程序|func_anno Main Entry.

端到端: interproscan(结构域) + eggnog-mapper(GO/KEGG 源) → 标准 GO/KEGG 表(衔接下游 R).
|End-to-end: interproscan (domains) + eggnog-mapper (GO/KEGG source) → standard
GO/KEGG tables for downstream R enrichment.

约束|Constraint: 不改 interproscan/eggnog_mapper 源码, 仅 import 调用(braker4ps 模式).

输入自动识别|Input auto-detection:
    -i 单文件 → 单样本(by-step, 不嵌套): output_dir/01_.../02_.../03_...
    -i 目录   → 多样本(by-sample, 每文件一子目录): output_dir/{sample}/01_.../...
    --by-sample 可强制单文件也嵌套(往同一 -o 多次跑不覆盖).
"""

import argparse
import logging
import os
import sys
from pathlib import Path

from .config import FuncAnnoConfig

# 蛋白序列扩展名|protein FASTA extensions
PROTEIN_EXTS = ("*.fa", "*.faa", "*.pep", "*.fasta")


def find_protein_files(input_dir: str) -> list:
    """扫描目录下蛋白序列文件|Find protein FASTA files in dir."""
    d = Path(input_dir)
    files = []
    for ext in PROTEIN_EXTS:
        files.extend(d.glob(ext))
    return sorted(set(files))


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
        epilog="示例|Examples:\n"
               "  单样本|single:  biopytools func-anno -i proteins.fa -o out/ -t 24\n"
               "  多样本|multi:   biopytools func-anno -i proteins_dir/ -o out/ -t 24",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="蛋白序列 FASTA(单文件→by-step) 或目录(多样本→by-sample)"
                        "|Protein FASTA (single→by-step) or dir (multi→by-sample)")
    parser.add_argument("-o", "--output-dir", required=True, help="输出目录|Output dir")
    parser.add_argument("-t", "--threads", type=int, default=12, help="线程数|Threads (default 12)")
    parser.add_argument("-s", "--sample-name", default=None,
                        help="样本名/输出前缀(默认输入文件名, 仅单文件模式生效)"
                        "|Sample name (single-file mode only)")
    parser.add_argument("--by-sample", action="store_true",
                        help="强制 by-sample(单文件也建 sample 子目录, 往同一 -o 多次跑不覆盖)"
                        "|Force by-sample layout")
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
    parser.add_argument("--kegg-exclude-keywords", default=None,
                        help="KEGG 通路 name 黑名单(逗号分隔子串, None=内置植物无关词 cancer/estrogen 等)"
                        "|KEGG name blacklist (None=built-in)")
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
    candidates = [d / f"{sample}.tsv", d / f"{sample}.proteins.tsv"]
    for c in candidates:
        if c.exists():
            return str(c)
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

    if cfg.ips_result:
        found = _find_ips_tsv(cfg.ips_result, cfg.sample_name, logger)
        if found:
            logger.info(f"复用 IPS 结果(跳过 IPS)|Reuse IPS (skip): {found}")
            return found
        logger.warning(f"未在 {cfg.ips_result} 找到 IPS TSV, 将重跑|"
                       f"No IPS TSV in {cfg.ips_result}, will rerun")

    if ips_tsv.exists():
        logger.info(f"IPS 已完成(断点续传)|IPS done (resume): {ips_tsv}")
        return str(ips_tsv)

    logger.info("-" * 70)
    logger.info("阶段1: interproscan 结构域注释|Phase 1: InterProScan domains")
    logger.info("-" * 70)
    cfg.ips_dir.mkdir(parents=True, exist_ok=True)

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
    # eggnog 实际输出在 output_dir/01_emapper/ 下(见 eggnog_mapper.utils.build_emapper_args)
    annotations = cfg.eggnog_dir / "01_emapper" / f"{cfg.sample_name}.emapper.annotations"

    if cfg.eggnog_result:
        if os.path.exists(cfg.eggnog_result):
            logger.info(f"复用 eggnog 结果(跳过 eggnog)|Reuse eggnog (skip): {cfg.eggnog_result}")
            return cfg.eggnog_result
        raise RuntimeError(f"指定的 eggnog 结果不存在|eggnog result not found: {cfg.eggnog_result}")

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
    # 建表是纯解析(秒级), 不做断点续传, 每次重跑重建(确保过滤参数生效).
    # |Tables are pure-parse (seconds); always rebuild to apply latest filter params.
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
                         go_dict=go_dict, kegg_db=kegg_db,
                         kegg_exclude_keywords=cfg.kegg_exclude_keywords,
                         kegg_exclude_categories=cfg.kegg_exclude_categories, logger=logger)
    logger.info("=" * 70)
    logger.info(f"GO 表|GO table: {stats['go_rows']} 行|rows "
                f"(term 缺失|missing: {stats['go_miss_term']})")
    logger.info(f"KEGG 表|KEGG table: {stats['kegg_rows']} 行|rows "
                f"(term 缺失|missing: {stats['kegg_miss_term']})")


def run_one_sample(input_file: str, args, by_sample: bool, sample_name=None):
    """跑单个样本三阶段|Run one sample through all 3 phases."""
    cfg = FuncAnnoConfig(
        input_file=input_file,
        output_dir=args.output_dir,
        threads=args.threads,
        sample_name=sample_name or args.sample_name,
        by_sample=by_sample,
        ips_result=args.ips_result,
        eggnog_result=args.eggnog_result,
        skip_ips=args.skip_ips,
        skip_eggnog=args.skip_eggnog,
        kegg_map=args.kegg_map,
        kegg_exclude_keywords=args.kegg_exclude_keywords,
        data_dir=args.data_dir,
        mode=args.mode,
        emapper_path=args.emapper_path,
    )
    cfg.validate()
    logger = _setup_logger(str(cfg.logs_dir))

    logger.info("=" * 70)
    logger.info(f"func_anno: IPS + eggnog → GO/KEGG 表|End-to-end | "
                f"sample={cfg.sample_name} | {'by-sample' if by_sample else 'by-step'}")
    logger.info("=" * 70)

    run_ips_phase(cfg, logger)
    annotations = run_eggnog_phase(cfg, logger)
    run_table_phase(cfg, annotations, logger)

    logger.info("=" * 70)
    logger.info(f"样本完成|Sample done: {cfg.sample_name} → {cfg.sample_dir}")
    logger.info("=" * 70)
    return cfg


def main():
    """主入口: 自动识别单/多样本|Main entry: auto-detect single/multi-sample."""
    args = parse_arguments()
    input_path = os.path.expanduser(args.input)

    try:
        if os.path.isdir(input_path):
            # 多样本: 目录 → by-sample|multi-sample: dir → by-sample
            protein_files = find_protein_files(input_path)
            if not protein_files:
                print(f"错误: 目录下无蛋白文件({', '.join(PROTEIN_EXTS)})|"
                      f"No protein files in: {input_path}", file=sys.stderr)
                sys.exit(1)
            print(f"多样本模式|Multi-sample (by-sample): {len(protein_files)} 个文件|files")
            failed = []
            for pf in protein_files:
                try:
                    run_one_sample(str(pf), args, by_sample=True)
                except Exception as e:
                    # 单样本失败不阻断其他|one failure does not block others
                    print(f"样本失败|Sample failed: {pf.name}: {e}", file=sys.stderr)
                    failed.append(pf.name)
            if failed:
                print(f"完成, 但 {len(failed)} 个样本失败|Done, but {len(failed)} failed: "
                      f"{failed}", file=sys.stderr)
            sys.exit(0)
        else:
            # 单样本: 单文件 → by-step(除非 --by-sample)|single → by-step
            run_one_sample(input_path, args, by_sample=args.by_sample)
            sys.exit(0)

    except SystemExit:
        raise
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
