"""
NLR-Annotator主程序模块|NLR-Annotator Main Module
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

from .config import NLRAnnotatorConfig
from .utils import (
    NLRLogger,
    build_conda_command,
    clean_output,
    collect_input_files,
    collect_result_files,
    filter_contained_calls,
    generate_summary,
    generate_summary_gff,
    gff_path_for,
    removed_tsv_path,
    tsv_to_gff,
)


def build_command(config: NLRAnnotatorConfig, input_file: str, output_file: str) -> list:
    """构建java命令|Build java command list"""
    # 用build_conda_command包装java(支持conda env的java,激活原死代码)|Wrap java via build_conda_command
    args = ['-jar', config.jar_path]
    args += ['-i', input_file]
    args += ['-x', config.mot_file]
    args += ['-y', config.store_file]
    args += ['-t', str(config.threads)]
    args += ['-n', str(config.num_seqs_per_thread)]
    args += ['-o', output_file]

    if config.output_bed:
        bed_path = output_file.rsplit('.', 1)[0] + '.bed'
        args += ['-b', bed_path]
    if config.output_motifs:
        motifs_path = output_file.rsplit('.', 1)[0] + '_motifs.bed'
        args += ['-m', motifs_path]
    if config.output_alignment:
        align_path = output_file.rsplit('.', 1)[0] + '_alignment.fa'
        args += ['-a', align_path]

    args += ['-distanceWithinMotifCombination', str(config.distance_within_motif_combination)]
    args += ['-distanceForElongating', str(config.distance_for_elongating)]
    args += ['-distanceBetweenMotifCombinations', str(config.distance_between_motif_combinations)]

    return build_conda_command(config.java_path, args)


def _is_step_completed(output_file: str) -> bool:
    """检查步骤是否已完成|Check if step is completed"""
    return Path(output_file).exists() and os.path.getsize(output_file) > 0


def _run_single(config: NLRAnnotatorConfig, input_file: str, sample_name: str,
                 sample_output_dir: Path, logger: logging.Logger) -> str:
    """
    运行单个样本的NLR-Annotator|Run NLR-Annotator for a single sample

    Args:
        config: 配置对象|Config object
        input_file: 输入文件路径|Input file path
        sample_name: 样本名|Sample name
        sample_output_dir: 样本输出目录|Sample output directory
        logger: 日志器|Logger

    Returns:
        输出TSV文件路径|Output TSV file path
    """
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    output_file = str(sample_output_dir / f"{sample_name}_nlr_annotator.tsv")

    # 断点续传：输出已存在则跳过java(冗余包含过滤在下方对新旧结果统一应用)
    # |Checkpoint resume: skip java if output exists (containment filtering below applies to both fresh & existing results)
    if _is_step_completed(output_file):
        logger.info(f"跳过已完成样本|Skipping completed sample: {sample_name}")
    else:
        # 结果重建前清理陈旧留档,避免与新结果不一致|Purge stale audit file before regenerating results
        stale_removed = removed_tsv_path(output_file)
        if os.path.exists(stale_removed):
            os.remove(stale_removed)
            logger.info(f"清理陈旧留档|Purged stale audit file: {stale_removed}")

        cmd = build_command(config, input_file, output_file)
        logger.info(f"处理样本|Processing sample: {sample_name}")
        logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"NLR-Annotator运行失败(退出码{result.returncode})|NLR-Annotator failed (exit {result.returncode})")
                if result.stderr:
                    logger.error(f"错误信息|Error message: {result.stderr}")
                sys.exit(1)
        except FileNotFoundError as e:
            logger.error(f"命令执行失败|Command execution failed: {e}")
            sys.exit(1)

        logger.info(f"样本完成|Sample done: {sample_name}")

        # 清洗输出文件：加表头、去重motif|Clean output: add header, deduplicate motifs
        clean_output(output_file, logger)

    # 冗余包含过滤:幂等,断点旧结果重跑同命令即可原地过滤(无需重跑java)
    # |Containment filtering: idempotent; rerun the same command to filter checkpointed results in place (no java rerun)
    if config.filter_contained:
        filter_contained_calls(output_file, logger)

    # GFF3以过滤后TSV为唯一数据源,断点跳过java时也重新生成(幂等覆盖,旧结果原地补GFF)
    # |GFF3 uses the filtered TSV as its single source; regenerated even on checkpoint skip (idempotent overwrite, backfills old results)
    tsv_to_gff(output_file, gff_path_for(output_file), logger)

    return output_file


def run_merge_only(input_path: str, output_path, logger: logging.Logger):
    """
    纯合并模式:合并已有各基因组NLR结果为总表(不执行NLR-Annotator)
    |Merge-only: merge existing per-genome NLR results into a summary (skip java)

    用于已有批处理结果但未生成汇总文件的情况(如运行中途被杀)。
    |Used when batch results exist but no summary was produced (e.g. job killed mid-run).

    Args:
        input_path: 结果文件/目录(支持by-sample子目录或平铺)|Result file/dir
        output_path: 输出目录(写 nlr_annotator_summary.tsv)|Output directory
        logger: 日志器|Logger
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    sample_results = collect_result_files(input_path, logger)
    logger.info(f"开始合并|Start merging: {len(sample_results)} 个基因组|genome(s)")
    generate_summary(sample_results, output_path, logger)
    generate_summary_gff(sample_results, output_path, logger)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='NLR-Annotator: 从DNA/CDS序列预测NLR基因|Predict NLR genes from DNA/CDS sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入DNA/CDS FASTA文件或目录|Input DNA/CDS FASTA file or directory')
    parser.add_argument('-o', '--output-dir', default='./output',
                        help='输出目录|Output directory (default: ./output)')
    parser.add_argument('--sample-suffix', default='*.fa',
                        help='目录模式下文件匹配后缀|File match suffix for directory mode (default: *.fa)')

    parser.add_argument('--merge-only', action='store_true',
                        help='只合并已有结果TSV(*_nlr_annotator.tsv),不运行NLR-Annotator'
                             '|Merge existing result TSVs only, skip NLR-Annotator')

    parser.add_argument('--no-filter-contained', action='store_true',
                        help='关闭被包含冗余调用过滤(默认开启:剔除同序列上被完整基因完全包含的短片段调用,'
                             '被剔除记录留档为*_nlr_annotator_removed.tsv)'
                             '|Disable contained-call filtering (default ON: drop calls fully contained in another call on the same sequence, archived to *_removed.tsv)')

    parser.add_argument('--jar-path', default='',
                        help='NLR-Annotator JAR文件路径|NLR-Annotator JAR file path')
    parser.add_argument('--mot-file', default='',
                        help='mot.txt配置文件路径|mot.txt config file path')
    parser.add_argument('--store-file', default='',
                        help='store.txt配置文件路径|store.txt config file path')
    parser.add_argument('--java-path', default='java',
                        help='Java解释器路径(默认系统java;conda env用~/miniforge3/envs/xxx/bin/java)|Java interpreter path')

    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('--num-seqs-per-thread', type=int, default=1000,
                        help='每线程处理序列数|Sequences per thread (default: 1000)')

    parser.add_argument('--output-bed', action='store_true',
                        help='输出BED文件|Output BED file')
    parser.add_argument('--output-motifs', action='store_true',
                        help='输出motifs BED文件|Output motifs BED file')
    parser.add_argument('--output-alignment', action='store_true',
                        help='输出motif比对FASTA|Output motif alignment FASTA')

    parser.add_argument('--distance-within-motif-combination', type=int, default=500,
                        help='motif组合内距离|Distance within motif combination (default: 500)')
    parser.add_argument('--distance-for-elongating', type=int, default=2500,
                        help='延伸距离|Distance for elongating (default: 2500)')
    parser.add_argument('--distance-between-motif-combinations', type=int, default=50000,
                        help='motif组合间距离|Distance between motif combinations (default: 50000)')

    args = parser.parse_args()

    # 先创建输出目录和日志管理器（用于记录验证错误）|Create output dir and logger first
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_manager = NLRLogger(output_dir)
    logger = log_manager.get_logger()

    # 纯合并模式:跳过NLR-Annotator执行,直接合并已有结果(不校验JAR)
    # |Merge-only: skip NLR-Annotator, merge existing results directly (no JAR check)
    if args.merge_only:
        logger.info("纯合并模式|Merge-only mode: 跳过NLR-Annotator执行|skip NLR-Annotator execution")
        try:
            run_merge_only(args.input, output_dir, logger)
        except ValueError as e:
            logger.error(f"合并失败|Merge failed: {e}")
            sys.exit(1)
        logger.info("合并完成|Merge done")
        return

    config = NLRAnnotatorConfig(
        input_path=args.input,
        output_dir=output_dir,
        sample_suffix=args.sample_suffix,
        jar_path=args.jar_path,
        mot_file=args.mot_file,
        store_file=args.store_file,
        java_path=args.java_path,
        threads=args.threads,
        num_seqs_per_thread=args.num_seqs_per_thread,
        output_bed="1" if args.output_bed else "",
        output_motifs="1" if args.output_motifs else "",
        output_alignment="1" if args.output_alignment else "",
        distance_within_motif_combination=args.distance_within_motif_combination,
        distance_for_elongating=args.distance_for_elongating,
        distance_between_motif_combinations=args.distance_between_motif_combinations,
        filter_contained=not args.no_filter_contained,
    )

    try:
        config.validate()
    except ValueError as e:
        logger.error(f"配置错误|Configuration error:\n{e}")
        sys.exit(1)

    # 收集输入文件|Collect input files
    try:
        input_files = collect_input_files(config.input_path, config.sample_suffix, logger)
    except ValueError as e:
        logger.error(f"收集输入文件失败|Failed to collect input files: {e}")
        sys.exit(1)
    is_batch = len(input_files) > 1

    logger.info(f"共{len(input_files)}个样本待处理|Total {len(input_files)} sample(s) to process")

    # 逐样本运行|Run per sample
    sample_results = []  # [(sample_name, output_tsv_path)]
    for input_file, sample_name in input_files:
        if is_batch:
            sample_output_dir = config.output_path / sample_name
            sample_log_dir = sample_output_dir / "99_logs"
            sample_log_dir.mkdir(parents=True, exist_ok=True)
            sample_log_manager = NLRLogger(str(sample_output_dir), f"{sample_name}_nlr_annotator.log")
            sample_logger = sample_log_manager.get_logger()
        else:
            sample_output_dir = config.output_path
            sample_logger = logger

        output_tsv = _run_single(config, input_file, sample_name, sample_output_dir, sample_logger)
        sample_results.append((sample_name, output_tsv))

    # 目录模式：生成汇总文件|Directory mode: generate summary file
    if is_batch:
        generate_summary(sample_results, config.output_path, logger)
        generate_summary_gff(sample_results, config.output_path, logger)

    logger.info(f"全部完成|All done: {len(input_files)} 个样本|sample(s)")


if __name__ == "__main__":
    main()
