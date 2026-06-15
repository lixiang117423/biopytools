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
from .utils import NLRLogger, clean_output, collect_input_files, generate_summary


def build_command(config: NLRAnnotatorConfig, input_file: str, output_file: str) -> list:
    """构建java命令|Build java command list"""
    cmd = ['java', '-jar', config.jar_path]
    cmd += ['-i', input_file]
    cmd += ['-x', config.mot_file]
    cmd += ['-y', config.store_file]
    cmd += ['-t', str(config.threads)]
    cmd += ['-n', str(config.num_seqs_per_thread)]
    cmd += ['-o', output_file]

    if config.output_gff:
        gff_path = output_file.rsplit('.', 1)[0] + '.gff'
        cmd += ['-g', gff_path]
    if config.output_bed:
        bed_path = output_file.rsplit('.', 1)[0] + '.bed'
        cmd += ['-b', bed_path]
    if config.output_motifs:
        motifs_path = output_file.rsplit('.', 1)[0] + '_motifs.bed'
        cmd += ['-m', motifs_path]
    if config.output_alignment:
        align_path = output_file.rsplit('.', 1)[0] + '_alignment.fa'
        cmd += ['-a', align_path]

    cmd += ['-distanceWithinMotifCombination', str(config.distance_within_motif_combination)]
    cmd += ['-distanceForElongating', str(config.distance_for_elongating)]
    cmd += ['-distanceBetweenMotifCombinations', str(config.distance_between_motif_combinations)]

    return cmd


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
    output_file = str(sample_output_dir / f"{sample_name}.nlr_annotator.tsv")

    # 断点续传：检查输出文件是否已存在|Checkpoint resume: check if output already exists
    if _is_step_completed(output_file):
        logger.info(f"跳过已完成样本|Skipping completed sample: {sample_name}")
        return output_file

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

    return output_file


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='NLR-Annotator: 从CDS序列预测NLR基因|Predict NLR genes from CDS sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入CDS FASTA文件或目录|Input CDS FASTA file or directory')
    parser.add_argument('-o', '--output-dir', default='./output',
                        help='输出目录|Output directory (default: ./output)')
    parser.add_argument('--sample-suffix', default='*.cds.fa',
                        help='目录模式下文件匹配后缀|File match suffix for directory mode (default: *.cds.fa)')

    parser.add_argument('--jar-path', default='',
                        help='NLR-Annotator JAR文件路径|NLR-Annotator JAR file path')
    parser.add_argument('--mot-file', default='',
                        help='mot.txt配置文件路径|mot.txt config file path')
    parser.add_argument('--store-file', default='',
                        help='store.txt配置文件路径|store.txt config file path')

    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('--num-seqs-per-thread', type=int, default=1000,
                        help='每线程处理序列数|Sequences per thread (default: 1000)')

    parser.add_argument('--output-gff', action='store_true',
                        help='输出GFF文件|Output GFF file')
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

    config = NLRAnnotatorConfig(
        input_path=args.input,
        output_dir=output_dir,
        sample_suffix=args.sample_suffix,
        jar_path=args.jar_path,
        mot_file=args.mot_file,
        store_file=args.store_file,
        threads=args.threads,
        num_seqs_per_thread=args.num_seqs_per_thread,
        output_gff="1" if args.output_gff else "",
        output_bed="1" if args.output_bed else "",
        output_motifs="1" if args.output_motifs else "",
        output_alignment="1" if args.output_alignment else "",
        distance_within_motif_combination=args.distance_within_motif_combination,
        distance_for_elongating=args.distance_for_elongating,
        distance_between_motif_combinations=args.distance_between_motif_combinations,
    )

    try:
        config.validate()
    except ValueError as e:
        logger.error(f"配置错误|Configuration error:\n{e}")
        sys.exit(1)

    # 收集输入文件|Collect input files
    input_files = collect_input_files(config.input_path, config.sample_suffix, logger)
    is_batch = len(input_files) > 1

    logger.info(f"共{len(input_files)}个样本待处理|Total {len(input_files)} sample(s) to process")

    # 逐样本运行|Run per sample
    sample_results = []  # [(sample_name, output_tsv_path)]
    for input_file, sample_name in input_files:
        if is_batch:
            sample_output_dir = config.output_path / sample_name
            sample_log_dir = sample_output_dir / "99_logs"
            sample_log_dir.mkdir(parents=True, exist_ok=True)
            sample_log_manager = NLRLogger(str(sample_output_dir), f"{sample_name}.nlr_annotator.log")
            sample_logger = sample_log_manager.get_logger()
        else:
            sample_output_dir = config.output_path
            sample_logger = logger

        output_tsv = _run_single(config, input_file, sample_name, sample_output_dir, sample_logger)
        sample_results.append((sample_name, output_tsv))

    # 目录模式：生成汇总文件|Directory mode: generate summary file
    if is_batch:
        generate_summary(sample_results, config.output_path, logger)

    logger.info(f"全部完成|All done: {len(input_files)} 个样本|sample(s)")


if __name__ == "__main__":
    main()
