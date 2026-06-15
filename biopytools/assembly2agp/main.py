#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Assembly to AGP Converter|Assembly转AGP格式转换器

将assembly文件转换为AGP格式，并生成染色体列表文件|Convert assembly files to AGP format and generate chromosome list
Author: BioTools Development Team
Version: 1.0.0
"""

import argparse
import logging
import sys
import os
from pathlib import Path
import pandas as pd
from collections import OrderedDict


# ==================== 日志配置|Logging Configuration ====================

def setup_logger(name, log_file=None, level=logging.INFO):
    """
    配置标准化的日志系统|Configure standardized logging system

    Args:
        name: logger名称|Logger name
        log_file: 日志文件路径(可选)|Log file path (optional)
        level: 日志级别|Log level

    Returns:
        logger对象|Logger object
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # 日志格式|Log format
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # stdout handler - INFO级别|stdout handler - INFO level
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # 文件handler(如果指定)|File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


# ==================== AGP转换核心功能|AGP Conversion Core Functions ====================

AGP_HEADER = ["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID",
              "Contig_start", "Contig_end", "Orientation"]


def seq2agp(seq, scaffoldname, idx2scffold, gap=100):
    """
    将序列转换为AGP格式|Convert sequence to AGP format

    Args:
        seq: 序列列表|Sequence list
        scaffoldname: scaffold名称|Scaffold name
        idx2scffold: 索引到scaffold的映射字典|Index to scaffold mapping dictionary
        gap: gap大小|Gap size

    Returns:
        DataFrame: AGP格式的DataFrame|DataFrame in AGP format
    """
    all_agp = []
    for i in seq:
        if i < 0:
            item = idx2scffold[abs(i)]
            oritention = "-"
        else:
            oritention = "+"
            item = idx2scffold[abs(i)]
        Chromosome = scaffoldname
        Start = 1
        End = int(item[2])
        Order = 1
        Tag = "W"
        Contig_ID = item[0]
        Contig_start = Start
        Contig_end = End
        Orientation = oritention
        temp_data = [Chromosome, Start, End, Order, Tag, Contig_ID,
                     Contig_start, Contig_end, Orientation]
        gap_item = [Chromosome, Start, End, Order, "U", gap, "scaffold",
                    "yes", "proximity_ligation"]
        all_agp.append(temp_data)
        all_agp.append(gap_item)
    agp = pd.DataFrame(data=all_agp[:-1], columns=AGP_HEADER)
    agp.iloc[0, 3] = 1
    agp.iloc[0, 1] = int(agp.iloc[0, 6])
    agp.iloc[0, 2] = int(agp.iloc[0, 7])
    for i in range(1, len(agp)):
        if agp.iloc[i, 4] == "W":
            agp.iloc[i, 3] = agp.iloc[i - 1, 3] + 1
            agp.iloc[i, 1] = int(agp.iloc[i - 1, 2])
            agp.iloc[i, 2] = int(agp.iloc[i - 1, 2]) + int(agp.iloc[i, 7])
        else:
            agp.iloc[i, 3] = agp.iloc[i - 1, 3] + 1
            agp.iloc[i, 1] = int(agp.iloc[i - 1, 2])
            agp.iloc[i, 2] = int(agp.iloc[i - 1, 2]) + int(agp.iloc[i, 5])
    return agp


def assembly2agp(assembly, agp_output, gap=100, num_chromosomes=None, logger=None):
    """
    将assembly文件转换为AGP格式|Convert assembly file to AGP format

    Args:
        assembly: assembly文件路径|Assembly file path
        agp_output: 输出AGP文件路径|Output AGP file path
        gap: gap大小|Gap size
        num_chromosomes: 染色体数量限制（None表示全部输出）|Chromosome number limit (None for all)
        logger: 日志对象|Logger object

    Returns:
        tuple: (agp文件路径, scaffold信息列表)|(agp file path, scaffold info list)
    """
    idx2scffold = {}
    count = 0
    agp_list = []
    all_agp = pd.DataFrame(data=[], columns=AGP_HEADER)

    if logger:
        logger.info(f"读取assembly文件|Reading assembly file: {assembly}")

    with open(assembly, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                items = line[1:].split(" ")
                idx2scffold[int(items[1])] = items
            else:
                items = line.split(" ")
                seq = [int(x) for x in items]
                count += 1
                tmp_agp = seq2agp(seq, f"scaffold_{count}", idx2scffold, gap=gap)
                agp_list.append([tmp_agp.iloc[-1, 2], tmp_agp])

    # 按长度排序（从大到小）|Sort by length (descending)
    agp_list.sort(key=lambda i: i[0], reverse=True)

    # 如果指定了染色体数量，只保留前N个scaffolds|If chromosome number specified, keep only top N scaffolds
    if num_chromosomes is not None:
        if logger:
            logger.info(f"限制输出为前 {num_chromosomes} 个scaffolds|Limiting output to top {num_chromosomes} scaffolds")
        agp_list = agp_list[:num_chromosomes]

    # 重命名scaffold|Rename scaffolds
    for i in range(len(agp_list)):
        agp_list[i][1].iloc[:, 0] = f"scaffold_{i + 1}"
        all_agp = pd.concat([all_agp, agp_list[i][1]])

    # 保存AGP文件|Save AGP file
    all_agp.to_csv(agp_output, sep="\t", header=False, index=False)

    if logger:
        logger.info(f"生成AGP文件|Generated AGP file: {agp_output}")
        logger.info(f"输出scaffold总数|Total scaffolds in output: {len(agp_list)}")

    return agp_output, agp_list


def generate_chr_list(agp_file, chr_list_output, num_chromosomes, logger=None):
    """
    从AGP文件生成染色体列表文件|Generate chromosome list file from AGP

    Args:
        agp_file: AGP文件路径|AGP file path
        chr_list_output: 输出chr.list文件路径|Output chr.list file path
        num_chromosomes: 染色体数量|Number of chromosomes
        logger: 日志对象|Logger object

    Returns:
        chr_list文件路径|chr.list file path
    """
    if logger:
        logger.info(f"从AGP生成chr.list文件|Generating chr.list file from AGP: {agp_file}")
        logger.info(f"染色体数量|Number of chromosomes: {num_chromosomes}")

    # 读取AGP文件|Read AGP file
    agp_df = pd.read_csv(agp_file, sep="\t", header=None, names=AGP_HEADER)

    # 获取每个scaffold的长度|Get length of each scaffold
    scaffold_lengths = {}
    for chromosome in agp_df['Chromosome'].unique():
        chrom_df = agp_df[agp_df['Chromosome'] == chromosome]
        max_end = chrom_df['End'].max()
        scaffold_lengths[chromosome] = max_end

    # 按长度排序（从大到小）|Sort by length (descending)
    sorted_scaffolds = sorted(scaffold_lengths.items(),
                             key=lambda x: x[1], reverse=True)

    # 选择前N个scaffolds作为染色体|Select top N scaffolds as chromosomes
    selected_chromosomes = sorted_scaffolds[:num_chromosomes]

    # 写入chr.list文件|Write to chr.list file
    with open(chr_list_output, 'w') as f:
        for chrom_name, chrom_length in selected_chromosomes:
            f.write(f"{chrom_name}\t{chrom_length}\n")

    if logger:
        logger.info(f"生成chr.list文件|Generated chr.list file: {chr_list_output}")
        logger.info(f"染色体1长度|Chromosome 1 length: {selected_chromosomes[0][1]:,} bp")
        logger.info(f"染色体{num_chromosomes}长度|Chromosome {num_chromosomes} length: {selected_chromosomes[-1][1]:,} bp")

        # 输出所有染色体的长度统计|Output length statistics for all chromosomes
        total_length = sum(length for _, length in selected_chromosomes)
        logger.info(f"前 {num_chromosomes} 个染色体总长度|Total length of top {num_chromosomes} chromosomes: {total_length:,} bp")

    return chr_list_output


# ==================== 参数解析|Argument Parsing ====================

def parse_arguments():
    """
    解析命令行参数|Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Assembly to AGP Converter - Convert assembly files to AGP format and generate chromosome list|Assembly转AGP格式转换器 - 将assembly文件转换为AGP格式并生成染色体列表',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -a corrected_asm.FINAL.assembly -p output_prefix -n 12
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', '--assembly', required=True,
                         help='Input assembly file path')
    required.add_argument('-p', '--prefix', required=True,
                         help='Output prefix (will be used for both AGP and chr.list files)')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-o', '--output-dir', default='.',
                         help='Output directory path (default: current directory)')
    optional.add_argument('-g', '--gap', type=int, default=100,
                         help='Gap size between scaffolds in bp (default: 100)')
    optional.add_argument('-n', '--num-chromosomes', type=int, required=True,
                         help='Number of chromosomes to include in chr.list file')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('logging options')
    log_group.add_argument('-v', '--verbose', action='count', default=0,
                          help='Verbose mode (-v: INFO, -vv: DEBUG)')
    log_group.add_argument('--quiet', action='store_true',
                          help='Quiet mode (only ERROR)')
    log_group.add_argument('--log-file',
                          help='Log file path')

    # 执行控制|Execution control
    exec_group = parser.add_argument_group('execution options')
    exec_group.add_argument('-f', '--force', action='store_true',
                           help='Force overwrite existing files')

    # 版本信息|Version information
    parser.add_argument('-V', '--version', action='version',
                       version='%(prog)s 1.0.0')

    return parser.parse_args()


# ==================== 主函数 ====================

def main():
    """
    主函数|Main function
    """
    import time
    start_time = time.time()

    # 解析参数|Parse arguments
    args = parse_arguments()

    # 设置日志级别|Set log level
    if args.quiet:
        log_level = logging.ERROR
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    elif args.verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    # 初始化日志|Initialize logger
    logger = setup_logger(__name__, args.log_file, log_level)

    # 输出程序信息|Output program information
    logger.info("=" * 60)
    logger.info("Program: Assembly to AGP Converter")
    logger.info("Version: 1.0.0")
    logger.info("=" * 60)
    logger.info(f"Input assembly: {args.assembly}")
    logger.info(f"Output prefix: {args.prefix}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Gap size: {args.gap} bp")
    logger.info(f"Number of chromosomes: {args.num_chromosomes}")

    # 验证染色体数量|Validate chromosome number
    if args.num_chromosomes <= 0:
        logger.critical("Number of chromosomes must be greater than 0")
        sys.exit(1)

    try:
        # 检查输入文件|Check input file
        if not os.path.exists(args.assembly):
            logger.critical(f"Input assembly file not found: {args.assembly}")
            sys.exit(1)

        # 创建输出目录|Create output directory
        output_dir = Path(args.output_dir)
        if not output_dir.exists():
            logger.info(f"Creating output directory: {args.output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)

        # 确定输出文件路径|Determine output file paths
        prefix_clean = args.prefix.replace(".assembly", "").replace(".agp", "")
        agp_output = os.path.join(args.output_dir, f"{prefix_clean}.agp")
        chr_list_output = os.path.join(args.output_dir, f"{prefix_clean}.chr.list")

        # 检查输出文件是否存在|Check if output files exist
        if not args.force:
            if os.path.exists(agp_output):
                logger.error(f"AGP output file already exists: {agp_output}")
                logger.error("Use --force to overwrite")
                sys.exit(1)
            if os.path.exists(chr_list_output):
                logger.error(f"Chr.list output file already exists: {chr_list_output}")
                logger.error("Use --force to overwrite")
                sys.exit(1)

        # 转换assembly为AGP|Convert assembly to AGP
        logger.info("=" * 60)
        logger.info("STEP 1: Assembly to AGP Conversion")
        logger.info("=" * 60)
        agp_file, agp_list = assembly2agp(
            assembly=args.assembly,
            agp_output=agp_output,
            gap=args.gap,
            num_chromosomes=args.num_chromosomes,
            logger=logger
        )

        # 生成chr.list文件|Generate chr.list file
        logger.info("=" * 60)
        logger.info("STEP 2: Generate Chromosome List")
        logger.info("=" * 60)
        chr_list_file = generate_chr_list(
            agp_file=agp_file,
            chr_list_output=chr_list_output,
            num_chromosomes=args.num_chromosomes,
            logger=logger
        )

        # 输出总结信息|Output summary information
        elapsed_time = time.time() - start_time
        logger.info("=" * 60)
        logger.info("Conversion Summary")
        logger.info("=" * 60)
        logger.info(f"Total runtime: {elapsed_time:.2f} seconds")
        logger.info(f"AGP output: {agp_file}")
        logger.info(f"Chr.list output: {chr_list_file}")
        logger.info("Conversion completed successfully")

    except KeyboardInterrupt:
        logger.warning("Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.critical(f"Conversion failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
