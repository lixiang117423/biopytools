"""
BRAKER3基因组注释主程序模块|BRAKER3 Genome Annotation Main Module
"""

import argparse
import sys
import os
import tempfile
import re
from pathlib import Path


def clean_protein_sequences(input_fasta, logger=None):
    """
    清理蛋白质序列文件，移除非法字符
    Clean protein sequence file, remove invalid characters

    Args:
        input_fasta: 输入蛋白质FASTA文件|Input protein FASTA file
        logger: 日志器|Logger object

    Returns:
        str: 清理后的FASTA文件路径|Path to cleaned FASTA file
    """
    # 创建临时文件（放在输出目录中，避免被删除）|Create temp file in output directory
    # 使用输入文件的目录作为输出位置
    input_dir = os.path.dirname(input_fasta)
    input_basename = os.path.basename(input_fasta)
    cleaned_path = os.path.join(input_dir, f"cleaned_{input_basename}")

    removed_chars = {'dots': 0, 'digits': 0, 'others': 0, 'total_sequences': 0}

    try:
        with open(input_fasta, 'r') as infile, open(cleaned_path, 'w') as outfile:
            current_seq = []
            header = None

            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    # 写入前一个序列（如果有）|Write previous sequence if exists
                    if header is not None:
                        # 清理序列：移除除标准20个氨基酸和X之外的所有字符
                        # Clean sequence: remove all characters except standard 20 amino acids and X
                        cleaned_seq = ''.join(current_seq)

                        # 统计移除的字符|Count removed characters
                        removed_chars['dots'] += cleaned_seq.count('.')
                        removed_chars['digits'] += len(re.sub(r'[0-9]', '', cleaned_seq)) - len(cleaned_seq)

                        cleaned_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWXY]', '', cleaned_seq.upper())

                        removed_chars['total_sequences'] += 1

                        outfile.write(f"{header}\n{cleaned_seq}\n")

                    # 保存新header|Save new header
                    header = line
                    current_seq = []
                else:
                    current_seq.append(line)

            # 写入最后一个序列|Write last sequence
            if header is not None:
                cleaned_seq = ''.join(current_seq)
                removed_chars['dots'] += cleaned_seq.count('.')
                cleaned_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWXY]', '', cleaned_seq.upper())
                removed_chars['total_sequences'] += 1
                outfile.write(f"{header}\n{cleaned_seq}\n")

        if logger:
            logger.info(f"蛋白质序列清理统计|Cleaning statistics:")
            logger.info(f"  序列数|Sequences: {removed_chars['total_sequences']}")
            logger.info(f"  移除的 '.' 字符|Removed dots: {removed_chars['dots']}")
            logger.info(f"  输出文件|Output: {cleaned_path}")

        return cleaned_path

    except Exception as e:
        # 清理临时文件|Clean up temp file
        if os.path.exists(cleaned_path):
            os.remove(cleaned_path)
        raise RuntimeError(f"清理蛋白质文件失败|Failed to clean protein file: {e}")


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='BRAKER3基因组注释工具|BRAKER3 Genome Annotation Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法
  %(prog)s --genome genome.fa --species my_oomycete --prot_seq proteins.fa

  # 使用所有证据类型
  %(prog)s --genome genome.fa --species my_oomycete --prot_seq proteins.fa \\
           --isoseq isoseq.fa --rnaseq_dirs /path/to/rnaseq1,/path/to/rnaseq2 \\
           --rnaseq_sets set1,set2 --threads 40

  # 使用真菌模式
  %(prog)s --genome genome.fa --species my_oomycete --prot_seq proteins.fa --fungus

  # 跳过重复序列屏蔽
  %(prog)s --genome genome.fa --species my_oomycete --bam rnaseq.bam --skip-repeat
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('--genome', '-g',
                       required=True,
                       help='基因组FASTA文件|Genome FASTA file')

    parser.add_argument('--species', '-s',
                       required=True,
                       help='物种名称|Species name (用于BRAKER输出命名)|used for BRAKER output naming')

    # 输入数据|Input data
    parser.add_argument('--prot_seq', '-p',
                       help='近缘物种蛋白质序列文件或文件夹|Protein sequences file or directory')

    parser.add_argument('--isoseq', '-l',
                       help='三代全长转录本文件夹|Long-read transcript directory')

    parser.add_argument('--rnaseq_dirs',
                       help='二代RNA-seq数据目录列表，逗号分隔|Comma-separated list of RNA-seq directories')

    parser.add_argument('--read1_pattern',
                       default='_1.clean.fq.gz',
                       help='R1文件模式|R1 file pattern (default: _1.clean.fq.gz)')

    parser.add_argument('--read2_pattern',
                       default='_2.clean.fq.gz',
                       help='R2文件模式|R2 file pattern (default: _2.clean.fq.gz)')

    # 输出配置|Output configuration
    parser.add_argument('--output_dir', '-o',
                       default='./braker_output',
                       help='输出目录|Output directory (default: ./braker_output)')

    # 流程参数|Pipeline parameters
    parser.add_argument('--threads', '-t',
                       type=int,
                       default=12,
                       help='线程数|Number of threads (default: 12)')

    parser.add_argument('--fungus',
                       action='store_true',
                       help='使用真菌模式|Use fungus mode (suitable for oomycetes)')

    # Singularity配置|Singularity configuration
    parser.add_argument('--singularity_image',
                       default="~/software/singularity/braker3_devel.sif",
                       help='Singularity镜像路径|Singularity image path')

    parser.add_argument('--no_singularity',
                       action='store_true',
                       help='不使用Singularity镜像|Do not use Singularity image')

    # 步骤控制|Step control
    parser.add_argument('--skip_repeat',
                       action='store_true',
                       help='跳过重复序列屏蔽|Skip repeat masking step')

    parser.add_argument('--skip_long_reads',
                       action='store_true',
                       help='跳过三代转录本处理|Skip long-read processing step')

    parser.add_argument('--skip_short_reads',
                       action='store_true',
                       help='跳过二代RNA-seq处理|Skip short-read processing step')

    # BRAKER3特定参数|BRAKER3 specific parameters
    parser.add_argument('--busco_lineage',
                       help='BUSCO谱系|BUSCO lineage')

    parser.add_argument('--utr',
                       action='store_true',
                       help='预测UTR|Predict UTR regions')

    parser.add_argument('--training_genes',
                       help='训练基因集文件|Training gene set file')

    parser.add_argument('--use_existing',
                       action='store_true',
                       help='使用已有参数|Use existing parameters')

    # repeat_refine 参数(repeat库过滤+证据还原)|repeat_refine params
    parser.add_argument('--skip_repeat_filter',
                       action='store_true',
                       help='跳过repeat库过滤(方案1,默认开启)|Skip repeat library filtering')
    parser.add_argument('--skip_rescue',
                       action='store_true',
                       help='跳过证据还原(方案2,默认开启)|Skip evidence-based rescue')
    parser.add_argument('--pfam_db',
                       default="~/database/eggnog/pfam/Pfam-A.hmm",
                       help='Pfam-A HMM 库路径|Pfam-A HMM DB path')
    parser.add_argument('--te_domain_evalue',
                       type=float,
                       default=1e-5,
                       help='TE domain hmmscan E-value 阈值|TE domain E-value cutoff')
    parser.add_argument('--filter_min_orf_len',
                       type=int,
                       default=30,
                       help='过滤用最小ORF长度(aa)|Min ORF length (aa) for filter')
    parser.add_argument('--rescue_min_cds_len',
                       type=int,
                       default=100,
                       help='rescue蛋白证据最小覆盖长度(bp)|Min CDS overlap (bp) for rescue')
    parser.add_argument('--rescue_min_identity',
                       type=float,
                       default=70,
                       help='rescue蛋白最小identity(%%)|Min protein identity (%%) for rescue')
    parser.add_argument('--rescue_min_depth',
                       type=int,
                       default=5,
                       help='rescue RNA-seq最小覆盖度|Min RNA-seq depth for rescue')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    # 创建日志器用于主程序输出|Create logger for main program output
    from .utils import BrakerLogger
    log_file = os.path.join(args.output_dir, "logs", "braker_main.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger_manager = BrakerLogger(log_file)
    logger = logger_manager.get_logger()

    try:
        # 导入模块|Import modules
        from .config import BrakerConfig
        from .pipeline import BrakerPipeline
        from .utils import find_long_reads_in_directory, find_protein_files_in_directory

        # 处理RNA-seq目录列表|Process RNA-seq directory list
        rnaseq_dirs = None
        if args.rnaseq_dirs:
            rnaseq_dirs = [d.strip() for d in args.rnaseq_dirs.split(',')]

        # 处理三代转录组：如果是文件夹，自动识别文件|Process long-reads: auto-detect if directory
        isoseq_file = args.isoseq
        if args.isoseq and os.path.isdir(args.isoseq):
            logger.info(f"检测到三代转录组为目录，自动识别文件|Detected long-read directory, auto-detecting files...")
            isoseq_file = find_long_reads_in_directory(args.isoseq, logger=logger)
            if not isoseq_file:
                logger.error(f"未找到三代转录组文件|No long-read files found in: {args.isoseq}")
                sys.exit(1)
            logger.info(f"找到三代转录组文件|Found long-read file: {isoseq_file}")

        # 处理蛋白质序列：如果是文件夹，自动识别并合并|Process proteins: auto-detect and merge if directory
        prot_seq_file = args.prot_seq
        if args.prot_seq and os.path.isdir(args.prot_seq):
            logger.info(f"检测到蛋白质序列为目录，自动识别并合并文件|Detected protein directory, auto-detecting and merging files...")
            prot_seq_file = find_protein_files_in_directory(args.prot_seq, logger=logger)
            if not prot_seq_file:
                logger.error(f"未找到蛋白质序列文件|No protein files found in: {args.prot_seq}")
                sys.exit(1)
            logger.info(f"蛋白质文件准备完成|Protein files ready: {prot_seq_file}")

        # 清理蛋白质序列，移除非法字符（如 .）|Clean protein sequences, remove invalid characters (e.g., .)
        if prot_seq_file:
            logger.info(f"清理蛋白质序列|Cleaning protein sequences...")
            prot_seq_file = clean_protein_sequences(prot_seq_file, logger=logger)
            logger.info(f"蛋白质序列清理完成|Protein sequences cleaned: {prot_seq_file}")

        # 创建配置|Create configuration
        config = BrakerConfig(
            genome=args.genome,
            species=args.species,
            prot_seq=prot_seq_file,
            isoseq=isoseq_file,
            rnaseq_dirs=rnaseq_dirs,
            rnaseq_sets=None,  # 已废弃，不再使用|Deprecated, no longer used
            read1_pattern=args.read1_pattern,
            read2_pattern=args.read2_pattern,
            use_singularity=not args.no_singularity,
            singularity_image=args.singularity_image,
            output_dir=args.output_dir,
            threads=args.threads,
            use_fungus=args.fungus,
            skip_repeat=args.skip_repeat,
            skip_long_reads=args.skip_long_reads,
            skip_short_reads=args.skip_short_reads,
            busco_lineage=args.busco_lineage,
            utr=args.utr,
            training_genes=args.training_genes,
            use_existing=args.use_existing,
            skip_repeat_filter=args.skip_repeat_filter,
            skip_rescue=args.skip_rescue,
            pfam_db=args.pfam_db,
            te_domain_evalue=args.te_domain_evalue,
            filter_min_orf_len=args.filter_min_orf_len,
            rescue_min_cds_len=args.rescue_min_cds_len,
            rescue_min_identity=args.rescue_min_identity,
            rescue_min_depth=args.rescue_min_depth
        )

        # 验证配置|Validate configuration
        config.validate()

        # 创建流程|Create pipeline
        pipeline = BrakerPipeline(config)

        # 运行流程|Run pipeline
        final_annotation = pipeline.run_pipeline()

        logger.info(f"注释完成|Annotation completed!")
        logger.info(f"最终结果|Final result: {final_annotation}")
        sys.exit(0)

    except ValueError as e:
        logger.error(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except RuntimeError as e:
        logger.error(f"运行错误|Runtime error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"未预期的错误|Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
