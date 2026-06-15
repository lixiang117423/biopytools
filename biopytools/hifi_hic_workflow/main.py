"""
HiFi+Hi-C工作流主入口|HiFi+Hi-C Workflow Main Entry

完整的植物基因组组装流程：HiFi组装 → Hi-C挂载 → 染色体重命名 → Hi-C热图
Complete plant genome assembly workflow: HiFi assembly → Hi-C scaffolding → Chromosome rename → Hi-C heatmap
"""

import argparse
import sys
import os
import logging
from datetime import datetime
from pathlib import Path

# 添加当前目录到sys.path|Add current directory to sys.path
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

try:
    from config import HifiHicWorkflowConfig
    from workflow import HifiHicWorkflow
except ImportError:
    # 如果相对导入失败，尝试绝对导入|Try absolute import if relative import fails
    sys.path.append(str(current_dir.parent))
    from biopytools.hifi_hic_workflow.config import HifiHicWorkflowConfig
    from biopytools.hifi_hic_workflow.workflow import HifiHicWorkflow


def setup_logger(output_dir: str, verbose: bool = False) -> logging.Logger:
    """
    设置日志记录器|Setup logger

    Args:
        output_dir: 输出目录|Output directory
        verbose: 是否显示详细日志|Whether to show verbose logs

    Returns:
        logging.Logger: 日志记录器|Logger
    """
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"workflow_{timestamp}.log"

    # 创建logger|Create logger
    logger = logging.getLogger("hifi_hic_workflow")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # 日志格式|Log format
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 控制台处理器|Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # 文件处理器|File handler
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='HiFi+Hi-C基因组组装与挂载流程|HiFi+Hi-C Genome Assembly and Scaffolding Workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本流程|Basic workflow
  %(prog)s --hifi hifi_reads.fq.gz --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz \\
           --ref reference.fa -o ./workflow_output

  # 使用NGS polish|With NGS polish
  %(prog)s --hifi hifi_reads.fq.gz --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz \\
           --ref reference.fa --ngs ngs_data_dir -o ./workflow_output

  # 跳过某些步骤|Skip certain steps
  %(prog)s --hifi hifi_reads.fq.gz --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz \\
           --ref reference.fa -o ./workflow_output --skip-heatmap
        '''
    )

    # ==================== 必需参数 | Required Parameters ====================
    required_group = parser.add_argument_group('必需参数|Required parameters')
    required_group.add_argument(
        '--hifi',
        required=True,
        help='HiFi reads文件|HiFi reads file'
    )
    required_group.add_argument(
        '--hic-r1',
        required=True,
        help='Hi-C R1文件|Hi-C R1 file'
    )
    required_group.add_argument(
        '--hic-r2',
        required=True,
        help='Hi-C R2文件|Hi-C R2 file'
    )
    required_group.add_argument(
        '--ref', '--reference',
        required=True,
        help='参考基因组FASTA文件（仅用于命名）|Reference genome FASTA file (for naming only)'
    )
    required_group.add_argument(
        '-o', '--output',
        required=True,
        help='输出目录|Output directory'
    )

    # ==================== 全局参数 | Global Parameters ====================
    global_group = parser.add_argument_group('全局参数|Global parameters')
    global_group.add_argument(
        '-p', '--prefix',
        default='genome_sample',
        help='样本前缀|Sample prefix (default: genome_sample)'
    )
    global_group.add_argument(
        '-t', '--threads',
        type=int,
        default=64,
        help='线程数|Number of threads (default: 64)'
    )

    # ==================== 流程控制 | Workflow Control ====================
    control_group = parser.add_argument_group('流程控制|Workflow control')
    control_group.add_argument(
        '--skip-hifi-hic',
        action='store_true',
        help='跳过HiFi组装|Skip HiFi assembly'
    )
    control_group.add_argument(
        '--skip-haphic',
        action='store_true',
        help='跳过Hi-C挂载|Skip Hi-C scaffolding'
    )
    control_group.add_argument(
        '--skip-rename',
        action='store_true',
        help='跳过重命名|Skip renaming'
    )
    control_group.add_argument(
        '--skip-heatmap',
        action='store_true',
        help='跳过热图|Skip heatmap'
    )
    control_group.add_argument(
        '--no-resume',
        action='store_true',
        help='禁用断点续传|Disable resume mode'
    )
    control_group.add_argument(
        '--force',
        action='store_true',
        help='强制重新运行所有步骤|Force rerun all steps'
    )

    # ==================== Step 1: HiFi组装参数 | Step 1: HiFi Assembly Parameters ====================
    hifi_group = parser.add_argument_group('HiFi组装参数|HiFi assembly parameters')
    hifi_group.add_argument(
        '--genome-size',
        default='1.45g',
        help='预估基因组大小|Estimated genome size (default: 1.45g)'
    )
    hifi_group.add_argument(
        '--n-hap',
        type=int,
        default=2,
        help='倍性|Ploidy (default: 2)'
    )
    hifi_group.add_argument(
        '--purge-level',
        type=int,
        help='Purge level (0=no purging, 1=light, 2/3=aggressive)'
    )
    hifi_group.add_argument(
        '--hom-cov',
        type=int,
        help='Homozygous read coverage (--hom-cov)'
    )

    # NGS Polish参数|NGS polish parameters
    ngs_group = parser.add_argument_group('NGS Polish参数|NGS polish parameters')
    ngs_group.add_argument(
        '--use-ngs-polish',
        action='store_true',
        help='启用NGS polish|Enable NGS polish'
    )
    ngs_group.add_argument(
        '--ngs-data',
        help='NGS二代数据目录|NGS second-generation data directory'
    )
    ngs_group.add_argument(
        '--ngs-high-cov',
        type=float,
        default=95.0,
        help='高质量contig覆盖度阈值|High quality contig coverage threshold (default: 95.0)'
    )
    ngs_group.add_argument(
        '--ngs-pattern',
        default='_1.clean.fq.gz',
        help='NGS文件匹配模式|NGS file matching pattern (default: _1.clean.fq.gz)'
    )

    # ==================== Step 2: HapHiC参数 | Step 2: HapHiC Parameters ====================
    haphic_group = parser.add_argument_group('HapHiC挂载参数|HapHiC scaffolding parameters')
    haphic_group.add_argument(
        '--nchrs',
        type=int,
        help='染色体数量（如不指定，从reference统计）|Number of chromosomes (count from reference if not specified)'
    )

    # 工具路径|Tool paths
    tool_group = parser.add_argument_group('HapHiC工具路径|HapHiC tool paths')
    tool_group.add_argument(
        '--haphic-bin',
        help='HapHiC可执行文件路径|HapHiC executable path'
    )
    tool_group.add_argument(
        '--bwa-bin',
        help='BWA可执行文件路径|BWA executable path'
    )
    tool_group.add_argument(
        '--samtools-bin',
        help='Samtools可执行文件路径|Samtools executable path'
    )

    # ==================== Step 3: 染色体重命名参数 | Step 3: Chromosome Rename Parameters ====================
    rename_group = parser.add_argument_group('染色体重命名参数|Chromosome rename parameters')
    rename_group.add_argument(
        '--rename-keep-all',
        action='store_true',
        default=True,
        help='保留所有序列（chrNN + scaffolds）|Keep all sequences (default: True)'
    )

    # 参考基因组命名参数|Reference genome naming parameters
    naming_group = parser.add_argument_group('参考基因组命名参数|Reference genome naming parameters')
    naming_group.add_argument(
        '--naming-min-identity',
        type=float,
        default=80.0,
        help='最小序列一致性|Min sequence identity %% (default: 80.0)'
    )
    naming_group.add_argument(
        '--naming-min-coverage',
        type=float,
        default=80.0,
        help='最小覆盖度|Min coverage %% (default: 80.0)'
    )
    naming_group.add_argument(
        '--naming-minimap2-preset',
        default='asm5',
        choices=['asm5', 'asm10', 'asm20'],
        help='minimap2预设|minimap2 preset (default: asm5)'
    )

    # ==================== Step 4: Hi-C热图参数 | Step 4: Hi-C Heatmap Parameters ====================
    heatmap_group = parser.add_argument_group('Hi-C热图参数|Hi-C heatmap parameters')
    heatmap_group.add_argument(
        '--hicpro-enzyme',
        default='MboI',
        help='HiCPro限制性内切酶|HiCPro restriction enzyme (default: MboI)'
    )
    heatmap_group.add_argument(
        '--heatmap-resolution',
        type=int,
        default=100000,
        help='热图分辨率|Heatmap resolution in bp (default: 100000)'
    )
    heatmap_group.add_argument(
        '--heatmap-colormap',
        default='YlOrRd',
        help='热图颜色方案|Heatmap color scheme (default: YlOrRd)'
    )
    heatmap_group.add_argument(
        '--heatmap-format',
        default='pdf',
        help='热图输出格式|Heatmap output format (default: pdf)'
    )

    # ==================== 其他参数 | Other Parameters ====================
    other_group = parser.add_argument_group('其他参数|Other parameters')
    other_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='显示详细日志|Show verbose logs'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 设置日志|Setup logging
        logger = setup_logger(args.output, args.verbose)

        logger.info("=" * 80)
        logger.info("HiFi+Hi-C基因组组装与挂载流程|HiFi+Hi-C Genome Assembly and Scaffolding Workflow")
        logger.info("=" * 80)

        # 创建配置|Create configuration
        config = HifiHicWorkflowConfig(
            # 必需参数|Required parameters
            hifi_reads=args.hifi,
            hic_r1=args.hic_r1,
            hic_r2=args.hic_r2,
            reference_genome=args.ref,
            output_dir=args.output,

            # 全局参数|Global parameters
            prefix=args.prefix,
            threads=args.threads,

            # 流程控制|Workflow control
            skip_hifi_hic=args.skip_hifi_hic,
            skip_haphic=args.skip_haphic,
            skip_rename=args.skip_rename,
            skip_heatmap=args.skip_heatmap,
            resume=not args.no_resume,
            force_rerun=args.force,

            # Step 1: HiFi组装参数|Step 1: HiFi assembly parameters
            genome_size=args.genome_size,
            n_hap=args.n_hap,
            purge_level=args.purge_level,
            hom_cov=args.hom_cov,

            # NGS Polish参数|NGS polish parameters
            use_ngs_polish=args.use_ngs_polish,
            ngs_data=args.ngs_data,
            ngs_high_cov=args.ngs_high_cov,
            ngs_pattern=args.ngs_pattern,

            # Step 2: HapHiC参数|Step 2: HapHiC parameters
            nchrs=args.nchrs,
            haphic_bin=args.haphic_bin or 'haphic',
            bwa_bin=args.bwa_bin or 'bwa',
            samtools_bin=args.samtools_bin or 'samtools',

            # Step 3: 染色体重命名参数|Step 3: Chromosome rename parameters
            rename_keep_all=args.rename_keep_all,
            naming_min_identity=args.naming_min_identity,
            naming_min_coverage=args.naming_min_coverage,
            naming_minimap2_preset=args.naming_minimap2_preset,

            # Step 4: Hi-C热图参数|Step 4: Hi-C heatmap parameters
            hicpro_restriction_enzyme=args.hicpro_enzyme,
            heatmap_resolution=args.heatmap_resolution,
            heatmap_color_map=args.heatmap_colormap,
            heatmap_format=args.heatmap_format,
        )

        # 验证配置|Validate configuration
        config.validate()

        # 运行工作流|Run workflow
        workflow = HifiHicWorkflow(config, logger)
        success = workflow.run_workflow()

        if success:
            logger.info("\n" + "=" * 80)
            logger.info("工作流成功完成！|Workflow completed successfully!")
            logger.info("=" * 80)
            return 0
        else:
            logger.error("\n" + "=" * 80)
            logger.error("工作流执行失败|Workflow execution failed")
            logger.error("=" * 80)
            return 1

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        return 130
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        return 1
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
