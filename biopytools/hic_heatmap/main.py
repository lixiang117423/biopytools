"""Hi-C热图分析主函数|Hi-C heatmap analysis main function"""

import argparse
import sys
from pathlib import Path

from .config import HiCProConfig
from .utils import HiCLogger
from .hicpro_pipeline import HiCProPipeline
from ..common.paths import expand_path


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Hi-C全基因组热图分析工具 (HiCPro + PlotHiC)|Hi-C whole genome heatmap analysis tool (HiCPro + PlotHiC)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--genome',
        required=True,
        help='基因组FASTA文件|Genome FASTA file'
    )

    parser.add_argument(
        '-g', '--genome-id',
        required=True,
        help='基因组ID（用于输出文件命名，如hg19, mm10）|Genome ID (for output file naming, e.g., hg19, mm10)'
    )

    parser.add_argument(
        '-o', '--output-dir',
        default='./hic_output',
        help='输出目录|Output directory (default: ./hic_output)'
    )

    # 输入文件选项|Input file options
    input_group = parser.add_argument_group('输入文件选项（必需）|Input file options (required)')

    input_group.add_argument(
        '-1', '--fastq-r1',
        required=True,
        help='R1测序文件|R1 sequencing file'
    )

    input_group.add_argument(
        '-2', '--fastq-r2',
        required=True,
        help='R2测序文件|R2 sequencing file'
    )

    # 可选参数|Optional parameters
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=64,
        help='线程数|Threads (default: 64)'
    )

    parser.add_argument(
        '--max-memory',
        type=int,
        default=200,
        help='HiC-Pro最大内存限制（GB）|Maximum memory limit for HiC-Pro in GB (default: 200)'
    )

    # HiCPro参数|HiCPro parameters
    hicpro_group = parser.add_argument_group('HiCPro参数|HiCPro parameters')

    hicpro_group.add_argument(
        '--restriction-enzyme',
        default='MboI',
        help='限制性内切酶名称|Restriction enzyme name (default: MboI). Options: MboI, HindIII, NcoI, EcoRI, BamHI'
    )

    hicpro_group.add_argument(
        '--bowtie2-idx',
        help='Bowtie2索引路径（默认自动生成）|Bowtie2 index path (auto-generated if not specified)'
    )

    hicpro_group.add_argument(
        '--bin-sizes',
        default='20000 40000 150000 500000 1000000',
        help='Contact map bin大小（空格分隔）|Contact map bin sizes, space-separated (default: "20000 40000 150000 500000 1000000")'
    )

    # PlotHiC参数|PlotHiC parameters
    plot_group = parser.add_argument_group('PlotHiC参数|PlotHiC parameters')

    plot_group.add_argument(
        '--resolution',
        type=int,
        default=100000,
        help='热图分辨率（bp）|Heatmap resolution in bp (default: 100000, 100kb)'
    )

    plot_group.add_argument(
        '--color-map',
        default='YlOrRd',
        help='热图颜色方案|Heatmap color scheme (default: YlOrRd, PlotHiC default)'
    )

    plot_group.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='图像分辨率|Image resolution in DPI (default: 300)'
    )

    plot_group.add_argument(
        '--format',
        default='pdf',
        help='输出格式|Output format (default: pdf, options: pdf, png, svg, etc.)'
    )

    plot_group.add_argument(
        '--bar-max',
        type=int,
        default=1,
        help='颜色条最大值（log变换后）|Color bar maximum value after log transform (default: 1)'
    )

    # 工具路径参数|Tool path parameters
    tool_group = parser.add_argument_group('工具路径参数|Tool path parameters')

    tool_group.add_argument(
        '--hicpro-sif',
        default='~/software/singularity/hicpro_latest.sif',
        help='HiCPro singularity镜像路径|HiCPro singularity image path'
    )

    tool_group.add_argument(
        '--singularity-exec',
        default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
        help='Singularity可执行文件路径|Singularity executable path'
    )

    tool_group.add_argument(
        '--plothic-path',
        default='~/miniforge3/envs/plothic_v.1.0.0/bin/plothic',
        help='PlotHiC可执行文件路径|PlotHiC executable path'
    )

    # 流程控制参数|Process control parameters
    control_group = parser.add_argument_group('流程控制参数|Process control parameters')

    control_group.add_argument(
        '--force',
        action='store_true',
        help='强制重新运行所有步骤|Force rerun all steps'
    )

    control_group.add_argument(
        '--verbose',
        action='store_true',
        help='显示详细日志|Show verbose logs'
    )

    control_group.add_argument(
        '--quiet',
        action='store_true',
        help='仅显示错误日志|Show error logs only'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = HiCProConfig(
            genome=args.genome,
            genome_id=args.genome_id,
            output_dir=args.output_dir,
            fastq_r1=args.fastq_r1,
            fastq_r2=args.fastq_r2,
            threads=args.threads,
            max_memory_gb=args.max_memory,
            restriction_enzyme=args.restriction_enzyme,
            bowtie2_idx=args.bowtie2_idx,
            bin_sizes=args.bin_sizes,
            resolution=args.resolution,
            color_map=args.color_map,
            dpi=args.dpi,
            output_format=args.format,
            bar_max=args.bar_max,
            skip_existing=not args.force,
            force_hicpro=args.force,
            hicpro_sif=expand_path(args.hicpro_sif) if args.hicpro_sif else None,
            singularity_exec=expand_path(args.singularity_exec),
            plothic_path=expand_path(args.plothic_path)
        )

        # 验证配置|Validate configuration
        config.validate()

        # 设置日志级别|Setup log level
        if args.verbose:
            log_level = "DEBUG"
        elif args.quiet:
            log_level = "ERROR"
        else:
            log_level = "INFO"

        log_file = Path(config.output_dir) / "hic_heatmap.log"
        logger_manager = HiCLogger(log_file=str(log_file), log_level=log_level)
        logger = logger_manager.get_logger()

        # 运行分析流程|Run analysis pipeline
        pipeline = HiCProPipeline(config, logger)
        success = pipeline.run()

        if success:
            logger.info("分析成功完成|Analysis completed successfully")
            sys.exit(0)
        else:
            logger.error("分析失败|Analysis failed")
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(130)

    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
