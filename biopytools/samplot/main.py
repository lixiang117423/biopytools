"""
Samplot Main Module
Samplot主模块
"""

import os
import sys
import time
import argparse
import glob

from .config import SamplotPlotConfig, SamplotVcfConfig
from .utils import SamplotLogger, SamplotRunner


class SamplotPlotter:
    """Samplot plot分析器|Samplot plot analyzer"""

    def __init__(self, **kwargs):
        """
        初始化分析器|Initialize analyzer

        Args:
            **kwargs: SamplotPlotConfig配置参数|SamplotPlotConfig parameters
        """
        self.config = SamplotPlotConfig(**kwargs)
        self.config.validate()

        self.logger_manager = SamplotLogger(self.config.output_dir, log_prefix="samplot_plot")
        self.logger = self.logger_manager.get_logger()
        self.runner = SamplotRunner(self.logger)

    def run(self) -> bool:
        """
        运行分析|Run analysis

        Returns:
            是否成功|Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("Samplot plot 结构变异可视化|Samplot plot SV visualization")
            self.logger.info("=" * 60)

            # 构建命令|Build command
            cmd = self.runner.build_plot_command(self.config)
            success, stdout, stderr = self.runner.run_command(
                cmd, "samplot plot"
            )

            if not success:
                return False

            # 检查输出|Check output
            elapsed = time.time() - start_time
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("分析总结|Analysis Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"染色体|Chromosome: {self.config.chrom}")
            self.logger.info(f"区域|Region: {self.config.start}-{self.config.end}")
            if self.config.sv_type:
                self.logger.info(f"SV类型|SV type: {self.config.sv_type}")
            self.logger.info(f"BAM文件数|BAM file count: {len(self.config.bams)}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"运行时间|Elapsed time: {elapsed:.2f}s")
            self.logger.info("分析完成|Analysis complete")
            self.logger.info("=" * 60)

            return True

        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {str(e)}", exc_info=True)
            return False


class SamplotVcfPlotter:
    """Samplot vcf批量分析器|Samplot vcf batch analyzer"""

    def __init__(self, **kwargs):
        """
        初始化分析器|Initialize analyzer

        Args:
            **kwargs: SamplotVcfConfig配置参数|SamplotVcfConfig parameters
        """
        self.config = SamplotVcfConfig(**kwargs)
        self.config.validate()

        self.logger_manager = SamplotLogger(self.config.output_dir, log_prefix="samplot_vcf")
        self.logger = self.logger_manager.get_logger()
        self.runner = SamplotRunner(self.logger)

    def run(self) -> bool:
        """
        运行分析|Run analysis

        Returns:
            是否成功|Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("Samplot vcf 批量结构变异可视化|Samplot vcf batch SV visualization")
            self.logger.info("=" * 60)

            # 构建命令|Build command
            cmd = self.runner.build_vcf_command(self.config)
            success, stdout, stderr = self.runner.run_command(
                cmd, "samplot vcf"
            )

            if not success:
                return False

            # 统计输出|Count output files
            output_files = self._count_output_files()

            elapsed = time.time() - start_time
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("分析总结|Analysis Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"VCF文件|VCF file: {self.config.vcf}")
            self.logger.info(f"BAM文件数|BAM file count: {len(self.config.bams)}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"生成图片数|Generated plots: {output_files}")
            self.logger.info(f"运行时间|Elapsed time: {elapsed:.2f}s")
            self.logger.info("分析完成|Analysis complete")
            self.logger.info("=" * 60)

            return True

        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {str(e)}", exc_info=True)
            return False

    def _count_output_files(self) -> int:
        """统计输出图片数量|Count output image files"""
        ext = self.config.output_type
        pattern = os.path.join(self.config.output_dir, f"*.{ext}")
        return len(glob.glob(pattern))


def main_plot():
    """samplot plot主函数|samplot plot main function"""
    parser = argparse.ArgumentParser(
        description='Samplot plot 结构变异可视化|Samplot plot SV visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -b sample.bam -c chr1 -s 1000 -e 5000 -t DEL -o output.png
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-b', '--bams', required=True, nargs='+',
                          help='[FILE] BAM/CRAM文件路径|BAM/CRAM file paths')
    required.add_argument('-c', '--chrom', required=True,
                          help='[STR] 染色体名称|Chromosome name')
    required.add_argument('-s', '--start', required=True, type=int,
                          help='[INT] 起始位置|Start position')
    required.add_argument('-e', '--end', required=True, type=int,
                          help='[INT] 结束位置|End position')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-t', '--sv-type', default=None,
                          help='[STR] SV类型(DEL/DUP/INV/BND)|SV type')
    optional.add_argument('-o', '--output-file', default=None,
                          help='[FILE] 输出文件名|Output file name')
    optional.add_argument('--output-dir', default='.',
                          help='[DIR] 输出目录|Output directory')
    optional.add_argument('-r', '--reference', default=None,
                          help='[FILE] 参考基因组(CRAM必须)|Reference genome (required for CRAM)')
    optional.add_argument('-d', '--max-depth', type=int, default=1,
                          help='[INT] 最大正常pair数|Max normal pairs to plot')
    optional.add_argument('-w', '--window', type=int, default=None,
                          help='[INT] 窗口大小|Window size')
    optional.add_argument('-z', '--z', type=int, default=4,
                          help='[INT] 标准差倍数|Number of stdevs from mean')
    optional.add_argument('-H', '--plot-height', type=int, default=None,
                          help='[INT] 图高|Plot height')
    optional.add_argument('-W', '--plot-width', type=int, default=8,
                          help='[INT] 图宽|Plot width')
    optional.add_argument('--dpi', type=int, default=300,
                          help='[INT] DPI|Dots per inch')
    optional.add_argument('--long-read', type=int, default=1000,
                          help='[INT] 长读长最小长度|Min length for long-read')
    optional.add_argument('--coverage-only', action='store_true',
                          help='[FLAG] 仅显示覆盖度|Show only coverage')
    optional.add_argument('--same-yaxis-scales', action='store_true',
                          help='[FLAG] 统一Y轴|Use same Y-axis scales')
    optional.add_argument('-n', '--titles', nargs='+', default=None,
                          help='[STR] 样本标题列表|Sample title list')
    optional.add_argument('--samplot-path',
                          default='~/miniforge3/envs/samplot_v.1.3.0/bin/samplot',
                          help='[FILE] samplot可执行文件路径|samplot binary path')

    args = parser.parse_args()

    # 创建输出目录|Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # 设置日志|Setup logging
    logger_manager = SamplotLogger(args.output_dir, log_prefix="samplot_plot")
    logger = logger_manager.get_logger()

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Samplot plot 结构变异可视化|Samplot plot SV visualization")
        logger.info("=" * 60)

        config = SamplotPlotConfig(
            bams=args.bams,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
            sv_type=args.sv_type,
            output_file=args.output_file,
            output_dir=args.output_dir,
            reference=args.reference,
            max_depth=args.max_depth,
            window=args.window,
            z=args.z,
            plot_height=args.plot_height,
            plot_width=args.plot_width,
            dpi=args.dpi,
            long_read=args.long_read,
            coverage_only=args.coverage_only,
            same_yaxis_scales=args.same_yaxis_scales,
            titles=args.titles,
            samplot_path=args.samplot_path,
        )

        config.validate()

        runner = SamplotRunner(logger)
        cmd = runner.build_plot_command(config)
        success, stdout, stderr = runner.run_command(cmd, "samplot plot")

        if not success:
            logger.error("samplot plot 运行失败|samplot plot failed")
            return 1

        elapsed = time.time() - start_time
        logger.info("")
        logger.info("=" * 60)
        logger.info("分析完成|Analysis complete")
        logger.info(f"运行时间|Elapsed time: {elapsed:.2f}s")
        logger.info("=" * 60)

        return 0

    except KeyboardInterrupt:
        logger.warning("用户中断操作|Operation interrupted by user")
        return 130
    except ValueError as e:
        logger.error(f"配置错误|Configuration error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"分析失败|Analysis failed: {str(e)}", exc_info=True)
        return 1


def main_vcf():
    """samplot vcf主函数|samplot vcf main function"""
    parser = argparse.ArgumentParser(
        description='Samplot vcf 批量结构变异可视化|Samplot vcf batch SV visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s --vcf variants.vcf -b sample1.bam sample2.bam -d output/
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-b', '--bams', required=True, nargs='+',
                          help='[FILE] BAM/CRAM文件路径|BAM/CRAM file paths')
    required.add_argument('--vcf', required=True,
                          help='[FILE] VCF文件路径|VCF file path')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-d', '--output-dir', default='samplot-out',
                          help='[DIR] 输出目录|Output directory')
    optional.add_argument('-O', '--output-type', default='png',
                          choices=['png', 'pdf', 'eps', 'jpg'],
                          help='[STR] 输出格式|Output format')
    optional.add_argument('-t', '--threads', type=int, default=1,
                          help='[INT] 线程数|Number of threads')
    optional.add_argument('--downsample', type=int, default=1,
                          help='[INT] 下采样数|Downsample count')
    optional.add_argument('--min-bp', type=int, default=20,
                          help='[INT] 最小SV长度(bp)|Min SV length in bp')
    optional.add_argument('--max-mb', type=int, default=None,
                          help='[INT] 最大SV长度(MB)|Max SV length in MB')
    optional.add_argument('--sample-ids', nargs='+', default=None,
                          help='[STR] 样本ID列表|Sample ID list')
    optional.add_argument('--plot-all', action='store_true',
                          help='[FLAG] 绘制所有样本|Plot all samples')
    optional.add_argument('--min-call-rate', type=float, default=None,
                          help='[FLOAT] 最小call rate|Min call rate')
    optional.add_argument('--max-hets', type=int, default=None,
                          help='[INT] 最大杂合数|Max heterozygotes')
    optional.add_argument('--min-entries', type=int, default=6,
                          help='[INT] 最小样本数|Min entries to plot')
    optional.add_argument('--max-entries', type=int, default=10,
                          help='[INT] 最大样本数|Max entries to plot')
    optional.add_argument('--gff3', default=None,
                          help='[FILE] GFF3注释文件|GFF3 annotation file')
    optional.add_argument('--samplot-path',
                          default='~/miniforge3/envs/samplot_v.1.3.0/bin/samplot',
                          help='[FILE] samplot可执行文件路径|samplot binary path')

    args = parser.parse_args()

    # 创建输出目录|Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # 设置日志|Setup logging
    logger_manager = SamplotLogger(args.output_dir, log_prefix="samplot_vcf")
    logger = logger_manager.get_logger()

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Samplot vcf 批量结构变异可视化|Samplot vcf batch SV visualization")
        logger.info("=" * 60)

        config = SamplotVcfConfig(
            bams=args.bams,
            vcf=args.vcf,
            output_dir=args.output_dir,
            output_type=args.output_type,
            threads=args.threads,
            downsample=args.downsample,
            min_bp=args.min_bp,
            max_mb=args.max_mb,
            sample_ids=args.sample_ids,
            plot_all=args.plot_all,
            min_call_rate=args.min_call_rate,
            max_hets=args.max_hets,
            min_entries=args.min_entries,
            max_entries=args.max_entries,
            gff3=args.gff3,
            samplot_path=args.samplot_path,
        )

        config.validate()

        runner = SamplotRunner(logger)
        cmd = runner.build_vcf_command(config)
        success, stdout, stderr = runner.run_command(cmd, "samplot vcf")

        if not success:
            logger.error("samplot vcf 运行失败|samplot vcf failed")
            return 1

        elapsed = time.time() - start_time
        logger.info("")
        logger.info("=" * 60)
        logger.info("分析完成|Analysis complete")
        logger.info(f"运行时间|Elapsed time: {elapsed:.2f}s")
        logger.info("=" * 60)

        return 0

    except KeyboardInterrupt:
        logger.warning("用户中断操作|Operation interrupted by user")
        return 130
    except ValueError as e:
        logger.error(f"配置错误|Configuration error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"分析失败|Analysis failed: {str(e)}", exc_info=True)
        return 1
