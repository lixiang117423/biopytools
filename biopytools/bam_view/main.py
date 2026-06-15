"""
BAM比对可视化主程序模块|BAM Alignment Visualization Main Module
"""

import argparse
import sys
from .config import BamViewConfig
from .utils import BamViewLogger, CommandRunner
from .calculator import BamViewCalculator


class BamViewer:
    """BAM比对可视化主类|BAM Alignment Visualization Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = BamViewConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = BamViewLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化计算器|Initialize calculator
        self.calculator = BamViewCalculator(self.config, self.logger, self.cmd_runner)

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            self.logger.info("开始BAM比对可视化分析|Starting BAM alignment visualization analysis")

            # 打印配置信息|Print configuration information
            self._print_config()

            # 生成可视化|Generate visualization
            success = self.calculator.generate_visualization()

            if success:
                self.logger.info("BAM比对可视化分析完成|BAM alignment visualization analysis completed")
                return True
            else:
                self.logger.error("BAM比对可视化分析失败|BAM alignment visualization analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False

    def _print_config(self):
        """打印配置信息|Print configuration information"""
        self.logger.info("=" * 60)
        self.logger.info("配置信息|Configuration Information")
        self.logger.info("=" * 60)
        self.logger.info(f"BAM文件|BAM file: {self.config.bam_file}")
        self.logger.info(f"参考序列|Reference: {self.config.reference}")
        self.logger.info(f"区域|Region: {self.config.region}")
        self.logger.info(f"输出格式|Output format: {self.config.output_format}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"最大reads深度|Max read depth: {self.config.max_read_depth}")
        self.logger.info(f"最大宽度|Max width: {self.config.max_width}")

        if self.config.vcf_file:
            self.logger.info(f"VCF文件|VCF file: {self.config.vcf_file}")

        if self.config.bed_file:
            self.logger.info(f"BED文件|BED file: {self.config.bed_file}")

        if self.config.highlight_intervals:
            self.logger.info(f"高亮区间|Highlight intervals: {', '.join(self.config.highlight_intervals)}")

        if self.config.aux_tags:
            self.logger.info(f"辅助标签|Aux tags: {', '.join(self.config.aux_tags)}")

        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='BAM比对可视化工具|BAM Alignment Visualization Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-b', '--bam',
                       required=True,
                       help='BAM文件路径|BAM file path')
    parser.add_argument('-r', '--reference',
                       required=True,
                       help='参考序列FASTA文件路径|Reference FASTA file path')
    parser.add_argument('-g', '--region',
                       required=True,
                       help='可视化区域(格式: chr:start-end)|Visualization region (format: chr:start-end)')

    # 软件配置|Software configuration
    parser.add_argument('--alignoth-path',
                       default='~/miniforge3/envs/alignoth/bin/alignoth',
                       help='alignoth软件路径|alignoth software path')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir',
                       default='./bam_view_output',
                       help='输出目录|Output directory')
    parser.add_argument('-f', '--output-format',
                       default='html',
                       choices=['html', 'json', 'svg', 'pdf'],
                       help='输出格式|Output format')

    # 可视化参数|Visualization parameters
    parser.add_argument('-d', '--max-read-depth',
                       type=int,
                       default=500,
                       help='最大reads显示深度|Maximum read depth to display')
    parser.add_argument('-w', '--max-width',
                       type=int,
                       default=1024,
                       help='最大宽度|Maximum width')
    parser.add_argument('--mismatch-display-min-percent',
                       type=float,
                       default=1.0,
                       help='显示错配的最小百分比|Minimum percentage of mismatches to display')

    # 高亮选项|Highlight options
    parser.add_argument('-v', '--vcf',
                       help='VCF文件路径(高亮变异位点)|VCF file path (highlight variants)')
    parser.add_argument('--bed',
                       help='BED文件路径(高亮区域)|BED file path (highlight regions)')
    parser.add_argument('-H', '--highlight',
                       dest='highlights',
                       action='append',
                       help='高亮区间(可多次使用, 格式: name:start-end)|'
                            'Highlight interval (can be used multiple times, format: name:start-end)')

    # 其他选项|Other options
    parser.add_argument('-x', '--aux-tag',
                       dest='aux_tags',
                       action='append',
                       help='辅助标签(可多次使用)|Auxiliary tag (can be used multiple times)')
    parser.add_argument('--no-embed-js',
                       action='store_true',
                       help='不嵌入JavaScript(仅HTML格式)|Do not embed JavaScript (HTML format only)')
    parser.add_argument('--plot-all',
                       action='store_true',
                       help='绘制所有reads|Plot all reads')

    args = parser.parse_args()

    # 创建可视化器并运行|Create viewer and run
    viewer = BamViewer(
        bam_file=args.bam,
        reference=args.reference,
        region=args.region,
        alignoth_path=args.alignoth_path,
        output_dir=args.output_dir,
        output_format=args.output_format,
        max_read_depth=args.max_read_depth,
        max_width=args.max_width,
        mismatch_display_min_percent=args.mismatch_display_min_percent,
        vcf_file=args.vcf,
        bed_file=args.bed,
        highlight_intervals=args.highlights,
        aux_tags=args.aux_tags,
        no_embed_js=args.no_embed_js,
        plot_all=args.plot_all
    )

    success = viewer.run_analysis()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
