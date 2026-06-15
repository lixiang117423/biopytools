"""
GenomeSyn2主程序模块|GenomeSyn2 Main Module
"""

import argparse
import sys
from .config import GenomeSyn2Config
from .utils import GenomeSyn2Logger
from .calculator import GenomeSyn2Calculator


class GenomeSyn2Runner:
    """GenomeSyn2运行器|GenomeSyn2 Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = GenomeSyn2Config(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = GenomeSyn2Logger(
            output_dir=self.config.output_path,
            log_name="genomesyn2.log"
        )
        self.logger = self.logger_manager.get_logger()

        # 检查Perl模块|Check Perl modules
        self._check_perl_modules()

        # 初始化计算器|Initialize calculator
        self.calculator = GenomeSyn2Calculator(self.config, self.logger)

    def _check_perl_modules(self):
        """检查必需的Perl模块|Check required Perl modules"""
        from .utils import check_perl_module

        required_modules = [
            ('Bio::SeqIO', 'Bioperl'),
            ('SVG', 'SVG'),
        ]

        missing_modules = []
        for module_name, display_name in required_modules:
            if not check_perl_module(self.config.perl_path, module_name):
                missing_modules.append(display_name)

        if missing_modules:
            self.logger.warning(
                f"部分Perl模块可能未安装，可能影响运行|"
                f"Some Perl modules may not be installed: {', '.join(missing_modules)}"
            )

    def run(self):
        """运行分析|Run analysis"""
        try:
            self.logger.info("GenomeSyn2开始运行|GenomeSyn2 starting")

            success = self.calculator.run_analysis()

            if success:
                self.logger.info("GenomeSyn2运行完成|GenomeSyn2 completed successfully")
                return 0
            else:
                self.logger.error("GenomeSyn2运行失败|GenomeSyn2 run failed")
                return 1

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return 1


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="GenomeSyn2比较基因组学可视化工具|GenomeSyn2 Comparative Genomics Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基因组比对|Genome alignment
  %(prog)s --align mummer --genome ./genome_dir/ --outdir ./output/ --thread 12

  # 蛋白质比对|Protein alignment
  %(prog)s --align blastp --genome ./genome_dir/ --gene ./gene_dir/ --outdir ./output/

  # 从VCF计算SNP|Calculate SNPs from VCF
  %(prog)s --vcf variants.vcf --bin 50000

  # 绘制血统解析图|Plot ancestry deconvolution
  %(prog)s --identity SNP_identity.bed --density SNP_density.bed

  # 绘制共线性图|Plot synteny diagram
  %(prog)s --conf total.conf

  # 生成文件列表|Generate file list
  %(prog)s --type fa --path ./genome_dir/ --out genome.info.tsv
        """
    )

    # ========== 比对模式参数 | Alignment Mode Parameters ==========
    parser.add_argument(
        '--align',
        choices=['mummer', 'minimap2', 'blastp', 'mmseqs', 'diamond'],
        help='比对软件类型|Alignment software type (mummer/minimap2/blastp/mmseqs/diamond)'
    )

    parser.add_argument(
        '--genome',
        help='基因组文件目录|Genome files directory (文件名需按数字排序|Files must be numbered)'
    )

    parser.add_argument(
        '--gene',
        help='基因注释文件目录|Gene annotation files directory (蛋白质比对需要|required for protein alignment)'
    )

    parser.add_argument(
        '--outdir',
        help='输出目录|Output directory'
    )

    # ========== VCF模式参数 | VCF Mode Parameters ==========
    parser.add_argument(
        '--vcf',
        help='VCF文件路径|VCF file path for SNP analysis'
    )

    parser.add_argument(
        '--bin',
        type=int,
        default=50000,
        help='Bin大小(用于SNP分析)|Bin size for SNP analysis (default: 50000)'
    )

    parser.add_argument(
        '--identity',
        help='SNP一致性文件|SNP identity BED file'
    )

    parser.add_argument(
        '--density',
        help='SNP密度文件|SNP density BED file'
    )

    # ========== 绘图模式参数 | Plotting Mode Parameters ==========
    parser.add_argument(
        '--conf',
        help='配置文件路径|Configuration file path'
    )

    parser.add_argument(
        '--anno',
        action='store_true',
        help='显示注释配置选项|Show annotation configuration options'
    )

    # ========== 文件生成模式参数 | File Generation Mode Parameters ==========
    parser.add_argument(
        '--type',
        choices=['fa', 'prot', 'anno'],
        help='文件类型(用于生成文件列表)|File type for generating file list (fa/prot/anno)'
    )

    parser.add_argument(
        '--path',
        help='文件路径|File path for generating list'
    )

    parser.add_argument(
        '--out',
        help='输出文件名|Output file name'
    )

    # ========== 通用参数 | Common Parameters ==========
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (default: 12)'
    )

    parser.add_argument(
        '--perl-path',
        default='~/miniforge3/envs/genomesyn2/bin/perl',
        help=argparse.SUPPRESS  # 隐藏参数|Hidden parameter
    )

    parser.add_argument(
        '--genomesyn2-pl',
        default='~/miniforge3/envs/genomesyn2/bin/GenomeSyn2.pl',
        help=argparse.SUPPRESS  # 隐藏参数|Hidden parameter
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    # 转换参数为字典|Convert args to dictionary
    kwargs = vars(args)

    # 移除None值的参数|Remove None values
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # 特殊处理anno参数|Special handling for anno parameter
    if 'anno' in kwargs:
        if kwargs['anno']:
            # anno为True时，在原始main中需要特殊处理
            # 这里我们保持原样，让后续处理逻辑决定
            pass

    try:
        runner = GenomeSyn2Runner(**kwargs)
        exit_code = runner.run()
        sys.exit(exit_code)
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
