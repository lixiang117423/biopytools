"""
Fst计算主程序模块|Fst Calculation Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import FstConfig
from .utils import FstLogger
from .fst_calculator import FstCalculator
from .results_processor import FstResultsProcessor


class FstAnalyzer:
    """Fst分析主类|Main Fst Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = FstConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FstLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化Fst计算器|Initialize Fst calculator
        self.fst_calculator = FstCalculator(self.config, self.logger)

        # 初始化结果处理器|Initialize results processor
        self.results_processor = FstResultsProcessor(self.logger, self.config.output_path)

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            self.logger.info("开始Fst计算流程|Starting Fst calculation pipeline")

            # 运行Fst计算|Run Fst calculation
            output_files = self.fst_calculator.run_analysis()

            # 检查是否成功|Check if successful
            if not output_files:
                self.logger.error("Fst计算失败|Fst calculation failed")
                sys.exit(1)

            # 判断是bootstrap模式还是直接计算模式|Determine bootstrap or direct mode
            is_bootstrap = 'bootstrap_results' in output_files

            if is_bootstrap:
                # Bootstrap模式|Bootstrap mode
                self.logger.info("Bootstrap模式完成|Bootstrap mode completed")
                # 结果已经在fst_calculator中处理|Results already processed in fst_calculator
            else:
                # 直接计算模式|Direct calculation mode
                if 'fst' not in output_files:
                    self.logger.error("Fst计算失败|Fst calculation failed")
                    sys.exit(1)

                # 处理结果|Process results
                fst_file = output_files['fst']
                populations = self.fst_calculator.filtered_populations

                result_files = self.results_processor.process_results(fst_file, populations)

                # 合并输出文件列表|Merge output file lists
                output_files.update(result_files)

            # 打印输出文件列表|Print output file list
            self.logger.info("输出文件列表|Output files:")
            for file_type, file_path in output_files.items():
                self.logger.info(f"  {file_type}: {file_path}")

            self.logger.info("Fst计算流程完成|Fst calculation pipeline completed")

            return output_files

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            sys.exit(1)


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Fst计算工具 - 计算群体间遗传分化系数|Fst Calculation Tool - Calculate genetic differentiation coefficient between populations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i input.vcf -p pop.txt -o output_dir
  %(prog)s -i input.vcf -p pop.txt -o output_dir --maf 0.01 --geno 0.2
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--vcf-file', required=True,
                        help='VCF文件路径|VCF file path')
    parser.add_argument('-p', '--pop-file', required=True,
                        help='群体文件路径（样本ID + 群体标签）|Population file path (sample ID + population label)')
    parser.add_argument('-o', '--output-dir', default='./fst_output',
                        help='输出目录|Output directory (default: ./fst_output)')

    # 软件路径|Software path
    parser.add_argument('--plink-path', default=None,
                        help='PLINK软件路径|PLINK software path (default: auto-detect)')

    # 质控参数|Quality control parameters
    parser.add_argument('--enable-qc', action='store_true',
                        help='启用质控过滤（默认禁用）|Enable quality control filtering (disabled by default)')
    parser.add_argument('--maf', type=float, default=0.05,
                        help='最小等位基因频率阈值|Minor allele frequency threshold (default: 0.05)')
    parser.add_argument('--geno', type=float, default=0.1,
                        help='位点缺失率阈值|Genotype missing rate threshold (default: 0.1)')
    parser.add_argument('--mind', type=float, default=0.1,
                        help='样本缺失率阈值|Sample missing rate threshold (default: 0.1)')
    parser.add_argument('--hwe', type=float, default=1e-6,
                        help='Hardy-Weinberg平衡p值阈值|Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)')

    # 输出控制|Output control
    parser.add_argument('--no-keep-intermediate', action='store_true',
                        help='不保留中间文件|Do not keep intermediate files')

    # Bootstrap参数|Bootstrap parameters
    parser.add_argument('--enable-bootstrap', action='store_true',
                        help='启用bootstrap抽样|Enable bootstrap sampling')
    parser.add_argument('--bootstrap-iterations', type=int, default=100,
                        help='Bootstrap迭代次数|Bootstrap iterations (default: 100)')
    parser.add_argument('--min-samples', type=int, default=10,
                        help='最小样本数阈值，排除样本数少于此值的群体|Minimum sample count threshold, exclude populations with fewer samples (default: 10)')
    parser.add_argument('--exclude-pops', default=None,
                        help='手动指定要排除的群体（逗号分隔）|Manually specify populations to exclude (comma-separated)')

    # LD pruning参数|LD pruning parameters
    parser.add_argument('--no-ld-prune', action='store_true',
                        help='禁用LD pruning|Disable LD pruning')
    parser.add_argument('--ld-window', type=int, default=50,
                        help='LD pruning窗口大小|LD pruning window size (default: 50)')
    parser.add_argument('--ld-step', type=int, default=10,
                        help='LD pruning步长|LD pruning step size (default: 10)')
    parser.add_argument('--ld-r2', type=float, default=0.2,
                        help='LD pruning R2阈值|LD pruning R2 threshold (default: 0.2)')

    # SNP抽稀参数|SNP thinning parameters
    parser.add_argument('--thin', type=float, default=None,
                        help='SNP抽稀比例（0-1之间）|SNP thinning ratio (between 0-1)')

    # 并行参数|Parallel parameters
    parser.add_argument('--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 构建配置参数|Build configuration parameters
        config_params = {
            'vcf_file': args.vcf_file,
            'pop_file': args.pop_file,
            'output_dir': args.output_dir,
            'maf': args.maf,
            'geno': args.geno,
            'mind': args.mind,
            'hwe': args.hwe,
            'keep_intermediate': not args.no_keep_intermediate,
            # 质控控制|Quality control
            'enable_qc': args.enable_qc,
            # Bootstrap参数|Bootstrap parameters
            'enable_bootstrap': args.enable_bootstrap,
            'bootstrap_iterations': args.bootstrap_iterations,
            'min_samples': args.min_samples,
            'exclude_pops': args.exclude_pops,
            # LD pruning参数|LD pruning parameters
            'enable_ld_prune': not args.no_ld_prune,
            'ld_window_size': args.ld_window,
            'ld_step_size': args.ld_step,
            'ld_r2_threshold': args.ld_r2,
            # SNP抽稀参数|SNP thinning parameters
            'thin_threshold': args.thin,
            # 并行参数|Parallel parameters
            'threads': args.threads
        }

        # 添加可选参数|Add optional parameters
        if args.plink_path:
            config_params['plink_path'] = args.plink_path

        # 创建分析器并运行|Create analyzer and run
        analyzer = FstAnalyzer(**config_params)
        analyzer.run_analysis()

        sys.exit(0)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
