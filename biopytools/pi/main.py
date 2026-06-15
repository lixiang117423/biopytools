"""
Pi计算主模块|Pi Calculation Main Module
核苷酸多样性(pi)计算分析器，使用vcftools计算
Nucleotide diversity (pi) calculator, using vcftools
"""

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path

from ..common.paths import expand_path
from .config import PiConfig
from .utils import PiLogger, parse_population_file, get_software_version, ensure_tabix_index, parse_fai_file
from .vcftools_pi import VcftoolsPiCalculator
from .results_merger import PiResultsMerger


class PiAnalyzer:
    """Pi计算分析器|Pi Calculation Analyzer"""

    def __init__(
        self,
        vcf_file: str,
        pop_file: str,
        output_dir: str,
        genome: str,
        window_size: int = None,
        window_step: int = None,
        default_window_size: int = 100000,
        default_window_step: int = 100000,
        maf: float = 0.0,
        max_missing: float = 1.0,
        threads: int = 12,
        keep_intermediate: bool = False,
        vcftools_path: str = None
    ):
        """初始化Pi分析器|Initialize Pi Analyzer"""
        # 解析genome参数：支持直接传入.fai或.fasta/.fa（自动找.fai）
        # Parse genome param: accept .fai or .fasta/.fa (auto-find .fai)
        genome_fai = self._resolve_fai_path(genome)

        # 构建配置|Build config
        config_kwargs = {
            'vcf_file': vcf_file,
            'pop_file': pop_file,
            'output_dir': output_dir,
            'genome_fai': genome_fai,
            'window_size': window_size,
            'window_step': window_step,
            'default_window_size': default_window_size,
            'default_window_step': default_window_step,
            'maf': maf,
            'max_missing': max_missing,
            'threads': threads,
            'keep_intermediate': keep_intermediate,
        }

        if vcftools_path:
            config_kwargs['vcftools_path'] = vcftools_path

        self.config = PiConfig(**config_kwargs)
        self.config.validate()

        # 初始化日志|Initialize logger
        self.logger_manager = PiLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

    @staticmethod
    def _resolve_fai_path(genome: str) -> str:
        """
        解析基因组路径，返回.fai文件路径
        Resolve genome path, return .fai file path

        Args:
            genome: 参考基因组fasta路径或.fai路径|Reference genome fasta path or .fai path

        Returns:
            .fai文件路径|.fai file path
        """
        genome = expand_path(genome)
        if genome.endswith('.fai'):
            return genome
        # 尝试同目录下的.fai文件|Try .fai file in same directory
        fai_path = genome + '.fai'
        if os.path.exists(fai_path):
            return fai_path
        raise FileNotFoundError(
            f"未找到.fai索引文件|Failed to find .fai index file: {fai_path}，"
            f"请先用samtools faidx生成|Please generate with samtools faidx first"
        )

    def run(self):
        """
        运行完整的pi计算分析流程
        Run the complete pi calculation pipeline
        """
        start_time = datetime.now()
        self.logger.info("=" * 60)
        self.logger.info("Pi核苷酸多样性计算分析|Pi Nucleotide Diversity Calculation")
        self.logger.info("=" * 60)

        # 打印配置信息|Print configuration
        self._print_config()

        # 读取染色体长度|Read chromosome lengths
        chrom_lengths = parse_fai_file(self.config.genome_fai)
        total_bases = sum(chrom_lengths.values())
        self.logger.info(
            f"参考基因组|Reference genome: {len(chrom_lengths)}条染色体|chromosomes, "
            f"总碱基数|total bases: {total_bases:,}"
        )

        # 检查并创建VCF tabix索引|Check and create VCF tabix index
        if not ensure_tabix_index(self.config.vcf_file, self.config.vcftools_path, self.logger):
            raise RuntimeError(
                f"VCF索引创建失败，无法继续|Failed to create VCF index, cannot continue: "
                f"{self.config.vcf_file}"
            )

        # 解析群体文件并打印摘要|Parse population file and print summary
        pop_to_samples = parse_population_file(self.config.pop_file)
        self._print_population_summary(pop_to_samples)

        # 记录软件版本|Record software versions
        versions = {}
        vcftools_ver = get_software_version(self.config.vcftools_path, self.logger)
        self.logger.info(f"vcftools版本|vcftools version: {vcftools_ver}")
        versions['vcftools'] = vcftools_ver

        # 运行计算|Run calculations
        vcftools_results = []
        windowed_vcftools_results = []

        is_windowed = self.config.window_size is not None

        if is_windowed:
            # 仅滑窗模式|Windowed-only mode
            calculator = VcftoolsPiCalculator(self.config, self.logger, chrom_lengths)
            vcftools_results = calculator.calculate()
        else:
            # 全基因组模式|Genome-wide mode
            calculator = VcftoolsPiCalculator(self.config, self.logger, chrom_lengths)
            vcftools_results = calculator.calculate()

            # 同步运行滑窗计算|Also run windowed calculation with default params
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("同步运行滑窗计算|Running windowed pi calculation alongside genome-wide")
            self.logger.info("=" * 60)
            self.logger.info(
                f"  窗口大小|Window size: {self.config.default_window_size}bp, "
                f"步长|Step: {self.config.default_window_step}bp"
            )

            from copy import deepcopy
            windowed_config = deepcopy(self.config)
            windowed_config.window_size = self.config.default_window_size
            windowed_config.window_step = self.config.default_window_step
            # 滑窗结果输出到独立子目录，避免与全基因组结果冲突
            # Output windowed results to separate subdirectory
            windowed_config.output_path = self.config.output_path / '03_windowed'
            windowed_config.output_path.mkdir(parents=True, exist_ok=True)
            (windowed_config.output_path / '01_vcftools').mkdir(parents=True, exist_ok=True)

            windowed_calc = VcftoolsPiCalculator(windowed_config, self.logger, chrom_lengths)
            windowed_vcftools_results = windowed_calc.calculate()

            # 合并滑窗结果|Merge windowed results
            windowed_merger = PiResultsMerger(windowed_config, self.logger)
            windowed_merger.merge_and_write(windowed_vcftools_results)

        # 合并全基因组结果|Merge genome-wide results
        merger = PiResultsMerger(self.config, self.logger)
        output_file = merger.merge_and_write(vcftools_results)

        # 生成版本信息文件|Generate version info file
        self._generate_versions_yml(versions, start_time)

        # 完成|Complete
        end_time = datetime.now()
        runtime = int((end_time - start_time).total_seconds())
        self.logger.info("=" * 60)
        self.logger.info(
            f"分析完成|Analysis completed in {runtime}秒|seconds"
        )
        self.logger.info(f"结果文件|Results file: {output_file}")
        self.logger.info("=" * 60)

    def _print_config(self):
        """打印配置信息|Print configuration"""
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  VCF文件|VCF file: {self.config.vcf_file}")
        self.logger.info(f"  群体文件|Population file: {self.config.pop_file}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  参考基因组fai|Reference genome fai: {self.config.genome_fai}")
        if self.config.window_size:
            self.logger.info(f"  窗口大小|Window size: {self.config.window_size}bp")
            step = self.config.window_step if self.config.window_step else self.config.window_size
            self.logger.info(f"  窗口步长|Window step: {step}bp")
        else:
            self.logger.info(f"  计算模式|Calculation mode: 全基因组|genome-wide")
        if self.config.maf > 0:
            self.logger.info(f"  MAF阈值|MAF threshold: {self.config.maf}")
        if self.config.max_missing < 1.0:
            self.logger.info(f"  最大缺失率|Max missing rate: {self.config.max_missing}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")

    def _print_population_summary(self, pop_to_samples: dict):
        """打印群体摘要信息|Print population summary"""
        self.logger.info("群体信息摘要|Population summary:")
        total = 0
        for pop_name in sorted(pop_to_samples.keys()):
            n = len(pop_to_samples[pop_name])
            total += n
            self.logger.info(f"  {pop_name}: {n}个样本|samples")
        self.logger.info(f"  总计|Total: {total}个样本|samples, {len(pop_to_samples)}个群体|populations")

    def _generate_versions_yml(self, versions: dict, start_time: datetime):
        """生成software_versions.yml文件|Generate software_versions.yml file"""
        try:
            import yaml

            end_time = datetime.now()
            runtime_seconds = int((end_time - start_time).total_seconds())

            info = {
                'pipeline': {
                    'name': 'biopytools pi',
                    'version': '2.0.0'
                },
                'tools': {},
                'parameters': {
                    'window_size': self.config.window_size,
                    'maf': self.config.maf,
                    'max_missing': self.config.max_missing,
                    'threads': self.config.threads,
                },
                'execution': {
                    'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'runtime_seconds': runtime_seconds
                }
            }

            for tool_name, version in versions.items():
                if tool_name == 'vcftools':
                    info['tools']['vcftools'] = {
                        'version': version,
                        'path': self.config.vcftools_path
                    }

            output_file = Path(self.config.output_dir) / '00_pipeline_info' / 'software_versions.yml'
            with open(output_file, 'w') as f:
                yaml.dump(info, f, default_flow_style=False, allow_unicode=True)

            self.logger.info(f"版本信息已保存|Version info saved: {output_file}")

        except ImportError:
            self.logger.warning("yaml模块未安装，跳过版本信息生成|yaml module not installed, skipping version info generation")


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="核苷酸多样性(pi)计算工具|Nucleotide diversity (pi) calculation tool"
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--vcf-file', required=True,
                        help='VCF文件路径（bgzip压缩+tabix索引）|VCF file path (bgzip-compressed + tabix-indexed)')
    parser.add_argument('-p', '--pop-file', required=True,
                        help='群体文件路径（样本ID 群体名）|Population file path (sample_id population_name)')
    parser.add_argument('-g', '--genome', required=True,
                        help='参考基因组fasta文件路径（需有.fai索引）|Reference genome fasta path (requires .fai index)')
    parser.add_argument('-o', '--output-dir', default='./pi_output',
                        help='输出目录|Output directory')

    # 窗口参数|Window parameters
    parser.add_argument('-w', '--window-size', type=int, default=None,
                        help='窗口大小bp（不设置则全基因组计算）|Window size in bp (omit for genome-wide)')
    parser.add_argument('--window-step', type=int, default=None,
                        help='窗口步长bp（默认等于窗口大小）|Window step in bp (default=window_size)')
    parser.add_argument('--default-window-size', type=int, default=100000,
                        help='自动滑窗默认窗口大小bp|Default windowed window size in bp (default: 100000)')
    parser.add_argument('--default-window-step', type=int, default=100000,
                        help='自动滑窗默认步长bp|Default windowed step in bp (default: 100000)')

    # 质控参数|QC parameters
    parser.add_argument('--maf', type=float, default=0.0,
                        help='最小等位基因频率|Minor allele frequency (default: 0.0)')
    parser.add_argument('--max-missing', type=float, default=1.0,
                        help='最大缺失率|Maximum missing rate (default: 1.0)')

    # 其他参数|Other parameters
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help='保留中间文件|Keep intermediate files')
    parser.add_argument('--vcftools-path', default=None,
                        help='vcftools路径|vcftools path')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        analyzer = PiAnalyzer(
            vcf_file=args.vcf_file,
            pop_file=args.pop_file,
            output_dir=args.output_dir,
            genome=args.genome,
            window_size=args.window_size,
            window_step=args.window_step,
            default_window_size=args.default_window_size,
            default_window_step=args.default_window_step,
            maf=args.maf,
            max_missing=args.max_missing,
            threads=args.threads,
            keep_intermediate=args.keep_intermediate,
            vcftools_path=args.vcftools_path
        )
        analyzer.run()
        sys.exit(0)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
