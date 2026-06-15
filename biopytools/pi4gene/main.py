"""
Pi4Gene主模块|Pi4Gene Main Module
基因分组核苷酸多样性(pi)计算分析器
Nucleotide diversity (pi) calculator per gene group
"""

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path

from ..common.paths import expand_path
from .config import Pi4GeneConfig
from .utils import (
    Pi4GeneLogger, parse_id_file, get_software_version
)
from .calculator import Pi4GeneCalculator


class Pi4GeneAnalyzer:
    """Pi4Gene分析器|Pi4Gene Analyzer"""

    def __init__(
        self,
        input_file: str,
        id_file: str,
        output_dir: str,
        threads: int = 12,
        mafft_path: str = None
    ):
        """初始化Pi4Gene分析器|Initialize Pi4Gene Analyzer"""
        config_kwargs = {
            'input_file': input_file,
            'id_file': id_file,
            'output_dir': output_dir,
            'threads': threads,
        }

        if mafft_path:
            config_kwargs['mafft_path'] = mafft_path

        self.config = Pi4GeneConfig(**config_kwargs)
        self.config.validate()

        self.logger_manager = Pi4GeneLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

    def run(self):
        """运行完整的pi4gene分析流程|Run the complete pi4gene pipeline"""
        start_time = datetime.now()
        self.logger.info("=" * 60)
        self.logger.info("Pi4Gene 基因分组核苷酸多样性计算|Pi4Gene Nucleotide Diversity Calculation")
        self.logger.info("=" * 60)

        self._print_config()

        # 解析分组ID文件|Parse group ID file
        group_to_ids = parse_id_file(self.config.id_file)
        self.logger.info(
            f"检测到{len(group_to_ids)}个分组|"
            f"Detected {len(group_to_ids)} groups"
        )

        total_seqs = sum(len(ids) for ids in group_to_ids.values())
        self.logger.info(
            f"共{total_seqs}个序列ID|Total {total_seqs} sequence IDs"
        )

        # 打印分组摘要|Print group summary
        self._print_group_summary(group_to_ids)

        # 记录软件版本|Record software version
        mafft_ver = get_software_version(self.config.mafft_path, self.logger)
        self.logger.info(f"MAFFT版本|MAFFT version: {mafft_ver}")

        # 初始化计算器|Initialize calculator
        calculator = Pi4GeneCalculator(self.config, self.logger)

        # 步骤1: 提取序列|Step 1: Extract sequences
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("步骤1: 按分组提取序列|Step 1: Extract sequences by group")
        self.logger.info("=" * 60)
        group_fasta_map = calculator.extract_sequences_by_group(group_to_ids)

        if not group_fasta_map:
            self.logger.error("没有有效分组，无法继续|No valid groups, cannot continue")
            sys.exit(1)

        # 步骤2: MAFFT比对|Step 2: MAFFT alignment
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("步骤2: MAFFT多序列比对|Step 2: MAFFT multiple sequence alignment")
        self.logger.info("=" * 60)
        aligned_map = calculator.run_mafft(group_fasta_map)

        if not aligned_map:
            self.logger.error("没有成功比对的分组，无法继续|No successfully aligned groups, cannot continue")
            sys.exit(1)

        # 步骤3: 计算pi|Step 3: Calculate pi
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("步骤3: 计算核苷酸多样性(pi)|Step 3: Calculate pi")
        self.logger.info("=" * 60)
        results = calculator.calculate_pi(aligned_map)

        # 写入结果|Write results
        output_file = calculator.write_results(results)

        # 生成版本信息文件|Generate version info file
        self._generate_versions_yml({'mafft': mafft_ver}, start_time)

        # 完成|Complete
        end_time = datetime.now()
        runtime = int((end_time - start_time).total_seconds())
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(
            f"分析完成|Analysis completed in {runtime}秒|seconds"
        )
        self.logger.info(f"结果文件|Results file: {output_file}")
        self.logger.info("=" * 60)

    def _print_config(self):
        """打印配置信息|Print configuration"""
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  序列文件|Input file: {self.config.input_file}")
        self.logger.info(f"  分组ID文件|ID file: {self.config.id_file}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  MAFFT路径|MAFFT path: {self.config.mafft_path}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")

    def _print_group_summary(self, group_to_ids: dict):
        """打印分组摘要|Print group summary"""
        self.logger.info("分组信息摘要|Group summary:")
        for group_name in sorted(group_to_ids.keys()):
            n = len(group_to_ids[group_name])
            self.logger.info(f"  {group_name}: {n}条序列|sequences")

    def _generate_versions_yml(self, versions: dict, start_time: datetime):
        """生成software_versions.yml文件|Generate software_versions.yml file"""
        try:
            import yaml

            end_time = datetime.now()
            runtime_seconds = int((end_time - start_time).total_seconds())

            info = {
                'pipeline': {
                    'name': 'biopytools pi4gene',
                    'version': '1.0.0'
                },
                'tools': {},
                'parameters': {
                    'threads': self.config.threads,
                },
                'execution': {
                    'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'runtime_seconds': runtime_seconds
                }
            }

            for tool_name, version in versions.items():
                if tool_name == 'mafft':
                    info['tools']['mafft'] = {
                        'version': version,
                        'path': self.config.mafft_path
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
        description="基因分组核苷酸多样性(pi)计算工具|Nucleotide diversity (pi) calculation tool per gene group"
    )

    parser.add_argument('-i', '--input', required=True,
                        help='序列FASTA文件路径|Input sequence FASTA file path')
    parser.add_argument('-d', '--id-file', required=True,
                        help='分组ID文件路径（第一列分组，第二列序列ID）|Group ID file path (col1: group, col2: seq_id)')
    parser.add_argument('-o', '--output-dir', default='./pi4gene_output',
                        help='输出目录|Output directory')

    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('--mafft-path', default=None,
                        help='MAFFT路径|MAFFT path')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        analyzer = Pi4GeneAnalyzer(
            input_file=args.input,
            id_file=args.id_file,
            output_dir=args.output_dir,
            threads=args.threads,
            mafft_path=args.mafft_path
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
