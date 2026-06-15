"""
Kraken2主程序模块|Kraken2 Main Module
"""

import argparse
import sys
import os
from datetime import datetime

from .config import Kraken2Config
from .utils import Kraken2Logger, build_conda_command, CommandRunner
from .analyzer import Kraken2Analyzer


class Kraken2Pipeline:
    """Kraken2分析流程|Kraken2 Analysis Pipeline"""

    def __init__(self, **kwargs):
        self.config = Kraken2Config(**kwargs)
        self.config.validate()

        log_file = os.path.join(self.config.output_dir, '99_logs', 'kraken2_pipeline.log')
        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        self.logger_manager = Kraken2Logger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

        self.analyzer = Kraken2Analyzer(self.config, self.logger)

    def _save_pipeline_info(self):
        """保存流程信息|Save pipeline info"""
        info_dir = os.path.join(self.config.output_dir, '00_pipeline_info')
        os.makedirs(info_dir, exist_ok=True)

        info_file = os.path.join(info_dir, 'software_versions.yml')
        cmd_runner = CommandRunner(self.logger)

        # 获取Kraken2版本|Get Kraken2 version
        version = 'unknown'
        cmd = build_conda_command('kraken2', ['--version'])
        success, stdout, stderr = cmd_runner.run(cmd, "Kraken2版本检测|Kraken2 version check")
        if success and stdout:
            version = stdout.strip().split('\n')[0]

        info = [
            f"pipeline:",
            f"  name: biopytools kraken2",
            f"  version: 1.0.0",
            f"",
            f"tools:",
            f"  kraken2:",
            f"    version: \"{version}\"",
            f"",
            f"parameters:",
            f"  db_path: {self.config.db_path}",
            f"  threads: {self.config.threads}",
            f"  confidence: {self.config.confidence}",
            f"  bracken_level: {self.config.bracken_level}",
            f"  read_len: {self.config.read_len}",
            f"",
            f"execution:",
            f"  start_time: \"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\"",
        ]

        with open(info_file, 'w') as f:
            f.write('\n'.join(info) + '\n')

    def run(self):
        """运行流程|Run pipeline"""
        self.logger.info("Kraken2宏基因组分类流程启动|Kraken2 metagenomic classification pipeline started")
        self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"数据库路径|Database path: {self.config.db_path}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"Bracken: {'启用|enabled' if self.config.run_bracken else '禁用|disabled'}")

        self._save_pipeline_info()

        success = self.analyzer.run_pipeline()

        if success:
            self.logger.info("Kraken2分析流程全部完成|Kraken2 analysis pipeline completed successfully")
        else:
            self.logger.warning("Kraken2分析流程部分完成(存在失败样本)|Kraken2 analysis pipeline partially completed (some samples failed)")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Kraken2宏基因组分类工具|Kraken2 Metagenomic Classification Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input-dir', required=True,
                        help='输入FASTQ目录|Input FASTQ directory')
    parser.add_argument('-d', '--db-path', required=True,
                        help='Kraken2数据库路径|Kraken2 database path')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir', default='./kraken2_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads')
    parser.add_argument('--read-len', type=int, default=150,
                        help='读长(用于Bracken)|Read length (for Bracken)')
    parser.add_argument('--confidence', type=float, default=0.0,
                        help='Kraken2置信度阈值|Kraken2 confidence score threshold')
    parser.add_argument('--bracken-level', default='S',
                        choices=['D', 'P', 'C', 'O', 'F', 'G', 'S', 'S1'],
                        help='Bracken分类级别|Bracken taxonomic level')
    parser.add_argument('--bracken-threshold', type=int, default=10,
                        help='Bracken最小读数阈值|Bracken minimum read threshold')
    parser.add_argument('--no-bracken', action='store_true',
                        help='跳过Bracken分析|Skip Bracken analysis')
    parser.add_argument('--r1-suffix', default='_1.clean.fq.gz',
                        help='R1文件后缀|R1 file suffix')
    parser.add_argument('--r2-suffix', default='_2.clean.fq.gz',
                        help='R2文件后缀|R2 file suffix')

    args = parser.parse_args()

    pipeline = Kraken2Pipeline(
        input_dir=args.input_dir,
        db_path=args.db_path,
        output_dir=args.output_dir,
        threads=args.threads,
        read_len=args.read_len,
        confidence=args.confidence,
        bracken_level=args.bracken_level,
        bracken_threshold=args.bracken_threshold,
        run_bracken=not args.no_bracken,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
    )
    pipeline.run()


if __name__ == '__main__':
    main()
