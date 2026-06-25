"""
FASTQ配对修复主程序模块|FASTQ Pair Fixing Main Module
"""

import os
import sys
from typing import Dict, Tuple
from .config import FastqPairConfig
from .utils import FastqPairLogger, CommandRunner, PairFinder


class FastqPairFixer:
    """FASTQ配对修复主类|Main FASTQ Pair Fixing Class"""

    def __init__(self, **kwargs):
        """
        初始化配对修复器|Initialize pair fixer

        参数支持|Parameters supported:
            input_dir: 输入目录|Input directory (required)
            output_dir: 输出目录|Output directory (required)
            suffix1: R1文件后缀|R1 file suffix (default: "_1.fq.gz")
            suffix2: R2文件后缀|R2 file suffix (default: "_2.fq.gz")
            threads: 线程数|Number of threads (default: 8)
            seqkit_bin: seqkit二进制路径|seqkit binary path (default: "seqkit")
            dry_run: 仅显示命令不执行|Dry run mode (default: False)
            verbose: 详细输出|Verbose output (default: False)
            log_file: 日志文件路径|Log file path (optional)
            log_level: 日志级别|Log level (default: "INFO")
        """
        # 初始化配置|Initialize configuration
        self.config = FastqPairConfig(**kwargs)

        # 验证配置|Validate configuration
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FastqPairLogger(
            self.config.output_path,
            self.config.log_file,
            self.config.log_level
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化配对查找器|Initialize pair finder
        self.pair_finder = PairFinder(
            self.config.input_dir,
            self.config.suffix1,
            self.config.suffix2,
            self.logger
        )

        # 处理结果统计|Processing result statistics
        self.total_samples = 0
        self.success_count = 0
        self.failed_samples = []

    def process_sample(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
        """
        处理单个样本|Process single sample

        Args:
            sample_name: 样本名称|Sample name
            r1_file: R1文件路径|R1 file path
            r2_file: R2文件路径|R2 file path

        Returns:
            是否成功|Success status
        """
        self.logger.info(f"处理样本: {sample_name}|Processing sample: {sample_name}")

        # 根据工具选择处理方法|Choose processing method based on tool
        if self.config.tool == "seqkit":
            success = self._process_with_seqkit(sample_name, r1_file, r2_file)
        elif self.config.tool == "repair":
            success = self._process_with_repair(sample_name, r1_file, r2_file)
        else:
            self.logger.error(f"未知工具|Unknown tool: {self.config.tool}")
            return False

        return success

    def _process_with_seqkit(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
        """使用seqkit处理样本|Process sample with seqkit"""
        # 创建子目录用于组织输出文件|Create subdirectories for organizing output files
        paired_dir = os.path.join(self.config.output_dir, "paired")
        single_dir = os.path.join(self.config.output_dir, "single")
        os.makedirs(paired_dir, exist_ok=True)
        os.makedirs(single_dir, exist_ok=True)

        # 构建seqkit命令|Build seqkit command
        cmd = (
            f"{self.config.seqkit_bin} pair "
            f"-1 \"{r1_file}\" "
            f"-2 \"{r2_file}\" "
            f"-O \"{paired_dir}\" "
            f"-j {self.config.threads}"
        )

        if self.config.dry_run:
            self.logger.info(f"[DRY RUN] {cmd}")
            return True

        # 执行命令|Execute command
        success = self.cmd_runner.run(
            cmd,
            f"配对修复样本 {sample_name} (seqkit)|Pair fixing sample {sample_name} (seqkit)"
        )

        if success:
            self.logger.info(f"  样本 {sample_name} 完成|Sample {sample_name} completed")

            # 检查并移动输出文件|Check and move output files
            output_r1 = os.path.join(paired_dir, f"{sample_name}{self.config.suffix1}")
            output_r2 = os.path.join(paired_dir, f"{sample_name}{self.config.suffix2}")
            single_file = os.path.join(paired_dir, f"{sample_name}.single.fq.gz")
            target_single = os.path.join(single_dir, f"{sample_name}_single.fq.gz")

            if os.path.exists(output_r1) and os.path.exists(output_r2):
                self.logger.debug(f"  输出R1|Output R1: {output_r1}")
                self.logger.debug(f"  输出R2|Output R2: {output_r2}")

            # 移动single文件到single目录|Move single file to single directory
            if os.path.exists(single_file):
                import shutil
                try:
                    shutil.move(single_file, target_single)
                    single_size = os.path.getsize(target_single)
                    single_size_mb = single_size / 1024 / 1024
                    if single_size_mb > 0:
                        self.logger.info(f"  单条reads|Singleton reads: {target_single} ({single_size_mb:.2f} MB)")
                except Exception as e:
                    self.logger.warning(f"  无法移动single文件|Failed to move single file: {e}")
            else:
                self.logger.debug(f"  无单条reads|No singleton reads for {sample_name}")
        else:
            self.logger.error(f"  样本 {sample_name} 失败|Sample {sample_name} failed")

        return success

    def _process_with_repair(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
        """使用repair.sh处理样本|Process sample with repair.sh"""
        # 创建子目录用于组织输出文件|Create subdirectories for organizing output files
        paired_dir = os.path.join(self.config.output_dir, "paired")
        single_dir = os.path.join(self.config.output_dir, "single")
        os.makedirs(paired_dir, exist_ok=True)
        os.makedirs(single_dir, exist_ok=True)

        # 定义输出文件|Define output files
        # 输出文件名和输入文件名保持一致|Output file names match input file names
        output_r1 = os.path.join(paired_dir, f"{sample_name}{self.config.suffix1}")
        output_r2 = os.path.join(paired_dir, f"{sample_name}{self.config.suffix2}")
        output_single = os.path.join(single_dir, f"{sample_name}_single.fq.gz")

        # 构建repair.sh命令|Build repair.sh command
        # 使用conda run调用repair.sh，确保在正确的conda环境中执行
        # 输出文件名保持和输入文件名一致，方便后续流程处理
        cmd = [
            "conda", "run", "-n", self.config.repair_conda_env,
            "--no-capture-output", self.config.repair_sh,
            f"in={r1_file}",
            f"in2={r2_file}",
            f"out={output_r1}",
            f"out2={output_r2}",
            f"outs={output_single}",
            "repair=t",
            f"threads={self.config.threads}",
            f"-Xmx{self.config.repair_memory}"
        ]

        if self.config.dry_run:
            self.logger.info(f"[DRY RUN] {' '.join(cmd)}")
            return True

        # 执行命令|Execute command
        success = self.cmd_runner.run(
            cmd,
            f"配对修复样本 {sample_name} (repair.sh)|Pair fixing sample {sample_name} (repair.sh)"
        )

        if success:
            self.logger.info(f"  样本 {sample_name} 完成|Sample {sample_name} completed")

            # 检查输出文件|Check output files
            if os.path.exists(output_r1) and os.path.exists(output_r2):
                r1_size = os.path.getsize(output_r1)
                r2_size = os.path.getsize(output_r2)
                self.logger.debug(f"  输出R1|Output R1: {output_r1} ({r1_size/1024/1024/1024:.2f} GB)")
                self.logger.debug(f"  输出R2|Output R2: {output_r2} ({r2_size/1024/1024/1024:.2f} GB)")

            if os.path.exists(output_single):
                single_size = os.path.getsize(output_single)
                single_size_mb = single_size / 1024 / 1024
                if single_size_mb > 0:
                    self.logger.info(f"  单条reads|Singleton reads: {output_single} ({single_size_mb:.2f} MB)")
        else:
            self.logger.error(f"  样本 {sample_name} 失败|Sample {sample_name} failed")

        return success

    def process_all_samples(self):
        """处理所有样本|Process all samples"""
        self.logger.info("=" * 60)
        self.logger.info("开始FASTQ配对修复流程|Starting FASTQ pair fixing pipeline")
        self.logger.info("=" * 60)

        # 显示配置信息|Show configuration
        self._show_config()

        # 查找配对文件|Find paired files
        pairs = self.pair_finder.find_pairs()

        if not pairs:
            self.logger.error("未找到任何配对文件|No paired files found")
            return False

        # 处理每个样本|Process each sample
        self.total_samples = len(pairs)
        self.logger.info(f"开始处理 {self.total_samples} 个样本|Start processing {self.total_samples} samples")
        self.logger.info("-" * 60)

        for idx, (sample_name, (r1_file, r2_file)) in enumerate(pairs.items(), 1):
            self.logger.info(f"[{idx}/{self.total_samples}] 处理中...|Processing...")

            success = self.process_sample(sample_name, r1_file, r2_file)

            if success:
                self.success_count += 1
            else:
                self.failed_samples.append(sample_name)

        # 处理结果汇总|Processing summary
        self._show_summary()

        return self.success_count == self.total_samples

    def _show_config(self):
        """显示配置信息|Show configuration"""
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  配对文件|Paired files: {os.path.join(self.config.output_dir, 'paired')}/")
        self.logger.info(f"  单条reads|Singleton reads: {os.path.join(self.config.output_dir, 'single')}/")
        self.logger.info(f"  R1后缀|R1 suffix: {self.config.suffix1}")
        self.logger.info(f"  R2后缀|R2 suffix: {self.config.suffix2}")
        self.logger.info(f"  工具|Tool: {self.config.tool}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        if self.config.tool == "seqkit":
            self.logger.info(f"  Seqkit路径|Seqkit path: {self.config.seqkit_bin}")
        elif self.config.tool == "repair":
            self.logger.info(f"  Repair.sh路径|Repair.sh path: {self.config.repair_sh}")
            self.logger.info(f"  Repair内存|Repair memory: {self.config.repair_memory}")
        self.logger.info(f"  Dry run模式|Dry run mode: {self.config.dry_run}")
        self.logger.info("-" * 60)

    def _show_summary(self):
        """显示处理结果汇总|Show processing summary"""
        self.logger.info("=" * 60)
        self.logger.info("处理结果汇总|Processing Summary:")
        self.logger.info(f"  总样本数|Total samples: {self.total_samples}")
        self.logger.info(f"  成功|Success: {self.success_count}/{self.total_samples}")
        self.logger.info(f"  失败|Failed: {len(self.failed_samples)}/{self.total_samples}")

        if self.failed_samples:
            self.logger.warning(f"失败样本列表|Failed samples: {', '.join(self.failed_samples)}")

        self.logger.info("=" * 60)


def main():
    """命令行入口函数|Command line entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='FASTQ配对修复工具|FASTQ Pair Fixing Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法|Basic usage
  python -m biopytools.pair_fastq.main -i raw_data -o fixed_data

  # 指定线程数|Specify threads
  python -m biopytools.pair_fastq.main -i raw_data -o fixed_data -t 16

  # 自定义文件后缀|Custom suffix
  python -m biopytools.pair_fastq.main -i raw_data -o fixed_data --suffix1 ".R1.fastq.gz"

  # Dry run模式|Dry run mode
  python -m biopytools.pair_fastq.main -i raw_data -o fixed_data --dry-run
        '''
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入目录（包含FASTQ文件）|Input directory containing FASTQ files'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出目录|Output directory'
    )

    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (default: 12)'
    )

    parser.add_argument(
        '--suffix1',
        default='_1.fq.gz',
        help='R1文件后缀|R1 file suffix (default: "_1.fq.gz")'
    )

    parser.add_argument(
        '--suffix2',
        default='_2.fq.gz',
        help='R2文件后缀|R2 file suffix (default: "_2.fq.gz")'
    )

    parser.add_argument(
        '--tool',
        type=str,
        default='repair',
        choices=['seqkit', 'repair'],
        help='工具选择|Tool selection (default: repair)'
    )

    parser.add_argument(
        '--seqkit-bin',
        default='seqkit',
        help='seqkit二进制路径|seqkit binary path (default: "seqkit")'
    )

    parser.add_argument(
        '--repair-sh',
        default='repair.sh',
        help='repair.sh脚本名称|repair.sh script name (default: "repair.sh")'
    )

    parser.add_argument(
        '--repair-conda-env',
        default='bbmap_v.39.81',
        help='repair.sh的conda环境名称|conda environment name for repair.sh (default: "bbmap_v.39.81")'
    )

    parser.add_argument(
        '--repair-memory',
        default='300g',
        help='repair.sh内存参数|repair.sh memory parameter (default: "300g")'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='仅显示命令不执行|Show commands without executing'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='详细输出|Verbose output'
    )

    parser.add_argument(
        '--log-file',
        help='日志文件路径|Log file path'
    )

    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='日志级别|Log level (default: INFO)'
    )

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )

    args = parser.parse_args()

    # 创建并运行配对修复器|Create and run pair fixer
    kwargs = {
        'input_dir': args.input,
        'output_dir': args.output,
        'threads': args.threads,
        'suffix1': args.suffix1,
        'suffix2': args.suffix2,
        'seqkit_bin': args.seqkit_bin,
        'dry_run': args.dry_run,
        'verbose': args.verbose,
        'log_file': args.log_file,
        'log_level': args.log_level
    }

    fixer = FastqPairFixer(**kwargs)
    success = fixer.process_all_samples()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
