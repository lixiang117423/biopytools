"""
VCF按染色体合并主程序模块|VCF Merge by Chromosome Main Module
"""

import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

from .config import VCFMergerConfig
from .utils import (
    VCFMergerLogger,
    VCFFileFinder,
    VCFGroupByChromosome,
    BCFToolsChecker,
    VCFMerger as VCFFileMerger,
    VCFIndexer
)


class VCFMerger:
    """
    VCF按染色体合并主类|Main VCF Merger by Chromosome Class

    功能|Features:
    - 自动识别染色体编号并分组
    - 使用bcftools合并同染色体的VCF文件
    - 生成压缩的VCF文件(.gz)
    - 自动创建索引文件(.csi)
    - 详细的日志记录
    """

    def __init__(self, **kwargs):
        """
        初始化VCF合并器|Initialize VCF merger

        Args:
            **kwargs: 配置参数|Configuration parameters
                - input_dir: 输入VCF文件目录|Input directory containing VCF files
                - output_dir: 输出目录|Output directory
                - pattern: VCF文件名模式|VCF file name pattern (default: *.joint.vcf.gz)
                - threads: 线程数|Number of threads (default: 4)
                - create_index: 是否创建索引|Whether to create index files (default: True)
                - log_file: 日志文件路径|Log file path (optional)
                - log_level: 日志级别|Log level (default: INFO)
                - quiet: 静默模式|Quiet mode (default: False)
                - verbose: 详细输出级别|Verbose level
        """
        # 初始化配置|Initialize configuration
        self.config = VCFMergerConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_level = self.config.get_effective_log_level()
        self.logger_manager = VCFMergerLogger(
            self.config.log_file,
            log_level
        )
        self.logger = self.logger_manager.get_logger()

        # 存储结果|Store results
        self.chr_groups: Optional[Dict[str, List[Path]]] = None
        self.success_count: int = 0
        self.failed_chromosomes: List[str] = []

    def run(self):
        """
        运行完整的合并流程|Run complete merge pipeline

        Returns:
            bool: 是否成功完成|Whether completed successfully
        """
        start_time = time.time()

        try:
            # 输出程序信息|Output program information
            self._print_header()

            # 检查bcftools|Check bcftools
            if not BCFToolsChecker.check_bcftools(self.logger):
                sys.exit(1)

            # 步骤1: 查找VCF文件|Step 1: Find VCF files
            self._find_vcf_files()

            # 步骤2: 按染色体分组|Step 2: Group by chromosome
            self._group_by_chromosome()

            # 步骤3: 创建输出目录|Step 3: Create output directory
            self._create_output_directory()

            # 步骤4: 合并VCF文件|Step 4: Merge VCF files
            self._merge_vcf_files()

            # 输出总结|Output summary
            self._print_summary(start_time)

            # 返回结果|Return result
            return len(self.failed_chromosomes) == 0

        except KeyboardInterrupt:
            self.logger.warning("用户中断|Interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"错误|Error: {e}", exc_info=True)
            sys.exit(1)

    def _print_header(self):
        """输出程序头部信息|Print program header"""
        self.logger.info("=" * 60)
        self.logger.info("VCF文件按染色体合并工具|VCF Merge by Chromosome Tool")
        self.logger.info("Version: 1.0.0")
        self.logger.info("=" * 60)
        self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"文件模式|File pattern: {self.config.pattern}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"创建索引|Create index: {self.config.create_index}")
        self.logger.info("=" * 60)

    def _find_vcf_files(self):
        """步骤1: 查找VCF文件|Step 1: Find VCF files"""
        self.logger.info("=" * 60)
        self.logger.info("STEP 1: 查找VCF文件|Finding VCF files")
        self.logger.info("=" * 60)

        self.vcf_files = VCFFileFinder.find_vcf_files(
            self.config.input_dir,
            self.config.pattern,
            self.logger
        )

    def _group_by_chromosome(self):
        """步骤2: 按染色体分组|Step 2: Group by chromosome"""
        self.logger.info("=" * 60)
        self.logger.info("STEP 2: 按染色体分组|Grouping by chromosome")
        self.logger.info("=" * 60)

        self.chr_groups = VCFGroupByChromosome.group_vcf_by_chromosome(
            self.vcf_files,
            self.logger
        )

        if not self.chr_groups:
            self.logger.error("没有有效的VCF文件组|No valid VCF file groups")
            sys.exit(1)

        self.logger.info(
            f"识别出 {len(self.chr_groups)} 个染色体|"
            f"Identified {len(self.chr_groups)} chromosomes"
        )

        # 排序显示染色体|Sort and display chromosomes
        for chr_id in sorted(self.chr_groups.keys(), key=lambda x: (len(x), x)):
            count = len(self.chr_groups[chr_id])
            self.logger.info(f"  {chr_id}: {count} 个文件|files")

    def _create_output_directory(self):
        """步骤3: 创建输出目录|Step 3: Create output directory"""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(
            f"输出目录已创建|Output directory created: {self.config.output_dir}"
        )

    def _merge_vcf_files(self):
        """步骤4: 合并VCF文件|Step 4: Merge VCF files"""
        self.logger.info("=" * 60)
        self.logger.info("STEP 3: 合并VCF文件|Merging VCF files")
        self.logger.info("=" * 60)

        self.success_count = 0
        self.failed_chromosomes = []

        for chr_id in sorted(self.chr_groups.keys(), key=lambda x: (len(x), x)):
            vcf_list = self.chr_groups[chr_id]
            output_file = self.config.output_dir / f"{chr_id}.joint.merged.vcf.gz"

            self.logger.info(f"处理染色体 {chr_id}|Processing chromosome {chr_id}")

            # 合并VCF文件|Merge VCF files
            if VCFFileMerger.merge_vcf_files(
                vcf_list,
                output_file,
                self.config.threads,
                self.logger
            ):
                self.success_count += 1

                # 创建索引（如果启用）| Create index (if enabled)
                if self.config.create_index:
                    VCFIndexer.index_vcf_file(output_file, self.logger)
            else:
                self.failed_chromosomes.append(chr_id)

    def _print_summary(self, start_time: float):
        """输出总结信息|Print summary"""
        elapsed_time = time.time() - start_time

        self.logger.info("=" * 60)
        self.logger.info("合并完成|Merge completed")
        self.logger.info("=" * 60)
        self.logger.info(
            f"成功合并|Successfully merged: {self.success_count}/{len(self.chr_groups)} 个染色体"
        )
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"总耗时|Total runtime: {elapsed_time:.2f} 秒|seconds")

        if self.failed_chromosomes:
            self.logger.warning(
                f"失败的染色体|Failed chromosomes: {', '.join(self.failed_chromosomes)}"
            )
        else:
            self.logger.info("所有操作完成|All operations completed")

    # 公共方法|Public methods

    def get_chromosome_groups(self) -> Optional[Dict[str, List[Path]]]:
        """
        获取染色体分组信息|Get chromosome group information

        Returns:
            染色体分组字典|Chromosome group dictionary
        """
        return self.chr_groups

    def get_statistics(self) -> Dict[str, any]:
        """
        获取统计信息|Get statistics

        Returns:
            统计信息字典|Statistics dictionary
        """
        return {
            "total_chromosomes": len(self.chr_groups) if self.chr_groups else 0,
            "success_count": self.success_count,
            "failed_count": len(self.failed_chromosomes),
            "failed_chromosomes": self.failed_chromosomes,
            "input_files": len(self.vcf_files) if hasattr(self, "vcf_files") else 0
        }


def main():
    """
    主函数入口|Main function entry point

    此函数用于作为模块直接运行时的入口
    This function serves as the entry point when the module is run directly
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="按染色体合并VCF文件|Merge VCF files by chromosome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i /path/to/vcf_files -o /path/to/output
        '''
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input-dir",
        required=True,
        help="输入VCF文件目录|Input directory containing VCF files"
    )
    required.add_argument(
        "-o", "--output-dir",
        required=True,
        help="输出目录|Output directory"
    )

    # 可选参数|Optional arguments
    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="使用的线程数|Number of threads"
    )
    optional.add_argument(
        "--pattern",
        default="*.joint.vcf.gz",
        help="VCF文件名模式|VCF file name pattern"
    )
    optional.add_argument(
        "--no-index",
        action="store_true",
        help="不生成索引文件|Do not generate index files"
    )

    # 日志参数|Logging parameters
    logging_group = parser.add_argument_group("logging options")
    logging_group.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="详细输出模式|Verbose mode (-v: INFO, -vv: DEBUG)"
    )
    logging_group.add_argument(
        "--quiet",
        action="store_true",
        help="静默模式|Quiet mode (only ERROR)"
    )

    args = parser.parse_args()

    # 创建合并器并运行|Create merger and run
    merger = VCFMerger(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        pattern=args.pattern,
        threads=args.threads,
        create_index=not args.no_index,
        quiet=args.quiet,
        verbose=args.verbose
    )

    success = merger.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
