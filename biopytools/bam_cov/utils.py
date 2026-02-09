"""
BAM覆盖度统计工具函数模块|BAM Coverage Statistics Utility Functions Module
"""

import logging
import sys
import time
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple


class BAMCoverageLogger:
    """BAM覆盖度统计日志管理器|BAM Coverage Statistics Logger Manager"""

    def __init__(self, log_dir: Path, verbose: bool = False, quiet: bool = False):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.log_dir / f"bam_coverage_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"bam_coverage_{timestamp}")

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤|Log step"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


class SAMToolsHelper:
    """SAMTools辅助工具类|SAMTools Helper Class"""

    def __init__(self, logger, threads: int = 64):
        self.logger = logger
        self.threads = threads

    def check_samtools(self) -> bool:
        """检查samtools是否可用|Check if samtools is available"""
        try:
            result = subprocess.run(
                ['samtools', '--version'],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"samtools版本|samtools version: {result.stdout.split()[1]}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("samtools未安装或不在PATH中|samtools not installed or not in PATH")
            return False

    def get_chromosome_length(self, bam_file: str, chromosome: str) -> Optional[int]:
        """从BAM文件头部获取染色体长度|Get chromosome length from BAM header"""
        try:
            result = subprocess.run(
                ['samtools', 'view', '-H', bam_file],
                capture_output=True,
                text=True,
                check=True
            )

            for line in result.stdout.split('\n'):
                if line.startswith('@SQ') and f'SN:{chromosome}' in line:
                    # 提取LN字段|Extract LN field
                    for field in line.split('\t'):
                        if field.startswith('LN:'):
                            length = int(field.split(':')[1])
                            return length

            self.logger.warning(f"未找到染色体|Chromosome not found: {chromosome}")
            return None

        except subprocess.CalledProcessError as e:
            self.logger.error(f"读取BAM头部失败|Failed to read BAM header: {e}")
            return None

    def check_bam_index(self, bam_file: str) -> bool:
        """检查BAM索引是否存在|Check if BAM index exists"""
        index_file = f"{bam_file}.bai"
        if not Path(index_file).exists():
            # 尝试替代的索引文件名|Try alternative index filename
            bam_base = bam_file.replace('.bam', '')
            index_file = f"{bam_base}.bai"

        return Path(index_file).exists()

    def create_bam_index(self, bam_file: str) -> bool:
        """创建BAM索引|Create BAM index"""
        self.logger.info(f"创建BAM索引|Creating BAM index: {bam_file}")

        try:
            subprocess.run(
                ['samtools', 'index', '-@', str(self.threads), bam_file],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info("索引创建成功|Index created successfully")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"索引创建失败|Failed to create index: {e}")
            return False

    def extract_coverage(self, bam_file: str, chromosome: str, start: int,
                         end: Optional[int], min_mapq: int, min_baseq: int,
                         output_file: str) -> bool:
        """使用samtools depth提取覆盖度|Extract coverage using samtools depth"""
        region = f"{chromosome}:{start}-{end if end else ''}"

        cmd = [
            'samtools', 'depth',
            '-@', str(self.threads),
            '-r', region,
            '-Q', str(min_mapq),
            '-q', str(min_baseq),
            bam_file
        ]

        self.logger.debug(f"执行命令|Executing command: {' '.join(cmd)}")

        try:
            with open(output_file, 'w') as f:
                subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )

            # 检查输出文件|Check output file
            if Path(output_file).stat().st_size == 0:
                self.logger.warning(f"覆盖度文件为空|Coverage file is empty: {output_file}")
                return False

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"提取覆盖度失败|Failed to extract coverage: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class CoverageDataProcessor:
    """覆盖度数据处理类|Coverage Data Processor"""

    @staticmethod
    def read_coverage_file(file_path: str) -> Dict[Tuple[str, int], int]:
        """读取覆盖度文件|Read coverage file"""
        coverage_data = {}

        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    pos = int(parts[1])
                    depth = int(parts[2])
                    coverage_data[(chrom, pos)] = depth

        return coverage_data

    @staticmethod
    def merge_coverage_data(coverage_files: List[str], sample_names: List[str],
                           output_file: str, logger) -> bool:
        """合并多个样本的覆盖度数据|Merge coverage data from multiple samples"""
        logger.info(f"合并{len(coverage_files)}个样本的覆盖度数据|Merging coverage data from {len(coverage_files)} samples")

        try:
            # 步骤1: 获取所有唯一位置|Step 1: Get all unique positions
            logger.info("步骤1: 收集所有位置信息|Step 1: Collecting all positions")
            all_positions = set()

            coverage_data_list = []
            for file_path in coverage_files:
                data = CoverageDataProcessor.read_coverage_file(file_path)
                coverage_data_list.append(data)
                all_positions.update(data.keys())

            # 排序位置|Sort positions
            sorted_positions = sorted(all_positions, key=lambda x: (x[0], x[1]))
            logger.info(f"共有 {len(sorted_positions)} 个唯一位置|Total {len(sorted_positions)} unique positions")

            # 步骤2: 写入合并文件|Step 2: Write merged file
            logger.info("步骤2: 写入合并文件|Step 2: Writing merged file")
            with open(output_file, 'w') as out:
                # 写入表头|Write header
                header = ['Chrom', 'Pos'] + sample_names
                out.write('\t'.join(header) + '\n')

                # 写入数据|Write data
                for chrom, pos in sorted_positions:
                    row = [chrom, str(pos)]

                    for data in coverage_data_list:
                        depth = data.get((chrom, pos), 0)
                        row.append(str(depth))

                    out.write('\t'.join(row) + '\n')

            logger.info(f"合并完成，输出文件|Merged, output file: {output_file}")
            return True

        except Exception as e:
            logger.error(f"合并覆盖度数据失败|Failed to merge coverage data: {e}")
            return False

    @staticmethod
    def calculate_summary_statistics(coverage_file: str, output_file: str, logger) -> bool:
        """计算统计摘要|Calculate summary statistics"""
        logger.info("计算统计摘要|Calculating summary statistics")

        try:
            import statistics

            # 读取数据|Read data
            with open(coverage_file, 'r') as f:
                lines = f.readlines()

            if len(lines) <= 1:
                logger.warning("数据文件为空|Data file is empty")
                return False

            # 解析表头|Parse header
            header = lines[0].strip().split('\t')
            sample_names = header[2:]  # 跳过Chrom和Pos|Skip Chrom and Pos

            # 初始化统计数据|Initialize statistics
            stats = {}

            for sample_idx, sample_name in enumerate(sample_names, start=2):
                depths = []

                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) > sample_idx:
                        depth = int(parts[sample_idx])
                        depths.append(depth)

                if not depths:
                    continue

                # 计算统计量|Calculate statistics
                mean_depth = statistics.mean(depths)
                median_depth = statistics.median(depths)
                max_depth = max(depths)
                min_depth = min(depths)

                # 计算覆盖度比例|Calculate coverage ratios
                total_positions = len(depths)
                positions_gt0 = sum(1 for d in depths if d > 0)
                positions_gt10 = sum(1 for d in depths if d > 10)
                positions_gt30 = sum(1 for d in depths if d > 30)

                pct_gt0 = (positions_gt0 / total_positions) * 100 if total_positions > 0 else 0
                pct_gt10 = (positions_gt10 / total_positions) * 100 if total_positions > 0 else 0

                stats[sample_name] = {
                    'total_positions': total_positions,
                    'mean_coverage': mean_depth,
                    'median_coverage': median_depth,
                    'max_coverage': max_depth,
                    'min_coverage': min_depth,
                    'positions_0x': total_positions - positions_gt0,
                    'positions_gt0': positions_gt0,
                    'positions_gt10': positions_gt10,
                    'positions_gt30': positions_gt30,
                    'coverage_pct_gt0': pct_gt0,
                    'coverage_pct_gt10': pct_gt10
                }

            # 写入统计文件|Write statistics file
            with open(output_file, 'w') as out:
                # 表头|Header
                out.write('Sample\tTotal_Positions\tMean_Coverage\tMedian_Coverage\t')
                out.write('Max_Coverage\tMin_Coverage\tPositions_0X\tPositions_>0X\t')
                out.write('Positions_>10X\tPositions_>30X\tCoverage_%_>0X\tCoverage_%_>10X\n')

                # 数据行|Data rows
                for sample_name, stat in stats.items():
                    out.write(f"{sample_name}\t")
                    out.write(f"{stat['total_positions']}\t")
                    out.write(f"{stat['mean_coverage']:.2f}\t")
                    out.write(f"{stat['median_coverage']:.2f}\t")
                    out.write(f"{stat['max_coverage']}\t")
                    out.write(f"{stat['min_coverage']}\t")
                    out.write(f"{stat['positions_0x']}\t")
                    out.write(f"{stat['positions_gt0']}\t")
                    out.write(f"{stat['positions_gt10']}\t")
                    out.write(f"{stat['positions_gt30']}\t")
                    out.write(f"{stat['coverage_pct_gt0']:.2f}\t")
                    out.write(f"{stat['coverage_pct_gt10']:.2f}\n")

            logger.info(f"统计摘要已保存|Summary statistics saved: {output_file}")
            return True

        except Exception as e:
            logger.error(f"计算统计摘要失败|Failed to calculate summary: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
