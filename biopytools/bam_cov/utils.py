"""
BAM覆盖度统计工具函数模块|BAM Coverage Statistics Utility Functions Module
"""

import logging
import sys
import time
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Tuple


@dataclass
class BEDInterval:
    """BED区间数据类|BED Interval Data Class"""
    interval_id: int  # 行号编号(1-based)|Row number (1-based)
    chrom: str  # 染色体名称|Chromosome name
    start: int  # BED起始位置(0-based)|BED start (0-based)
    end: int  # BED终止位置(0-based, half-open)|BED end (0-based, half-open)
    name: Optional[str] = None  # 可选名称列|Optional name column


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
    def merge_bed_coverage_data(coverage_files: List[str], sample_names: List[str],
                                output_file: str, logger) -> bool:
        """合并BED模式多个样本的覆盖度数据|Merge BED mode coverage data from multiple samples

        输入文件格式: bed_start\\tbed_end\\tchrom\\tpos\\tdepth
        输出文件格式: Chrom\\tStart\\tEnd\\tPos\\tsample1\\tsample2...
        """
        logger.info(f"合并{len(coverage_files)}个样本的覆盖度数据|Merging coverage data from {len(coverage_files)} samples")

        try:
            all_keys = set()
            coverage_data_list = []

            for file_path in coverage_files:
                data = {}
                with open(file_path, 'r') as f:
                    for line in f:
                        if not line.strip():
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            key = (parts[2], int(parts[0]), int(parts[1]), int(parts[3]))
                            data[key] = int(parts[4])
                coverage_data_list.append(data)
                all_keys.update(data.keys())

            sorted_keys = sorted(all_keys, key=lambda x: (x[0], x[1], x[3]))
            logger.info(f"共有 {len(sorted_keys)} 个唯一位置|Total {len(sorted_keys)} unique positions")

            with open(output_file, 'w') as out:
                header = ['Chrom', 'Start', 'End', 'Pos'] + sample_names
                out.write('\t'.join(header) + '\n')

                for chrom, start, end, pos in sorted_keys:
                    row = [chrom, str(start), str(end), str(pos)]
                    for data in coverage_data_list:
                        depth = data.get((chrom, start, end, pos), 0)
                        row.append(str(depth))
                    out.write('\t'.join(row) + '\n')

            logger.info(f"合并完成，输出文件|Merged, output file: {output_file}")
            return True

        except Exception as e:
            logger.error(f"合并覆盖度数据失败|Failed to merge coverage data: {e}")
            return False

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

    @staticmethod
    def parse_bed_file(bed_path: str, logger) -> List[BEDInterval]:
        """解析BED文件|Parse BED file

        Args:
            bed_path: BED文件路径|BED file path
            logger: 日志记录器|Logger instance

        Returns:
            BEDInterval列表|List of BEDInterval objects
        """
        intervals = []
        skipped = 0
        interval_id = 0

        with open(bed_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    logger.warning(f"BED文件第{line_num}行格式错误，已跳过|Invalid BED line {line_num}, skipped")
                    skipped += 1
                    continue

                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    logger.warning(f"BED文件第{line_num}行坐标格式错误，已跳过|Invalid coordinate format at BED line {line_num}, skipped")
                    skipped += 1
                    continue

                if start < 0 or end <= start:
                    logger.warning(f"BED文件第{line_num}行区间无效(start={start}, end={end})，已跳过|Invalid interval at BED line {line_num} (start={start}, end={end}), skipped")
                    skipped += 1
                    continue

                interval_id += 1
                name = parts[3] if len(parts) >= 4 else None
                intervals.append(BEDInterval(
                    interval_id=interval_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name
                ))

        logger.info(f"BED文件解析完成|BED file parsed: {len(intervals)} 个有效区间|valid intervals")
        if skipped > 0:
            logger.warning(f"跳过|Skipped: {skipped} 行|lines")

        return intervals

    @staticmethod
    def calculate_interval_summary(depth_values: List[int]) -> Dict[str, float]:
        """计算单个区间单个样本的覆盖度统计|Calculate coverage statistics for a single interval of a single sample

        Args:
            depth_values: 所有位点的覆盖度值列表|List of coverage depth values for all positions

        Returns:
            统计指标字典|Dictionary of statistics metrics
        """
        import statistics

        if not depth_values:
            return {
                'mean_coverage': 0.0,
                'median_coverage': 0.0,
                'pct_gt0': 0.0,
                'pct_gt10': 0.0,
                'pct_gt30': 0.0,
                'total_positions': 0
            }

        total = len(depth_values)
        return {
            'mean_coverage': round(statistics.mean(depth_values), 2),
            'median_coverage': round(statistics.median(depth_values), 2),
            'pct_gt0': round((sum(1 for d in depth_values if d > 0) / total) * 100, 2),
            'pct_gt10': round((sum(1 for d in depth_values if d > 10) / total) * 100, 2),
            'pct_gt30': round((sum(1 for d in depth_values if d > 30) / total) * 100, 2),
            'total_positions': total
        }

    @staticmethod
    def read_depth_values(depth_file: str) -> List[int]:
        """从samtools depth输出文件读取覆盖度值列表|Read coverage depth values from samtools depth output file

        Args:
            depth_file: depth输出文件路径|Depth output file path

        Returns:
            覆盖度值列表|List of depth values
        """
        depth_values = []
        with open(depth_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    depth_values.append(int(parts[2]))
        return depth_values

    @staticmethod
    def write_bed_batch_summary(
        intervals: List[BEDInterval],
        sample_names: List[str],
        results: Dict[int, Dict[str, Dict[str, float]]],
        output_file: str,
        logger
    ) -> bool:
        """写入BED批量覆盖度统计摘要|Write BED batch coverage summary TSV

        Args:
            intervals: BED区间列表|List of BED intervals
            sample_names: 样本名列表|List of sample names
            results: 统计结果 {interval_id: {sample_name: {stat_name: value}}}|Statistics results
            output_file: 输出文件路径|Output file path
            logger: 日志记录器|Logger instance

        Returns:
            是否成功|Whether successful
        """
        try:
            # 构建表头|Build header
            header = ['Interval_ID', 'Chrom', 'Start', 'End', 'Length']
            for sample_name in sample_names:
                header.extend([
                    f'{sample_name}_Mean',
                    f'{sample_name}_%>0X',
                    f'{sample_name}_%>10X',
                    f'{sample_name}_%>30X'
                ])

            with open(output_file, 'w') as out:
                out.write('\t'.join(header) + '\n')

                for interval in intervals:
                    row = [
                        str(interval.interval_id),
                        interval.chrom,
                        str(interval.start),
                        str(interval.end),
                        str(interval.end - interval.start)
                    ]

                    for sample_name in sample_names:
                        stats = results.get(interval.interval_id, {}).get(sample_name, {})
                        row.append(str(stats.get('mean_coverage', 0.0)))
                        row.append(str(stats.get('pct_gt0', 0.0)))
                        row.append(str(stats.get('pct_gt10', 0.0)))
                        row.append(str(stats.get('pct_gt30', 0.0)))

                    out.write('\t'.join(row) + '\n')

            logger.info(f"BED批量汇总已保存|BED batch summary saved: {output_file}")
            return True

        except Exception as e:
            logger.error(f"写入BED批量汇总失败|Failed to write BED batch summary: {e}")
            return False
