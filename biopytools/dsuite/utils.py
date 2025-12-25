"""
Dsuite Utilities
Dsuite工具类
"""

import logging
import sys
import subprocess
import os
from datetime import datetime
from typing import Optional, Tuple, Dict


class DsuiteLogger:
    """Dsuite日志管理器 | Dsuite Logger Manager"""

    def __init__(self, output_dir: str):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            output_dir: 输出目录 | Output directory
        """
        self.output_dir = output_dir

        # 创建日志目录 | Create log directory
        os.makedirs(output_dir, exist_ok=True)

        # 设置日志文件名 | Set log file name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(output_dir, f"dsuite_analysis_{timestamp}.log")

        # 配置日志处理器 | Configure log handlers
        self._setup_logging(log_file)

    def _setup_logging(self, log_file: str):
        """
        设置日志 | Setup logging

        Args:
            log_file: 日志文件路径 | Log file path
        """
        # 创建logger | Create logger
        self.logger = logging.getLogger("Dsuite")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 创建formatter | Create formatter
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout处理器 - 只输出INFO及以上 | stdout handler - INFO and above only
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr处理器 - 输出WARNING及以上 | stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 文件处理器 - 输出所有级别 | File handler - all levels
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # 添加处理器 | Add handlers
        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class VCFStatsCollector:
    """VCF统计信息收集器 | VCF Statistics Collector"""

    def __init__(self, logger):
        """
        初始化统计器 | Initialize collector

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def collect_statistics(
        self,
        vcf_file: str,
        bcftools: str
    ) -> Dict[str, int]:
        """
        收集VCF文件统计信息 | Collect VCF file statistics

        Args:
            vcf_file: VCF文件路径 | VCF file path
            bcftools: bcftools命令 | bcftools command

        Returns:
            统计信息字典 | Statistics dictionary
        """
        self.logger.info("=" * 60)
        self.logger.info("VCF文件统计")
        self.logger.info("=" * 60)

        stats = {}

        try:
            # 统计样本数 | Count samples
            result = subprocess.run(
                [bcftools, 'query', '-l', vcf_file],
                capture_output=True,
                text=True,
                check=True
            )
            sample_count = len([line for line in result.stdout.strip().split('\n') if line])
            stats['sample_count'] = sample_count
            self.logger.info(f"样本数量: {sample_count}")

            # 统计总变异数 | Count total variants
            result = subprocess.run(
                f"{bcftools} view -H {vcf_file} | wc -l",
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            total_variants = int(result.stdout.strip())
            stats['total_variants'] = total_variants
            self.logger.info(f"总变异数: {total_variants}")

            # 统计染色体/scaffold数量 | Count chromosomes/scaffolds
            result = subprocess.run(
                f"{bcftools} view -h {vcf_file} | grep '^##contig' | wc -l",
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            chrom_count = int(result.stdout.strip())
            stats['chrom_count'] = chrom_count
            self.logger.info(f"染色体/Scaffold数: {chrom_count}")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"统计失败: {e.stderr}")
        except Exception as e:
            self.logger.error(f"统计时出错: {str(e)}")

        return stats

    def count_filtered_variants(
        self,
        vcf_file: str,
        bcftools: str,
        min_alleles: int,
        max_alleles: int,
        variant_type: str,
        verbose: bool = False
    ) -> int:
        """
        统计过滤后的变异数 | Count filtered variants

        Args:
            vcf_file: VCF文件路径 | VCF file path
            bcftools: bcftools命令 | bcftools command
            min_alleles: 最小等位基因数 | Min alleles
            max_alleles: 最大等位基因数 | Max alleles
            variant_type: 变异类型 | Variant type
            verbose: 是否输出详细日志 | Whether to output verbose logs

        Returns:
            过滤后的变异数 | Number of filtered variants
        """
        if verbose:
            self.logger.info("=" * 60)
            self.logger.info("应用过滤条件并统计")
            self.logger.info("=" * 60)
            self.logger.info(f"过滤参数:")
            self.logger.info(f"  - 最小等位基因数: {min_alleles}")
            self.logger.info(f"  - 最大等位基因数: {max_alleles}")
            self.logger.info(f"  - 变异类型: {variant_type}")

        try:
            # 使用bcftools过滤并统计 | Filter and count with bcftools
            cmd = f"{bcftools} view -m{min_alleles} -M{max_alleles} -v {variant_type} {vcf_file} | grep -v '^#' | wc -l"
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )

            count = int(result.stdout.strip())
            if verbose:
                self.logger.info(f"符合条件的变异数: {count}")

            return count

        except subprocess.CalledProcessError as e:
            self.logger.error(f"过滤统计失败: {e.stderr}")
            return 0


class DsuiteRunner:
    """Dsuite运行器 | Dsuite Runner"""

    def __init__(self, logger):
        """
        初始化运行器 | Initialize runner

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def run_dtrios(
        self,
        vcf_file: str,
        sets_file: str,
        output_prefix: str,
        numlines: int,
        dsuite_bin: str,
        bcftools: str,
        min_alleles: int = 2,
        max_alleles: int = 2,
        variant_type: str = "snps"
    ) -> bool:
        """
        运行Dsuite Dtrios | Run Dsuite Dtrios

        Args:
            vcf_file: VCF文件路径 | VCF file path
            sets_file: SETS文件路径 | SETS file path
            output_prefix: 输出前缀 | Output prefix
            numlines: VCF行数 | Number of VCF lines
            dsuite_bin: Dsuite可执行文件 | Dsuite binary
            bcftools: bcftools命令 | bcftools command
            min_alleles: 最小等位基因数 | Min alleles
            max_alleles: 最大等位基因数 | Max alleles
            variant_type: 变异类型 | Variant type

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("运行Dsuite Dtrios")
        self.logger.info("=" * 60)

        try:
            # 构建bcftools过滤命令 | Build bcftools filter command
            filter_cmd = [
                bcftools, 'view',
                f'-m{min_alleles}',
                f'-M{max_alleles}',
                '-v', variant_type,
                vcf_file
            ]

            # 构建Dsuite命令 | Build Dsuite command
            dsuite_cmd = [
                dsuite_bin, 'Dtrios',
                '-l', str(numlines),
                '-o', output_prefix,
                'stdin',
                sets_file
            ]

            self.logger.info(f"命令预览:")
            self.logger.info(f"  {' '.join(filter_cmd)} |")
            self.logger.info(f"  {' '.join(dsuite_cmd)}")

            # 执行管道命令 | Execute pipeline
            self.logger.info(f"开始分析...")

            # 启动bcftools进程 | Start bcftools process
            bcftools_proc = subprocess.Popen(
                filter_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            # 启动Dsuite进程 | Start Dsuite process
            dsuite_proc = subprocess.Popen(
                dsuite_cmd,
                stdin=bcftools_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # 关闭bcftools的stdout，让Dsuite能接收EOF
            bcftools_proc.stdout.close()

            # 等待Dsuite完成 | Wait for Dsuite to complete
            stdout, stderr = dsuite_proc.communicate()

            # 检查退出码 | Check exit codes
            bcftools_exit = bcftools_proc.wait()
            dsuite_exit = dsuite_proc.returncode

            if dsuite_exit != 0:
                self.logger.error(f"Dsuite运行失败 (退出码: {dsuite_exit})")
                self.logger.error(f"错误信息: {stderr}")
                return False

            if bcftools_exit != 0:
                self.logger.warning(f"bcftools警告 (退出码: {bcftools_exit})")

            self.logger.info("Dsuite运行成功")
            return True

        except Exception as e:
            self.logger.error(f"运行时出错: {str(e)}", exc_info=True)
            return False

    def check_output_files(self, output_prefix: str) -> Dict[str, bool]:
        """
        检查输出文件 | Check output files

        Args:
            output_prefix: 输出前缀 | Output prefix

        Returns:
            输出文件状态字典 | Output files status dictionary
        """
        self.logger.info("=" * 60)
        self.logger.info("结果验证")
        self.logger.info("=" * 60)

        expected_files = [
            f"{output_prefix}_BBAA.txt",
            f"{output_prefix}_Dmin.txt",
            f"{output_prefix}_tree.txt"
        ]

        file_status = {}

        for file_path in expected_files:
            exists = os.path.exists(file_path)
            file_status[file_path] = exists

            if exists:
                file_size = os.path.getsize(file_path)
                self.logger.info(f"  ✓ {os.path.basename(file_path)} (大小: {file_size:,} bytes)")
            else:
                self.logger.warning(f"  ✗ {os.path.basename(file_path)} 不存在")

        # 检查主要结果文件 | Check main result file
        main_result = f"{output_prefix}_BBAA.txt"
        if os.path.exists(main_result):
            with open(main_result, 'r') as f:
                trio_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
            self.logger.info(f"分析的三元组数量: {trio_count}")

        return file_status
