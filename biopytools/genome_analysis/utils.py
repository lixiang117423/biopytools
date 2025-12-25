"""
Genome Analysis Utility Classes
基因组分析工具类
"""

import logging
import sys
import os
import subprocess
import re
from typing import Tuple, Optional


class GenomeAnalysisLogger:
    """基因组分析日志管理器 | Genome Analysis Logger Manager"""

    def __init__(self, output_dir: str):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            output_dir: 输出目录 | Output directory
        """
        self.output_dir = output_dir
        log_file = os.path.join(output_dir, 'genome_analysis.log')
        self.setup_logging(log_file)

    def setup_logging(self, log_file: str):
        """设置日志 | Setup logging"""
        self.logger = logging.getLogger("GenomeAnalysis")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式 | Log format
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler | File handler
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class GenomeScopeRunner:
    """GenomeScope运行器 | GenomeScope Runner"""

    def __init__(self, logger):
        """
        初始化运行器 | Initialize runner

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def run_jellyfish_count(self, fastq_files: list, output_prefix: str,
                           kmer_size: int, hash_size: str, threads: int) -> bool:
        """
        运行Jellyfish count | Run Jellyfish count

        Args:
            fastq_files: FASTQ文件列表 | List of FASTQ files
            output_prefix: 输出文件前缀 | Output file prefix
            kmer_size: K-mer大小 | K-mer size
            hash_size: 哈希表大小 | Hash size
            threads: 线程数 | Number of threads

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 1/3: 运行Jellyfish Count")
        self.logger.info("=" * 60)
        self.logger.info(f"配置参数:")
        self.logger.info(f"  - K-mer大小: {kmer_size}")
        self.logger.info(f"  - 哈希表大小: {hash_size}")
        self.logger.info(f"  - 线程数: {threads}")
        self.logger.info(f"  - 输入文件数: {len(fastq_files)}")

        jf_output = f"{output_prefix}.jf"

        try:
            # 检查是否有压缩文件 | Check for compressed files
            gzip_files = [f for f in fastq_files if f.endswith('.gz')]
            plain_files = [f for f in fastq_files if not f.endswith('.gz')]

            if gzip_files:
                self.logger.info("检测到压缩文件，使用zcat解压缩")
                cmd = f"zcat {' '.join(gzip_files)} {' '.join(plain_files)} | " \
                      f"jellyfish count -C -m {kmer_size} -s {hash_size} -t {threads} " \
                      f"-o {jf_output} /dev/fd/0"
                self.logger.info(f"命令: {cmd}")

                process = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                stdout, stderr = process.communicate()

                if process.returncode != 0:
                    self.logger.error(f"Jellyfish count失败: {stderr.decode('utf-8')}")
                    return False
            else:
                cmd = [
                    'jellyfish', 'count',
                    '-C',
                    '-m', str(kmer_size),
                    '-s', hash_size,
                    '-t', str(threads),
                    '-o', jf_output
                ] + fastq_files
                self.logger.info(f"命令: {' '.join(cmd)}")

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True
                )

                if result.returncode != 0:
                    self.logger.error(f"Jellyfish count失败: {result.stderr}")
                    return False

            self.logger.info("Jellyfish count步骤成功完成")
            return True

        except Exception as e:
            self.logger.error(f"运行Jellyfish count时出错: {str(e)}", exc_info=True)
            return False

    def run_jellyfish_histo(self, jf_file: str, output_prefix: str, threads: int) -> bool:
        """
        运行Jellyfish histo | Run Jellyfish histo

        Args:
            jf_file: Jellyfish输出文件 | Jellyfish output file
            output_prefix: 输出文件前缀 | Output file prefix
            threads: 线程数 | Number of threads

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 2/3: 运行Jellyfish Histo")
        self.logger.info("=" * 60)

        histo_file = f"{output_prefix}.histo"

        try:
            cmd = [
                'jellyfish', 'histo',
                '-t', str(threads),
                jf_file
            ]
            self.logger.info(f"命令: {' '.join(cmd)}")

            with open(histo_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True
                )

            if result.returncode != 0:
                self.logger.error(f"Jellyfish histo失败: {result.stderr}")
                return False

            self.logger.info("Jellyfish histo步骤成功完成")
            return True

        except Exception as e:
            self.logger.error(f"运行Jellyfish histo时出错: {str(e)}", exc_info=True)
            return False

    def run_genomescope(self, histo_file: str, kmer_size: int, read_length: int,
                       output_dir: str, max_kmer_cov: int, r_script: str) -> Optional[float]:
        """
        运行GenomeScope | Run GenomeScope

        Args:
            histo_file: Histogram文件 | Histogram file
            kmer_size: K-mer大小 | K-mer size
            read_length: 读长 | Read length
            output_dir: 输出目录 | Output directory
            max_kmer_cov: 最大覆盖度 | Max coverage
            r_script: R脚本路径 | R script path

        Returns:
            k-mer coverage值或None | K-mer coverage value or None
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 3/3: 运行GenomeScope")
        self.logger.info("=" * 60)

        try:
            cmd = [
                'Rscript', r_script,
                histo_file,
                str(kmer_size),
                str(read_length),
                output_dir,
                str(max_kmer_cov)
            ]
            self.logger.info(f"命令: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"GenomeScope失败: {result.stderr}")
                return None

            # 从summary.txt中提取kcov | Extract kcov from summary.txt
            summary_file = os.path.join(output_dir, 'summary.txt')
            if os.path.exists(summary_file):
                kcov = self._extract_kcov(summary_file)
                if kcov:
                    self.logger.info(f"提取的k-mer coverage: {kcov}")
                    return kcov
                else:
                    self.logger.warning("未能从summary.txt中提取kcov值")
            else:
                self.logger.warning(f"未找到summary.txt文件: {summary_file}")

            self.logger.info("GenomeScope分析成功完成")
            return None

        except Exception as e:
            self.logger.error(f"运行GenomeScope时出错: {str(e)}", exc_info=True)
            return None

    def _extract_kcov(self, summary_file: str) -> Optional[float]:
        """
        从summary.txt提取kcov | Extract kcov from summary.txt

        Args:
            summary_file: summary.txt文件路径 | summary.txt file path

        Returns:
            kcov值或None | kcov value or None
        """
        try:
            with open(summary_file, 'r') as f:
                for line in f:
                    # 查找类似 "Kmer_coverage: 53.4471" 的行
                    if 'Kmer_coverage' in line or 'Kmer coverage' in line:
                        # 尝试多种模式
                        match = re.search(r'Kmer_coverage[:\s]+([\d.]+)', line)
                        if not match:
                            match = re.search(r'Kmer coverage[:\s]+([\d.]+)', line)
                        if not match:
                            match = re.search(r'([\d.]+)', line)

                        if match:
                            try:
                                return float(match.group(1))
                            except (ValueError, IndexError):
                                continue

            return None

        except Exception as e:
            self.logger.error(f"读取summary.txt时出错: {str(e)}")
            return None


class SmudgeplotRunner:
    """Smudgeplot运行器 | Smudgeplot Runner"""

    def __init__(self, logger):
        """
        初始化运行器 | Initialize runner

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def run_smudgeplot(self, fastk_table: str, output_prefix: str,
                       kcov: float, kmer_size: int, threads: int) -> bool:
        """
        运行Smudgeplot | Run Smudgeplot

        Args:
            fastk_table: FastK表文件 | FastK table file
            output_prefix: 输出前缀 | Output prefix
            kcov: K-mer coverage | K-mer coverage
            kmer_size: K-mer大小 | K-mer size
            threads: 线程数 | Number of threads

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 4/5: 运行Smudgeplot Hetmers")
        self.logger.info("=" * 60)
        self.logger.info(f"配置参数:")
        self.logger.info(f"  - K-mer coverage: {kcov}")
        self.logger.info(f"  - K-mer大小: {kmer_size}")

        try:
            # 步骤1: 运行smudgeplot hetmers
            hetmers_output = f"{output_prefix}_kmerpairs"
            cmd = [
                'smudgeplot', 'hetmers',
                '-L', str(int(kcov * 0.5)),  # 使用kcov的50%作为下限
                '-t', str(threads),
                '-o', hetmers_output,
                fastk_table
            ]
            self.logger.info(f"命令: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"Smudgeplot hetmers失败: {result.stderr}")
                return False

            # 步骤2: 运行smudgeplot plot
            smu_file = f"{hetmers_output}.smu"
            if not os.path.exists(smu_file):
                self.logger.error(f"未找到.smu文件: {smu_file}")
                return False

            self.logger.info("=" * 60)
            self.logger.info("步骤 5/5: 运行Smudgeplot Plot")
            self.logger.info("=" * 60)

            plot_output = output_prefix
            cmd = [
                'smudgeplot', 'plot',
                '-o', plot_output,
                smu_file
            ]
            self.logger.info(f"命令: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"Smudgeplot plot失败: {result.stderr}")
                return False

            self.logger.info("Smudgeplot分析成功完成")
            return True

        except Exception as e:
            self.logger.error(f"运行Smudgeplot时出错: {str(e)}", exc_info=True)
            return False
