"""
基因组分析工具类|Genome Analysis Utility Classes
"""

import logging
import sys
import os
import subprocess
import re
import fnmatch
import shutil
import tempfile
from pathlib import Path
from typing import Tuple, Optional, List


class GenomeAnalysisLogger:
    """基因组分析日志管理器|Genome Analysis Logger Manager"""

    def __init__(self, output_dir: str):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
        """
        self.output_dir = output_dir
        # 使用99_logs目录统一管理日志|Use 99_logs directory for centralized log management
        logs_dir = os.path.join(output_dir, '99_logs')
        os.makedirs(logs_dir, exist_ok=True)
        log_file = os.path.join(logs_dir, 'genomescope_pipeline.log')
        self.setup_logging(log_file)

    def setup_logging(self, log_file: str):
        """设置日志|Setup logging"""
        self.logger = logging.getLogger("GenomeAnalysis")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 不传播到父logger，避免重复|Don't propagate to parent logger to avoid duplicates

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别→超算.out文件|stdout handler - INFO level → job scheduler .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上→超算.err文件|stderr handler - WARNING+ → job scheduler .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别→本地日志文件|File handler - all levels → local log file
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


class GenomeScopeRunner:
    """GenomeScope运行器|GenomeScope Runner"""

    def __init__(self, logger, config):
        """
        初始化运行器|Initialize runner

        Args:
            logger: 日志对象|Logger object
            config: 配置对象|Configuration object
        """
        self.logger = logger
        self.config = config

    def run_jellyfish_count(self, fastq_files: list, output_prefix: str,
                           kmer_size: int, hash_size: str, threads: int) -> bool:
        """
        运行Jellyfish count|Run Jellyfish count

        Args:
            fastq_files: FASTQ文件列表|List of FASTQ files
            output_prefix: 输出文件前缀|Output file prefix
            kmer_size: K-mer大小|K-mer size
            hash_size: 哈希表大小|Hash size
            threads: 线程数|Number of threads

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤1/3: 运行Jellyfish Count|Step 1/3: Run Jellyfish Count")
        self.logger.info("=" * 60)
        self.logger.info(f"配置参数|Configuration parameters:")
        self.logger.info(f"  - K-mer大小|K-mer size: {kmer_size}")
        self.logger.info(f"  - 哈希表大小|Hash size: {hash_size}")
        self.logger.info(f"  - 线程数|Threads: {threads}")
        self.logger.info(f"  - 输入文件数|Input files: {len(fastq_files)}")

        jf_output = f"{output_prefix}.jf"

        # 检查输出文件是否已存在|Check if output file already exists
        if os.path.exists(jf_output):
            self.logger.info(f"输出文件已存在，跳过此步骤: {jf_output}")
            self.logger.info("如需重新运行，请删除现有文件")
            return True

        try:
            # 检查是否有压缩文件|Check for compressed files
            gzip_files = [f for f in fastq_files if f.endswith('.gz')]
            plain_files = [f for f in fastq_files if not f.endswith('.gz')]

            if gzip_files:
                self.logger.info("检测到压缩文件，使用zcat解压缩")
                cmd = f"zcat {' '.join(gzip_files)} {' '.join(plain_files)}|" \
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
        运行Jellyfish histo|Run Jellyfish histo

        Args:
            jf_file: Jellyfish输出文件|Jellyfish output file
            output_prefix: 输出文件前缀|Output file prefix
            threads: 线程数|Number of threads

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤2/3: 运行Jellyfish Histo")
        self.logger.info("=" * 60)

        histo_file = f"{output_prefix}.histo"

        # 检查输出文件是否已存在|Check if output file already exists
        if os.path.exists(histo_file):
            self.logger.info(f"输出文件已存在，跳过此步骤: {histo_file}")
            self.logger.info("如需重新运行，请删除现有文件")
            return True

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
                       output_dir: str, max_kmer_cov: int, ploidy: int) -> Optional[float]:
        """
        运行GenomeScope 2.0 (使用genomescope2命令行工具)|Run GenomeScope 2.0 (using genomescope2 CLI tool)

        Args:
            histo_file: Histogram文件|Histogram file
            kmer_size: K-mer大小|K-mer size
            read_length: 读长|Read length (未使用，保留兼容性)|Read length (unused, kept for compatibility)
            output_dir: 输出目录|Output directory
            max_kmer_cov: 最大覆盖度|Max coverage
            ploidy: 基因组倍性|Genome ploidy (1-6)

        Returns:
            k-mer coverage值或None|K-mer coverage value or None
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤3/3: 运行GenomeScope 2.0")
        self.logger.info("=" * 60)

        # 使用绝对路径|Use absolute paths
        histo_file_abs = os.path.abspath(histo_file)
        output_dir_abs = os.path.abspath(output_dir)

        # 检查GenomeScope是否已经运行完成|Check if GenomeScope is already completed
        model_file = os.path.join(output_dir_abs, 'model.txt')
        if os.path.exists(model_file):
            self.logger.info(f"GenomeScope输出已存在，跳过运行步骤: {model_file}")
            self.logger.info("如需重新运行，请删除现有输出目录")
            # 直接从已有结果中提取kcov|Extract kcov directly from existing results
            kcov = self._extract_kcov(output_dir_abs)
            if kcov:
                self.logger.info(f"从已有结果中提取k-mer coverage: {kcov}")
                return kcov
            else:
                self.logger.warning("已有结果中未能提取kcov值")

        try:
            self.logger.info(f"Histogram文件|Histogram file: {histo_file_abs}")
            self.logger.info(f"输出目录|Output directory: {output_dir_abs}")

            # 构建conda run命令|Build conda run command
            # 按照开发规范：使用conda run -n <env_name> <command>|Following dev guide: use conda run -n <env_name> <command>
            cmd = [
                'conda', 'run', '-n', self.config.genomescope_env,
                'genomescope2',
                '-i', histo_file_abs,
                '-o', output_dir_abs,
                '-k', str(kmer_size),
                '-p', str(ploidy)
            ]

            # 添加可选参数|Add optional parameters
            if max_kmer_cov > 0:
                cmd.extend(['-m', str(max_kmer_cov)])

            self.logger.info(f"GenomeScope 2.0分析|GenomeScope 2.0 analyzing: {histo_file_abs}")
            self.logger.info(f"参数|Parameters: k={kmer_size}, p={ploidy}, outdir={output_dir_abs}")
            self.logger.info(f"命令|Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            # 记录命令输出|Log command output
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if line.strip():
                        self.logger.info(line)

            if result.stderr:
                # 记录stderr但不作为错误处理|Log stderr but not treat as error
                for line in result.stderr.split('\n'):
                    if line.strip() and 'Loading' not in line and 'Setting' not in line and 'Mirrors' not in line:
                        # 过滤掉R启动信息|Filter out R startup messages
                        self.logger.debug(line)

            # 检查GenomeScope是否成功生成结果|Check if GenomeScope successfully generated results
            # 不依赖returncode，而是检查输出文件|Don't rely on returncode, check output files instead
            kcov = self._extract_kcov(output_dir_abs)
            if kcov:
                self.logger.info(f"成功提取k-mer coverage|Successfully extracted k-mer coverage: {kcov}")

                # 即使returncode非0，如果结果文件存在且有效，也认为成功
                # Even if returncode != 0, if result file exists and valid, consider it successful
                if result.returncode != 0:
                    self.logger.warning(f"命令返回码非0({result.returncode})，但结果文件有效|Command return code non-zero({result.returncode}), but result file is valid")

                return kcov
            else:
                self.logger.error("未能从GenomeScope输出中提取kcov值|Failed to extract kcov from GenomeScope output")
                if result.returncode != 0:
                    self.logger.error(f"GenomeScope 2.0失败|GenomeScope 2.0 failed: returncode={result.returncode}")
                    if result.stderr:
                        self.logger.error(f"错误信息|Error message: {result.stderr[:500]}")
                return None

        except Exception as e:
            self.logger.error(f"运行GenomeScope时出错|Error running GenomeScope: {str(e)}", exc_info=True)
            return None

    def _extract_kcov(self, output_dir: str) -> Optional[float]:
        """
        从GenomeScope输出中提取kcov|Extract kcov from GenomeScope output

        Args:
            output_dir: GenomeScope输出目录|GenomeScope output directory

        Returns:
            kcov值或None|kcov value or None
        """
        try:
            # 首先尝试从model.txt中提取|Try to extract from model.txt first
            model_file = os.path.join(output_dir, 'model.txt')
            if os.path.exists(model_file):
                with open(model_file, 'r') as f:
                    for line in f:
                        # 查找类似 "kmercov 1.705e+01" 的行
                        if line.strip().startswith('kmercov'):
                            # 使用正则表达式提取kmercov后面的数值
                            # 格式: kmercov 1.705e+01  2.005e-02  ...
                            match = re.search(r'kmercov\s+([\d.e+-]+)', line)
                            if match:
                                try:
                                    kcov_value = float(match.group(1))
                                    self.logger.info(f"从model.txt中提取k-mer coverage: {kcov_value}")
                                    return kcov_value
                                except (ValueError, IndexError):
                                    continue

            # 如果model.txt中没有找到，尝试从summary.txt中提取
            # If not found in model.txt, try to extract from summary.txt
            summary_file = os.path.join(output_dir, 'summary.txt')
            if os.path.exists(summary_file):
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
                                    kcov_value = float(match.group(1))
                                    self.logger.info(f"从summary.txt中提取k-mer coverage: {kcov_value}")
                                    return kcov_value
                                except (ValueError, IndexError):
                                    continue

            return None

        except Exception as e:
            self.logger.error(f"读取GenomeScope输出时出错: {str(e)}")
            return None


class SmudgeplotRunner:
    """Smudgeplot运行器|Smudgeplot Runner"""

    def __init__(self, logger):
        """
        初始化运行器|Initialize runner

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def run_fastk(self, fastq_files: list, fastk_table: str,
                  kmer_size: int, threads: int, memory: str = "16G",
                  sample_name: str = "sample") -> bool:
        """
        运行FastK|Run FastK

        Args:
            fastq_files: FASTQ文件列表|List of FASTQ files
            fastk_table: FastK表输出路径|FastK table output path
            kmer_size: K-mer大小|K-mer size
            threads: 线程数|Number of threads
            memory: 内存大小|Memory size (e.g., "16G")
            sample_name: 样本名称|Sample name (for temp directory)

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤3.5/5: 运行FastK生成FastK表")
        self.logger.info("=" * 60)
        self.logger.info(f"配置参数:")
        self.logger.info(f"  - K-mer大小: {kmer_size}")
        self.logger.info(f"  - 线程数: {threads}")
        self.logger.info(f"  - 内存: {memory}")
        self.logger.info(f"  - 输入文件数: {len(fastq_files)}")

        # 检查FastK表是否已存在|Check if FastK table already exists
        # FastK生成.ktab和.hist文件，检查.ktab文件即可|FastK generates .ktab and .hist files, check .ktab is enough
        fastk_table_ktab = f"{fastk_table}.ktab"
        if os.path.exists(fastk_table_ktab):
            self.logger.info(f"FastK表已存在，跳过此步骤: {fastk_table_ktab}")
            self.logger.info("如需重新运行，请删除现有文件")
            return True

        # FastK的-M参数需要纯数字，不带单位
        # 提取内存数值|Extract memory value
        memory_value = memory.rstrip('Gg')
        try:
            memory_int = int(memory_value)
        except ValueError:
            self.logger.warning(f"无法解析内存值 '{memory}'，使用默认值16")
            memory_int = 16

        # 检查是否有.gz文件，如果有则先解压到系统临时目录
        # Check if there are .gz files, if so decompress to system temp directory first
        import tempfile
        temp_base_dir = tempfile.mkdtemp(prefix=f"genomescope_fastk_{sample_name}_")
        self.logger.debug(f"创建临时目录: {temp_base_dir}")
        decompressed_files = []

        try:
            for fastq_file in fastq_files:
                if fastq_file.endswith('.gz'):
                    # 解压到临时目录|Decompress to temp directory
                    decompressed_file = os.path.join(temp_base_dir, os.path.basename(fastq_file[:-3]))
                    self.logger.info(f"解压文件到临时目录: {fastq_file} -> {decompressed_file}")
                    with open(decompressed_file, 'wb') as f_out:
                        result = subprocess.run(
                            ['zcat', fastq_file],
                            stdout=f_out,
                            stderr=subprocess.PIPE
                        )
                        if result.returncode != 0:
                            self.logger.error(f"解压失败: {result.stderr.decode()}")
                            return False
                    decompressed_files.append(decompressed_file)
                else:
                    decompressed_files.append(fastq_file)

            # 使用解压后的文件（如果有）|Use decompressed files (if any)
            fastq_files_to_use = decompressed_files if decompressed_files else fastq_files

            # 构建FastK命令
            # FastK的参数格式：-t threads -k kmer_size -M memory -T threads
            # 参数和值需要连在一起，如 -t12 -k21 -M16 -T12
            cmd = [
                'FastK',
                '-v',
                f'-t{threads}',
                f'-k{kmer_size}',
                f'-M{memory_int}',
                f'-T{threads}',
            ] + fastq_files_to_use + [f'-N{fastk_table}']

            self.logger.info(f"命令: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"FastK失败: {result.stderr}")
                # 清理临时文件|Clean up temp files
                shutil.rmtree(temp_base_dir, ignore_errors=True)
                return False

            self.logger.info("FastK步骤成功完成")

            # 清理临时文件|Clean up temporary files
            self.logger.info(f"清理临时目录: {temp_base_dir}")
            shutil.rmtree(temp_base_dir, ignore_errors=True)

            return True

        except Exception as e:
            self.logger.error(f"运行FastK时出错: {str(e)}", exc_info=True)
            # 确保清理临时文件|Ensure cleanup temp files
            if 'temp_base_dir' in locals() and os.path.exists(temp_base_dir):
                shutil.rmtree(temp_base_dir, ignore_errors=True)
            return False

    def run_smudgeplot(self, fastk_table: str, output_dir: str, sample_name: str,
                       kcov: float, kmer_size: int, threads: int) -> bool:
        """
        运行Smudgeplot|Run Smudgeplot

        Args:
            fastk_table: FastK表文件|FastK table file
            output_dir: 输出目录|Output directory
            sample_name: 样本名|Sample name
            kcov: K-mer coverage|K-mer coverage
            kmer_size: K-mer大小|K-mer size
            threads: 线程数|Number of threads

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤4/5: 运行Smudgeplot Hetmers")
        self.logger.info("=" * 60)
        self.logger.info(f"配置参数:")

        # 计算阈值|Calculate threshold
        threshold = int(kcov * 0.5)
        self.logger.info(f"  - K-mer coverage: {kcov}")
        self.logger.info(f"  - K-mer大小: {kmer_size}")
        self.logger.info(f"  - 错误k-mer阈值(-L): {threshold} (计算公式: int({kcov} * 0.5) = {threshold})")
        self.logger.info(f"    注: 根据Smudgeplot官方文档，低端阈值设为单倍体k-mer coverage的50%")

        # 检查Logex是否可用（smudgeplot的C后端）|Check if Logex is available (smudgeplot C backend)
        import shutil
        logex_path = shutil.which('Logex')
        smudgeplot_dir = None  # 保存smudgeplot目录，用于后续subprocess调用|Save smudgeplot dir for later subprocess calls

        # 先尝试找到smudgeplot的安装目录|First try to find smudgeplot installation directory
        if not smudgeplot_dir:
            smudgeplot_path = shutil.which('smudgeplot')
            if smudgeplot_path:
                # 解析符号链接，获取真实路径|Resolve symbolic link to get real path
                real_smudgeplot_path = os.path.realpath(smudgeplot_path)
                smudgeplot_dir = os.path.dirname(real_smudgeplot_path)
                self.logger.debug(f"找到smudgeplot目录: {smudgeplot_dir}")
            else:
                # 如果PATH中找不到，尝试直接从conda环境查找
                # If not in PATH, try to find in conda environment directly
                import os
                conda_prefix = os.environ.get('CONDA_PREFIX')
                if conda_prefix:
                    # 当前环境的bin目录|Current environment's bin directory
                    smudgeplot_in_conda = os.path.join(conda_prefix, 'bin', 'smudgeplot')
                    if os.path.exists(smudgeplot_in_conda):
                        smudgeplot_dir = os.path.join(conda_prefix, 'bin')
                        self.logger.debug(f"在当前conda环境中找到smudgeplot: {smudgeplot_dir}")
                    else:
                        # 尝试smudgeplot环境|Try smudgeplot environment
                        smudgeplot_env = os.path.join(os.path.dirname(conda_prefix), 'smudgeplot', 'bin', 'smudgeplot')
                        if os.path.exists(smudgeplot_env):
                            smudgeplot_dir = os.path.join(os.path.dirname(conda_prefix), 'smudgeplot', 'bin')
                            self.logger.debug(f"在smudgeplot环境中找到: {smudgeplot_dir}")

                # 如果还是找不到，尝试常见的conda路径|If still not found, try common conda paths
                if not smudgeplot_dir:
                    common_paths = [
                        '/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/smudgeplot',
                        os.path.expanduser('~/miniforge3/envs/smudgeplot/bin/smudgeplot'),
                        '/miniconda3/envs/smudgeplot/bin/smudgeplot',
                    ]
                    for path in common_paths:
                        if os.path.exists(path):
                            smudgeplot_dir = os.path.dirname(path)
                            self.logger.debug(f"在常见路径中找到smudgeplot: {smudgeplot_dir}")
                            break

        if not logex_path:
            self.logger.warning("Logex未在PATH中找到，正在尝试定位...")
            # 尝试在smudgeplot安装目录中查找|Try to find in smudgeplot installation directory
            try:
                if smudgeplot_dir:
                    logex_in_dir = os.path.join(smudgeplot_dir, 'Logex')
                    if os.path.exists(logex_in_dir):
                        self.logger.info(f"找到Logex: {logex_in_dir}")
                        logex_path = logex_in_dir
                    else:
                        self.logger.warning(f"在 {smudgeplot_dir} 中未找到Logex")
                else:
                    self.logger.warning("smudgeplot命令也未找到，无法定位Logex")
            except Exception as e:
                self.logger.warning(f"检测Logex时出错: {e}")

        if logex_path:
            self.logger.info(f"Logex可用|Logex available: {logex_path}")
        else:
            self.logger.error("Logex不可用！Smudgeplot的C后端未正确安装|Logex not available! Smudgeplot C backend not properly installed")
            self.logger.error("解决方案: 重新安装smudgeplot或手动编译C代码|Solution: Reinstall smudgeplot or compile C code manually")
            self.logger.error("参考: https://github.com/KamilSJaron/smudgeplot|Reference: https://github.com/KamilSJaron/smudgeplot")
            return False

        # 创建输出目录|Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # 输出文件前缀（使用样本名，smudgeplot会自动添加_smudgeplot后缀）|Output file prefix
        # 文件命名：样本名_smudgeplot.png（smudgeplot自动添加）|File naming: sample_smudgeplot.png (auto-added by smudgeplot)
        # 避免重复：不要用 sample_name.smudgeplot，否则会变成 sample.smudgeplot_smudgeplot.png |Avoid duplication
        output_prefix = os.path.join(output_dir, sample_name)
        hetmers_output = f"{output_prefix}.kmerpairs"
        smu_file = f"{hetmers_output}.smu"

        # 检查.smu文件是否已存在|Check if .smu file already exists
        if os.path.exists(smu_file):
            self.logger.info(f"Hetmers输出已存在，跳过此步骤: {smu_file}")
            self.logger.info("如需重新运行，请删除现有文件")
        else:
            try:
                # 步骤1: 运行smudgeplot hetmers
                cmd = [
                    'smudgeplot', 'hetmers',
                    '-L', str(threshold),  # 使用kcov的50%作为错误k-mer下限阈值
                    '-t', str(threads),
                    '-o', hetmers_output,
                    '--verbose',
                    fastk_table
                ]
                self.logger.info(f"命令: {' '.join(cmd)}")

                # 准备环境变量：如果smudgeplot_dir存在，临时添加到PATH
                # Prepare environment: if smudgeplot_dir exists, temporarily add to PATH
                env = os.environ.copy()
                if smudgeplot_dir:
                    env['PATH'] = smudgeplot_dir + ':' + env.get('PATH', '')

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    env=env  # 临时设置环境变量|Temporarily set environment variables
                )

                if result.returncode != 0:
                    self.logger.error(f"Smudgeplot hetmers失败: {result.stderr}")
                    return False

                # 记录命令输出|Log command output
                if result.stdout:
                    self.logger.debug(f"Smudgeplot hetmers stdout:\n{result.stdout}")
                if result.stderr:
                    self.logger.info(f"Smudgeplot hetmers stderr:\n{result.stderr}")

                self.logger.info("Hetmers步骤成功完成")

            except Exception as e:
                self.logger.error(f"运行Hetmers时出错: {str(e)}", exc_info=True)
                return False

        # 步骤2: 运行smudgeplot plot
        if not os.path.exists(smu_file):
            self.logger.error(f"未找到.smu文件: {smu_file}")
            self.logger.warning("可能的原因:")
            self.logger.warning("  1. 基因组为纯合二倍体（homozygous diploid），异质性低")
            self.logger.warning("  2. K-mer coverage 过低，无法检测到足够的异质 k-mer 对")
            self.logger.warning("  3. 阈值设置过高（当前阈值: {}）".format(threshold))
            self.logger.info("这是正常现象，对于纯合基因组 Smudgeplot 可能无法生成图形")
            return False

        self.logger.info("=" * 60)
        self.logger.info("步骤5/5: 运行Smudgeplot Plot")
        self.logger.info("=" * 60)

        # 检查smudgeplot输出是否已存在|Check if smudgeplot output already exists
        plot_output_prefix = output_prefix
        # 检查是否有生成的图形文件（PNG或PDF）|Check if plot files (PNG or PDF) are generated
        import glob
        existing_plots = glob.glob(f"{plot_output_prefix}*.png") + glob.glob(f"{plot_output_prefix}*.pdf")
        if existing_plots:
            self.logger.info(f"Smudgeplot输出已存在，跳过此步骤")
            self.logger.info(f"找到文件|Found files: {', '.join([os.path.basename(f) for f in existing_plots])}")
            self.logger.info("如需重新运行，请删除现有输出文件")
            self.logger.info("Smudgeplot分析成功完成")
            return True

        try:
            cmd = [
                'smudgeplot', 'all',
                '-o', plot_output_prefix,
                '-cov', str(kcov),  # 使用从GenomeScope提取的kcov值
                smu_file
            ]
            self.logger.info(f"命令: {' '.join(cmd)}")

            # 准备环境变量：如果smudgeplot_dir存在，临时添加到PATH
            # Prepare environment: if smudgeplot_dir exists, temporarily add to PATH
            env = os.environ.copy()
            if smudgeplot_dir:
                env['PATH'] = smudgeplot_dir + ':' + env.get('PATH', '')

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env  # 临时设置环境变量|Temporarily set environment variables
            )

            if result.returncode != 0:
                self.logger.error(f"Smudgeplot all失败: {result.stderr}")
                return False

            self.logger.info("Smudgeplot分析成功完成")
            return True

        except Exception as e:
            self.logger.error(f"运行Smudgeplot时出错: {str(e)}", exc_info=True)
            return False


def extract_sample_name(filename: str, pattern: str) -> str:
    """
    从文件名中提取样本名，支持通配符模式|Extract sample name from filename, supports wildcard patterns

    Args:
        filename: 文件名|Filename
        pattern: 文件名模式（可包含通配符）|Filename pattern (may contain wildcards)

    Returns:
        提取的样本名|Extracted sample name
    """
    # 如果模式不包含通配符，直接替换|If pattern doesn't contain wildcard, replace directly
    if '*' not in pattern and '?' not in pattern:
        return filename.replace(pattern, "")

    # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
    # 将 * 替换为 (.+) 捕获组以提取样本名|Replace * with (.+) capture group to extract sample name
    regex_pattern = fnmatch.translate(pattern)
    # 修改正则表达式，将通配符部分设为捕获组|Modify regex to make wildcard part a capture group
    regex_pattern = regex_pattern.replace(r'\.\*.*', r'(.+)\..*')
    regex_pattern = regex_pattern.replace(r'\Z', r'\Z')

    match = re.match(regex_pattern, filename)
    if match:
        return match.group(1)

    # 如果正则匹配失败，回退到简单替换|If regex match fails, fallback to simple replacement
    # 去除已知的后缀|Remove known suffix
    name = filename
    # 尝试去除模式中通配符后面的部分|Try to remove part after wildcard in pattern
    if '*' in pattern:
        suffix = pattern.split('*', 1)[1]
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


class SampleFinder:
    """样本查找器|Sample Finder for genome analysis"""

    def __init__(self, logger, read1_suffix: str = "*_1.clean.fq.gz"):
        """
        初始化样本查找器|Initialize sample finder

        Args:
            logger: 日志对象|Logger object
            read1_suffix: Read1文件后缀模式|Read1 file suffix pattern
        """
        self.logger = logger
        self.read1_suffix = read1_suffix

    def find_samples(self, input_path: str) -> List[Tuple[str, List[Path]]]:
        """
        查找样本及其FASTQ文件|Find samples and their FASTQ files

        Args:
            input_path: 输入路径(文件或目录)|Input path (file or directory)

        Returns:
            样本列表，每个元素为 (样本名, FASTQ文件列表)|List of samples (sample_name, fastq_files_list)
        """
        samples = []

        # 检查输入是文件还是目录|Check if input is file or directory
        if os.path.isfile(input_path):
            self.logger.info(f"单文件输入模式|Single file input mode: {input_path}")
            # 提取样本名|Extract sample name
            input_file = Path(input_path)
            sample_name = self._extract_name_from_file(input_file.name)
            samples.append((sample_name, [input_file]))
            return samples

        # 目录模式|Directory mode
        self.logger.info(f"目录输入模式，正在查找FASTQ文件|Directory input mode, searching FASTQ files: {input_path}")
        input_dir = Path(input_path)

        # 查找所有FASTQ文件|Find all FASTQ files
        fastq_files = list(input_dir.rglob('*.fastq')) + \
                      list(input_dir.rglob('*.fq')) + \
                      list(input_dir.rglob('*.fastq.gz')) + \
                      list(input_dir.rglob('*.fq.gz'))

        if not fastq_files:
            self.logger.error(f"未找到任何FASTQ文件|No FASTQ files found: {input_path}")
            return []

        # 按样品分组|Group by sample
        sample_dict = self._group_files_by_sample(fastq_files)

        # 转换为列表格式|Convert to list format
        samples = [(name, files) for name, files in sorted(sample_dict.items())]

        self.logger.info(f"找到 {len(samples)} 个样本|Found {len(samples)} samples:")
        for sample_name, files in samples:
            self.logger.info(f"  - {sample_name}: {len(files)} 个文件")

        return samples

    def _group_files_by_sample(self, fastq_files: List[Path]) -> dict:
        """
        将FASTQ文件按样品分组|Group FASTQ files by sample

        Args:
            fastq_files: FASTQ文件列表|List of FASTQ files

        Returns:
            样本字典 {样本名: [文件列表]}|Sample dictionary {sample_name: [file_list]}
        """
        sample_dict = {}

        for fastq_file in fastq_files:
            # 使用 read1_suffix 模式提取样品名|Use read1_suffix pattern to extract sample name
            sample_name = self._extract_name_using_suffix(fastq_file.name)

            if sample_name not in sample_dict:
                sample_dict[sample_name] = []
            sample_dict[sample_name].append(fastq_file)

        return sample_dict

    def _extract_name_using_suffix(self, filename: str) -> str:
        """
        使用 read1_suffix 模式提取样品名|Extract sample name using read1_suffix pattern

        Args:
            filename: 文件名|Filename

        Returns:
            样品名|Sample name
        """
        # 从 read1_suffix 中提取模式|Extract pattern from read1_suffix
        # 例如 *_1.clean.fq.gz → 需要匹配 *_1 或 *_2
        # For example *_1.clean.fq.gz → need to match *_1 or *_2

        # 构建可能的模式列表|Build list of possible patterns
        patterns = []
        patterns.append(self.read1_suffix)  # 原始模式，如 *_1.clean.fq.gz

        # 生成配对模式|Generate paired patterns
        if '_1' in self.read1_suffix:
            patterns.append(self.read1_suffix.replace('_1', '_2'))
        if '_R1' in self.read1_suffix:
            patterns.append(self.read1_suffix.replace('_R1', '_R2'))

        # 尝试每个模式|Try each pattern
        for pattern in patterns:
            sample_name = self._match_pattern(filename, pattern)
            if sample_name:
                # 清理可能的 .clean 等中间标识|Clean up intermediate identifiers like .clean
                sample_name = self._clean_sample_name(sample_name)
                if sample_name:
                    return sample_name

        # 回退到基础方法|Fallback to base method
        return self._extract_sample_base_name(filename)

    def _match_pattern(self, filename: str, pattern: str) -> Optional[str]:
        """
        使用模式匹配文件名并提取样品名|Match filename with pattern and extract sample name

        Args:
            filename: 文件名|Filename
            pattern: 模式（如 *_1.clean.fq.gz）|Pattern (e.g. *_1.clean.fq.gz)

        Returns:
            提取的样品名或None|Extracted sample name or None
        """
        if '*' not in pattern:
            # 没有通配符，直接替换|No wildcard, direct replacement
            if filename.endswith(pattern):
                return filename[:-len(pattern)]
            return None

        # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
        # 将 * 替换为捕获组 (.*)|Replace * with capture group (.*)
        regex_pattern = pattern.replace('.', r'\.')  # 转义点
        regex_pattern = regex_pattern.replace('*', '(.*)')  # 将 * 替换为捕获组
        regex_pattern = '^' + regex_pattern + '$'  # 添加首尾锚点

        match = re.match(regex_pattern, filename)
        if match:
            # 返回第一个捕获组（通配符部分）|Return first capture group (wildcard part)
            return match.group(1)

        return None

    def _clean_sample_name(self, name: str) -> str:
        """
        清理样品名，去除常见中间标识|Clean sample name, remove common intermediate identifiers

        Args:
            name: 原始样品名|Original sample name

        Returns:
            清理后的样品名|Cleaned sample name
        """
        # 去除常见的中间标识|Remove common intermediate identifiers
        # 这些标识通常出现在配对标识之后、扩展名之前
        # These identifiers usually appear between pair identifier and extension
        patterns_to_remove = [
            '\.clean$',   # .clean
            '\.trimmed$', # .trimmed
            '\.filtered$', # .filtered
        ]

        for pattern in patterns_to_remove:
            name = re.sub(pattern, '', name, flags=re.IGNORECASE)

        return name

    def _extract_sample_base_name(self, filename: str) -> str:
        """
        从文件名提取样品基础名称（去除配对标识）|Extract sample base name from filename (remove pair identifiers)

        Args:
            filename: 文件名|Filename

        Returns:
            样品基础名称|Sample base name
        """
        name = filename

        # 去除压缩文件扩展名|Remove compressed file extensions
        for ext in ['.fastq.gz', '.fq.gz']:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        else:
            # 尝试单个扩展名|Try single extension
            for ext in ['.fq', '.fastq', '.gz']:
                if name.endswith(ext):
                    name = name[:-len(ext)]
                    break

        # 去除配对标识|Remove pair identifiers
        # 常见模式：_1, _2, .R1, .R2, _R1, _R2, .r1, .r2, _r1, _r2
        # Common patterns: _1, _2, .R1, .R2, _R1, _R2, .r1, .r2, _r1, _r2
        pair_patterns = [
            '_1$', '_2$',          # _1, _2
            '_R1$', '_R2$',        # _R1, _R2
            '_r1$', '_r2$',        # _r1, _r2
            '\.R1\.', '\.R2\.',    # .R1., .R2. (在扩展名之前)
            '\.r1\.', '\.r2\.',    # .r1., .r2.
            'read1$', 'read2$',    # read1, read2
        ]

        for pattern in pair_patterns:
            if re.search(pattern, name, re.IGNORECASE):
                name = re.sub(pattern, '', name, flags=re.IGNORECASE)
                break

        # 清理可能的连续点或下划线|Clean up possible consecutive dots or underscores
        name = re.sub(r'_+', '_', name).strip('_')
        name = re.sub(r'\.+', '.', name).strip('.')

        return name

    def _extract_name_from_file(self, filename: str) -> str:
        """
        从文件名提取样本名（保留用于向后兼容）|Extract sample name from filename (kept for backward compatibility)

        Args:
            filename: 文件名|Filename

        Returns:
            样本名|Sample name
        """
        return self._extract_sample_base_name(filename)
