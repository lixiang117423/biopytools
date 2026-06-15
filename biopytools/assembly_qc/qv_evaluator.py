"""
QV质量值评估模块（完全独立实现）|QV Quality Value Evaluation Module (Fully Independent)

直接调用meryl/merqury.sh命令行工具
Directly calls meryl/merqury.sh command-line tools
"""

import os
import subprocess
import re
import glob
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from .utils import get_conda_env_from_path, build_conda_command


class QVEvaluator:
    """QV质量值评估器|QV Quality Value Evaluator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.working_dir = Path(self.config.qv_output_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

        # 从conda环境路径提取环境名|Extract env name from conda env path
        self.merqury_env_name = get_conda_env_from_path(self.config.conda_env_merqury)
        self.merqury_path = Path(os.path.expanduser(self.config.conda_env_merqury)).parent.parent

    def evaluate(self) -> Optional[Dict[str, Any]]:
        """运行QV评估|Run QV evaluation"""
        if not self.config.enable_qv or self.config.skip_qv:
            self.logger.info("跳过QV评估|Skipping QV evaluation")
            return None

        self.logger.info("开始QV质量值评估|Starting QV quality value evaluation")

        # 存储所有QV结果|Store all QV results
        all_results = {}

        # NGS数据QV评估|NGS data QV evaluation
        if self.config.enable_qv_ngs and self.config.ngs_reads:
            self.logger.info("评估NGS数据QV|Evaluating NGS data QV")
            ngs_reads_dir = self.config.get_qv_ngs_reads()
            if ngs_reads_dir:
                ngs_qv_file = self.working_dir / "qv_result_ngs.qv"
                if self.config.resume and ngs_qv_file.exists():
                    self.logger.info("NGS数据QV评估已完成，跳过|NGS QV evaluation already completed, skipping")
                    ngs_result = self._parse_results(ngs_qv_file, prefix="ngs_")
                else:
                    ngs_result = self._run_qv_pipeline(ngs_reads_dir, "ngs")
                if ngs_result:
                    all_results.update(ngs_result)

        # Long-read数据QV评估|Long-read data QV evaluation
        if self.config.enable_qv_long_read and self.config.long_reads:
            self.logger.info("评估Long-read数据QV|Evaluating long-read data QV")
            long_reads_dir = self.config.get_qv_long_read_reads()
            if long_reads_dir:
                long_qv_file = self.working_dir / "qv_result_long_read.qv"
                if self.config.resume and long_qv_file.exists():
                    self.logger.info("Long-read数据QV评估已完成，跳过|Long-read QV evaluation already completed, skipping")
                    long_result = self._parse_results(long_qv_file, prefix="long_read_")
                else:
                    long_result = self._run_qv_pipeline(long_reads_dir, "long_read")
                if long_result:
                    all_results.update(long_result)

        if not all_results:
            self.logger.error("未提供QV评估所需的reads数据|No reads data provided for QV evaluation")
            return None

        return all_results

    def _run_qv_pipeline(self, reads_dir: str, data_type: str) -> Optional[Dict[str, Any]]:
        """
        运行QV流程|Run QV pipeline

        Args:
            reads_dir: reads目录|Reads directory
            data_type: 数据类型（ngs或long_read）|Data type (ngs or long_read)
        """
        try:
            # 1. 构建Meryl数据库|Build Meryl database
            meryl_db = self._build_meryl_db(reads_dir, data_type)
            if not meryl_db:
                return None

            # 2. 计算QV值|Calculate QV value
            return self._calculate_qv(meryl_db, data_type)

        except Exception as e:
            self.logger.error(f"{data_type.upper()}数据QV评估异常|{data_type.upper()} data QV evaluation error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _build_meryl_db(self, reads_dir: str, data_type: str) -> Optional[str]:
        """
        构建Meryl k-mer数据库|Build Meryl k-mer database

        Args:
            reads_dir: reads目录|Reads directory
            data_type: 数据类型（ngs或long_read）|Data type (ngs or long_read)
        """
        self.logger.info(f"构建Meryl k-mer数据库 ({data_type.upper()})|Building Meryl k-mer database ({data_type.upper()})")

        # 查找FASTQ文件|Find FASTQ files
        fastq_files = self._discover_fastq_files(reads_dir)
        if not fastq_files:
            self.logger.error(f"未找到FASTQ文件|No FASTQ files found in: {reads_dir}")
            return None

        # 确定k值|Determine k-mer size
        k = self.config.qv_kmer_size
        if k is None:
            k = 21  # 默认值|Default value
            self.logger.info(f"使用默认k值|Using default k-mer size: {k}")
        else:
            self.logger.info(f"使用指定k值|Using specified k-mer size: {k}")

        meryl_db = self.working_dir / f"{data_type}_reads.meryl"

        # 检查是否已存在|Check if already exists
        if meryl_db.exists():
            self.logger.info(f"Meryl数据库已存在，跳过构建|Meryl database already exists, skipping: {meryl_db}")
            return str(meryl_db)

        # 构建meryl命令|Build meryl command
        # 检测是否有配对数据|Check if paired-end data
        pairs, single_end = self._find_paired_files(fastq_files)

        if pairs:
            # 使用配对数据|Use paired-end data
            r1, r2 = pairs[0]
            args = [
                "count",
                f"k={k}",
                "output", str(meryl_db),
                r1,
                r2
            ]
        elif single_end:
            # 使用单端数据|Use single-end data
            fastq_list = single_end[:5]  # 限制文件数量避免命令过长
            args = [
                "count",
                f"k={k}",
                "output", str(meryl_db),
            ] + fastq_list
        else:
            self.logger.error("未找到有效的FASTQ文件|No valid FASTQ files found")
            return None

        # 使用conda run调用meryl|Use conda run to call meryl
        cmd = build_conda_command(self.merqury_env_name, "meryl", args)

        # 设置环境变量|Set environment variables
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(self.config.qv_threads)

        # 执行命令|Execute command
        try:
            self.logger.info(f"运行Meryl|Running Meryl: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                cwd=str(self.working_dir),
                env=env,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                self.logger.error(f"Meryl数据库构建失败|Meryl database build failed")
                self.logger.error(f"STDERR: {result.stderr}")
                return None

            self.logger.info(f"Meryl数据库构建完成|Meryl database built successfully: {meryl_db}")
            return str(meryl_db)

        except Exception as e:
            self.logger.error(f"构建Meryl数据库异常|Error building Meryl database: {e}")
            return None

    def _calculate_qv(self, meryl_db: str, data_type: str) -> Optional[Dict[str, Any]]:
        """
        计算QV值|Calculate QV value

        Args:
            meryl_db: Meryl数据库路径|Meryl database path
            data_type: 数据类型（ngs或long_read）|Data type (ngs or long_read)
        """
        self.logger.info(f"计算QV值 ({data_type.upper()})|Calculating QV value ({data_type.upper()})")

        # 查找qv.sh脚本|Find qv.sh script
        qv_script_candidates = [
            self.merqury_path / "share" / "merqury" / "eval" / "qv.sh",
            self.merqury_path / "eval" / "qv.sh",
            Path(os.path.expanduser(self.config.conda_env_merqury)) / "share" / "merqury" / "eval" / "qv.sh",
        ]

        qv_script = None
        for candidate in qv_script_candidates:
            if candidate.exists():
                qv_script = candidate
                self.logger.debug(f"找到qv.sh脚本|Found qv.sh script: {qv_script}")
                break

        if not qv_script:
            self.logger.error(f"未找到qv.sh脚本|qv.sh script not found, searched paths:")
            for candidate in qv_script_candidates:
                self.logger.error(f"  - {candidate}")
            self.logger.info("尝试使用系统PATH中的qv.sh|Trying qv.sh from system PATH")
            qv_script_str = "qv.sh"
        else:
            qv_script_str = str(qv_script)

        output_prefix = self.working_dir / f"qv_result_{data_type}"

        # 设置MERQURY环境变量|Set MERQURY environment variable
        # qv.sh脚本需要这个变量来找到util.sh|qv.sh script needs this to find util.sh
        merqury_base = os.path.join(os.path.expanduser(self.config.conda_env_merqury), "share", "merqury")
        if not os.path.exists(merqury_base):
            # 尝试其他可能的位置|Try other possible locations
            merqury_base = os.path.join(os.path.expanduser(self.config.conda_env_merqury), "lib")

        env = os.environ.copy()
        env['MERQURY'] = merqury_base
        self.logger.debug(f"设置MERQURY环境变量|Set MERQURY env var: {merqury_base}")

        # 添加merqury的bin目录到PATH，使qv.sh能找到meryl-lookup
        # Add merqury bin directory to PATH so qv.sh can find meryl-lookup
        merqury_bin = os.path.join(os.path.expanduser(self.config.conda_env_merqury), "bin")
        env['PATH'] = f"{merqury_bin}:{env.get('PATH', '')}"
        self.logger.debug(f"设置PATH包含merqury bin|Set PATH to include merqury bin: {merqury_bin}")

        # 构建命令|Build command
        cmd = [
            "bash",
            qv_script_str,
            meryl_db,
            self.config.genome,
            str(output_prefix)
        ]

        try:
            self.logger.info(f"运行QV计算|Running QV calculation: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                cwd=str(self.working_dir),
                capture_output=True,
                text=True,
                check=False,
                env=env  # 使用设置好MERQURY的环境变量|Use env with MERQURY set
            )

            # 记录完整输出|Log full output for debugging
            if result.stdout:
                self.logger.debug(f"QV计算STDOUT|QV calculation STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"QV计算STDERR|QV calculation STDERR:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"QV计算失败|QV calculation failed, return code: {result.returncode}")
                # 检查是否有常见错误信息|Check for common error messages
                if "meryl" in result.stderr.lower():
                    self.logger.error("Meryl数据库相关错误|Meryl database related error")
                elif "genome" in result.stderr.lower():
                    self.logger.error("基因组文件相关错误|Genome file related error")
                return None

            # 解析结果|Parse results
            qv_file = str(output_prefix) + ".qv"

            # 列出工作目录内容用于调试|List working directory for debugging
            self.logger.debug(f"工作目录内容|Working directory contents:")
            try:
                for item in self.working_dir.iterdir():
                    self.logger.debug(f"  - {item.name}")
            except Exception as e:
                self.logger.debug(f"  无法列出目录|Cannot list directory: {e}")

            if os.path.exists(qv_file):
                return self._parse_results(Path(qv_file), prefix=f"{data_type}_")
            else:
                self.logger.error(f"QV输出文件未生成|QV output file not generated: {qv_file}")
                self.logger.info(f"请检查qv.sh脚本是否正确执行|Please verify qv.sh script executed correctly")
                # 尝试查找可能的输出文件|Try to find possible output files
                for ext in ['.qv', '.txt', '.log']:
                    possible_file = str(output_prefix) + ext
                    if os.path.exists(possible_file):
                        self.logger.info(f"发现可能的输出文件|Found possible output file: {possible_file}")
                return None

        except Exception as e:
            self.logger.error(f"QV计算异常|Error calculating QV: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _discover_fastq_files(self, reads_dir: str) -> List[str]:
        """发现FASTQ文件|Discover FASTQ files"""
        patterns = [
            "*.fq", "*.fq.gz",
            "*.fastq", "*.fastq.gz",
            "*_1.fq.gz", "*_2.fq.gz",
            "*_1.clean.fq.gz", "*_2.clean.fq.gz"
        ]

        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(reads_dir, pattern)))

        return sorted(list(set(files)))

    def _find_paired_files(self, fastq_files: List[str]) -> Tuple[List[tuple], List[str]]:
        """查找配对的FASTQ文件|Find paired FASTQ files"""
        pairs = []
        single = []

        processed = set()

        for f in fastq_files:
            if f in processed:
                continue

            basename = os.path.basename(f)

            # 尝试找到配对|Try to find pair
            if "_1." in basename or "_R1." in basename:
                r2 = f.replace("_1.", "_2.").replace("_R1.", "_R2.")
                if os.path.exists(r2):
                    pairs.append((f, r2))
                    processed.add(f)
                    processed.add(r2)
                    continue

            # 如果没有找到配对，当作单端处理
            if f not in processed:
                single.append(f)
                processed.add(f)

        return pairs, single

    def _parse_results(self, qv_file: Path, prefix: str = "") -> Optional[Dict[str, Any]]:
        """
        解析QV结果文件|Parse QV result file

        Args:
            qv_file: QV结果文件路径|QV result file path
            prefix: 结果键前缀（ngs_或long_read_）|Result key prefix (ngs_ or long_read_)
        """
        try:
            self.logger.info(f"解析QV结果|Parsing QV results: {qv_file}")

            with open(qv_file, 'r') as f:
                content = f.read()

            # QV结果格式：多列，最后一列是QV值
            # QV result format: multiple columns, last column is QV value
            lines = content.strip().split('\n')
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        # 格式：sequence_name N50/contigs total_bases QV_value error_rate
                        # Format: sequence_name N50/contigs total_bases QV_value error_rate
                        try:
                            qv_value = float(parts[3])
                            error_rate = float(parts[4]) if len(parts) > 4 else None

                            self.logger.info(f"QV质量值|QV value: {qv_value}")
                            if error_rate is not None:
                                self.logger.info(f"错误率|Error rate: {error_rate}")

                            return {
                                f'{prefix}qv_value': qv_value,
                                f'{prefix}error_rate': error_rate,
                                f'{prefix}sequence_name': parts[0],
                                f'{prefix}total_bases': parts[2],
                            }
                        except (ValueError, IndexError) as e:
                            self.logger.warning(f"解析QV行失败|Failed to parse QV line: {line}")
                            continue

            # 如果上述方法失败，尝试直接查找QV值
            # If above method fails, try to find QV value directly
            match = re.search(r'QV[:\s]+([\d.]+)', content, re.IGNORECASE)
            if match:
                qv_value = float(match.group(1))
                self.logger.info(f"QV质量值|QV value: {qv_value}")
                return {f'{prefix}qv_value': qv_value}

            self.logger.error("无法解析QV值|Cannot parse QV value")
            self.logger.debug(f"QV文件内容|QV file content:\n{content}")
            return None

        except Exception as e:
            self.logger.error(f"解析QV结果失败|Failed to parse QV results: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None
