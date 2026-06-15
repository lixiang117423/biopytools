"""
Cactus泛基因组分析核心模块|Cactus Pangenome Analysis Core Module
"""

import os
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from typing import Optional

from .config import CactusConfig
from .utils import (
    CactusLogger,
    validate_singularity,
    validate_cactus_sif,
    check_genome_files,
    create_absolute_path_seqfile,
    build_singularity_command,
    check_output_files,
    get_expected_output_files,
    get_genome_count,
    parse_cactus_version
)


class CactusPangenomeRunner:
    """Cactus泛基因组分析运行器|Cactus Pangenome Analysis Runner"""

    def __init__(self, config: CactusConfig, logger: Optional[CactusLogger] = None):
        """
        初始化运行器|Initialize runner

        Args:
            config: Cactus配置对象|Cactus configuration object
            logger: 日志对象|Logger object (可选|optional)
        """
        self.config = config
        self.logger_obj = logger
        self.logger = logger.get_logger() if logger else None
        self.absolute_seqfile = None  # 存储使用绝对路径的seqfile|Store seqfile with absolute paths

    def run(self) -> bool:
        """
        运行Cactus泛基因组分析|Run Cactus pangenome analysis

        Returns:
            是否成功|Whether successful
        """
        try:
            # 1. 设置日志|Setup logging
            if self.logger_obj is None:
                log_file = Path(self.config.output_dir) / f"{self.config.out_name}.log"
                self.logger_obj = CactusLogger(log_file, self.config.log_level)
                self.logger = self.logger_obj.get_logger()

            # 2. 打印开始信息|Print start information
            self._print_header()

            # 3. 验证配置|Validate configuration
            if not self._validate_environment():
                return False

            # 4. 检查是否已完成|Check if already completed
            if self._check_completion():
                self.logger.info("")
                self.logger.info("分析已完成，跳过|Analysis already completed, skipping")
                return True

            # 5. 构建并执行命令|Build and execute command
            success = self._run_cactus()

            # 6. 清理或保留jobstore|Cleanup or keep jobstore
            self._handle_jobstore(success)

            # 7. 打印结束信息|Print end information
            self._print_footer(success)

            return success

        except Exception as e:
            if self.logger:
                self.logger.error(f"分析失败|Analysis failed: {e}")
            return False

    def _print_header(self):
        """打印头部信息|Print header information"""
        if not self.logger:
            return

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("Cactus泛基因组分析|Cactus Pangenome Analysis")
        self.logger.info("=" * 60)
        self.logger.info(f"版本|Version: 1.0.0")
        self.logger.info(f"时间|Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  序列文件|Sequence file: {self.config.seqfile}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  参考基因组|Reference: {self.config.reference}")
        self.logger.info(f"  输出名称|Output name: {self.config.out_name}")
        self.logger.info(f"  输出格式|Output formats: {', '.join(self.config.output_formats)}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  最大内存|Max memory: {self.config.max_memory}")
        self.logger.info(f"  工作目录|Work directory: {self.config.work_dir}")
        self.logger.info(f"  Singularity|Singularity: {self.config.singularity_path}")
        self.logger.info(f"  Cactus SIF|Cactus SIF: {self.config.cactus_sif}")
        self.logger.info("")

    def _validate_environment(self) -> bool:
        """
        验证运行环境|Validate runtime environment

        Returns:
            是否验证通过|Whether validation passed
        """
        if not self.logger:
            return False

        self.logger.info("验证运行环境|Validating runtime environment")
        self.logger.info("-" * 60)

        # 1. 验证Singularity|Validate Singularity
        if not validate_singularity(self.config.singularity_path, self.logger):
            return False

        # 2. 验证Cactus SIF|Validate Cactus SIF
        if not validate_cactus_sif(self.config.cactus_sif, self.logger):
            return False

        # 3. 解析Cactus版本|Parse Cactus version
        parse_cactus_version(self.config.cactus_sif, self.config.singularity_path, self.logger)

        # 4. 创建使用绝对路径的seqfile|Create seqfile with absolute paths
        self.absolute_seqfile = create_absolute_path_seqfile(
            self.config.seqfile,
            self.config.output_dir,
            self.logger
        )

        # 5. 检查基因组文件|Check genome files
        all_exist, genome_files = check_genome_files(self.config.seqfile, self.logger)
        if not all_exist:
            self.logger.error("基因组文件检查失败|Genome file check failed")
            return False

        genome_count = get_genome_count(self.config.seqfile)
        self.logger.info(f" 基因组数量|Genome count: {genome_count}")

        self.logger.info("")
        return True

    def _check_completion(self) -> bool:
        """
        检查分析是否已完成|Check if analysis is already completed

        Returns:
            是否已完成|Whether already completed
        """
        if not self.logger:
            return False

        self.logger.info("检查输出文件|Checking output files")

        # 检查输出文件是否存在|Check if output files exist
        if check_output_files(self.config.output_dir, self.config.out_name, self.config.output_formats):
            self.logger.info("  所有输出文件已存在|All output files already exist")
            return True
        else:
            self.logger.info("  输出文件不完整，需要重新运行|Output files incomplete, need to rerun")
            return False

    def _run_cactus(self) -> bool:
        """
        运行Cactus分析|Run Cactus analysis

        Returns:
            是否成功|Whether successful
        """
        if not self.logger:
            return False

        self.logger.info("")
        self.logger.info("开始Cactus分析|Starting Cactus analysis")
        self.logger.info("-" * 60)

        # 1. 构建命令|Build command
        jobstore_path = self.config.get_jobstore_path()
        bind_args = self.config.get_bind_args()

        command = build_singularity_command(
            singularity_path=self.config.singularity_path,
            cactus_sif=self.config.cactus_sif,
            jobstore=jobstore_path,
            seqfile=self.absolute_seqfile if self.absolute_seqfile else self.config.seqfile,
            output_dir=self.config.output_dir,
            out_name=self.config.out_name,
            reference=self.config.reference,
            output_formats=self.config.output_formats,
            threads=self.config.threads,
            max_memory=self.config.max_memory,
            bind_args=bind_args,
            work_dir=self.config.work_dir,
            logger=self.logger
        )

        # 2. 打印命令|Print command
        self.logger.info("")
        self.logger.info(f"执行命令|Executing command:")
        self.logger.info(f"  {command}")
        self.logger.info("")

        # 3. 执行命令|Execute command
        try:
            # 使用subprocess运行命令|Run command using subprocess
            result = subprocess.run(
                command,
                shell=True,
                cwd=self.config.output_dir,
                timeout=None,  # Cactus可能需要很长时间|Cactus may take a long time
            )

            # 4. 检查返回码|Check return code
            if result.returncode == 0:
                self.logger.info("")
                self.logger.info("Cactus分析成功完成|Cactus analysis completed successfully")
                return True
            else:
                self.logger.error("")
                self.logger.error(f"Cactus分析失败，返回码|Cactus analysis failed, return code: {result.returncode}")

                # 如果有stderr输出，显示最后几行|If stderr output, show last few lines
                if result.stderr:
                    stderr_lines = result.stderr.decode('utf-8', errors='ignore').split('\n') if isinstance(result.stderr, bytes) else result.stderr.split('\n')
                    self.logger.error("错误信息|Error message (last 10 lines):")
                    for line in stderr_lines[-10:]:
                        if line.strip():
                            self.logger.error(f"  {line}")

                return False

        except subprocess.TimeoutExpired:
            self.logger.error("命令执行超时|Command execution timeout")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {e}")
            return False

    def _handle_jobstore(self, success: bool):
        """
        处理jobstore清理|Handle jobstore cleanup

        Args:
            success: 分析是否成功|Whether analysis succeeded
        """
        if not self.logger:
            return

        jobstore_path = self.config.get_jobstore_path()

        # 如果配置了清理且分析成功|If configured to cleanup and analysis succeeded
        if self.config.cleanup and success:
            self.logger.info("")
            self.logger.info(f"清理jobstore|Cleaning up jobstore: {jobstore_path}")
            try:
                if Path(jobstore_path).exists():
                    shutil.rmtree(jobstore_path)
                    self.logger.info("  jobstore已删除|jobstore deleted")
                else:
                    self.logger.info("  jobstore不存在，跳过|jobstore does not exist, skipping")
            except Exception as e:
                self.logger.warning(f"  jobstore清理失败|jobstore cleanup failed: {e}")
        else:
            self.logger.info("")
            self.logger.info(f"保留jobstore|Keeping jobstore: {jobstore_path}")

    def _print_footer(self, success: bool):
        """
        打印尾部信息|Print footer information

        Args:
            success: 分析是否成功|Whether analysis succeeded
        """
        if not self.logger:
            return

        self.logger.info("")
        self.logger.info("=" * 60)

        if success:
            self.logger.info("分析成功完成|Analysis completed successfully")
            self.logger.info("")
            self.logger.info("输出文件|Output files:")

            # 列出生成的文件|List generated files
            output_path = Path(self.config.output_dir)
            expected_files = get_expected_output_files(self.config.out_name, self.config.output_formats)

            for filename in expected_files:
                file_path = output_path / filename
                if file_path.exists():
                    file_size = file_path.stat().st_size
                    file_size_mb = file_size / (1024 * 1024)
                    self.logger.info(f"  [OK] {filename} ({file_size_mb:.2f} MB)")
                else:
                    self.logger.warning(f"  [MISSING] {filename} (未找到|not found)")

            # 日志文件|Log file
            log_file = output_path / f"{self.config.out_name}.log"
            if log_file.exists():
                self.logger.info(f"  [OK] {log_file.name} (日志|log)")
        else:
            self.logger.error("分析失败|Analysis failed")
            self.logger.error("请检查日志文件获取更多信息|Please check log file for more details")

        self.logger.info("=" * 60)
        self.logger.info("")
