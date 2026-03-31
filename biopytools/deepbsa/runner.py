"""
DeepBSA执行器模块|DeepBSA Runner Module
"""

import subprocess
from pathlib import Path
from typing import Dict, List

from .config import DeepBSAConfig
from .utils import clean_vcf_comments, run_deepbsa_method, check_parallel_results
from .merge_results import DeepBSAMerger


class DeepBSARunner:
    """DeepBSA执行器类|DeepBSA Runner Class"""

    def __init__(self, config: DeepBSAConfig, logger):
        """初始化执行器|Initialize runner

        Args:
            config: DeepBSA配置对象|DeepBSA configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger
        self.clean_file = None

    def preprocess_input(self) -> Path:
        """预处理输入文件|Preprocess input file

        Returns:
            Path: 处理后的文件路径|Processed file path
        """
        input_file = self.config.input_path

        # 检测文件类型|Detect file type
        # DeepBSA用简单的split('.')解析文件名，无法处理复杂的文件名
        # 例如：variation.filtered.snp.biallelic.vcf 会被解析为 file_type="filtered"
        # DeepBSA uses simple split('.') to parse filename, cannot handle complex filenames
        # e.g., variation.filtered.snp.biallelic.vcf will be parsed as file_type="filtered"

        file_ext = input_file.suffix.lower()

        # 如果是VCF文件（检查扩展名）|If VCF file (check extension)
        if file_ext in ['.vcf', '.vcf.gz', '.vcf.bz2']:
            # DeepBSA只能识别简单文件名（如name.vcf），无法识别复杂文件名
            # DeepBSA can only recognize simple filenames (e.g., name.vcf), not complex ones
            # 如果文件名包含多个点，创建符号链接|If filename contains multiple dots, create symlink
            if '.' in input_file.stem:
                import tempfile
                import shutil

                # 在输出目录创建临时符号链接|Create temporary symlink in output directory
                temp_link_name = self.config.output_path / f"deepbsa_input{file_ext}"

                self.logger.info(f"检测到复杂文件名|Detected complex filename: {input_file.name}")
                self.logger.info(f"创建临时符号链接|Creating temporary symlink: {temp_link_name}")

                # 删除已存在的符号链接|Remove existing symlink if any
                if temp_link_name.exists() or temp_link_name.is_symlink():
                    temp_link_name.unlink()

                # 创建符号链接|Create symlink
                temp_link_name.symlink_to(input_file.absolute())

                self.clean_file = temp_link_name
                return temp_link_name
            else:
                # 简单文件名，直接使用|Simple filename, use directly
                return input_file
        else:
            # 非VCF文件，直接使用|Non-VCF file, use directly
            return input_file

    def run_serial(self, input_file: Path) -> Dict[str, bool]:
        """串行运行所有方法|Run all methods serially

        Args:
            input_file: 输入文件|Input file

        Returns:
            dict: 每个方法的运行结果|Run result for each method
        """
        results = {}

        # 创建each目录存放各方法结果|Create 'each' directory for individual method results
        each_dir = self.config.output_path / "each"
        each_dir.mkdir(exist_ok=True)
        self.logger.info(f"创建各方法结果目录|Created directory for individual results: {each_dir}")

        for method in self.config.methods_list:
            # 为每个方法创建独立的工作目录|Create independent working directory for each method
            method_work_dir = each_dir / method
            method_work_dir.mkdir(exist_ok=True)

            # 创建Models符号链接（DL方法需要预训练模型）
            # Create Models symlink (DL method requires pre-trained models)
            models_src = self.config.deepbsa_script.parent / "Models"
            models_link = method_work_dir / "Models"
            if models_src.exists() and not models_link.exists():
                models_link.symlink_to(models_src)
                self.logger.debug(f"创建Models符号链接|Created Models symlink: {models_link}")

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info(f"运行方法|Running method: {method}")
            self.logger.info(f"工作目录|Working directory: {method_work_dir}")
            self.logger.info("=" * 60)

            success = run_deepbsa_method(
                method=method,
                input_file=input_file,
                deepbsa_script=self.config.deepbsa_script,
                output_dir=method_work_dir,
                conda_env_name=self.config.conda_env_name,
                logger=self.logger,
                parallel=False
            )

            results[method] = success

        return results

    def run_parallel(self, input_file: Path) -> Dict[str, bool]:
        """并行运行所有方法|Run all methods in parallel

        Args:
            input_file: 输入文件|Input file

        Returns:
            dict: 每个方法的运行结果|Run result for each method
        """
        self.logger.info("")
        self.logger.info("并行运行所有方法|Running all methods in parallel...")

        # 创建each目录存放各方法结果|Create 'each' directory for individual method results
        each_dir = self.config.output_path / "each"
        each_dir.mkdir(exist_ok=True)
        self.logger.info(f"创建各方法结果目录|Created directory for individual results: {each_dir}")

        processes = []

        for method in self.config.methods_list:
            # 为每个方法创建独立的工作目录和日志文件
            # Create independent working directory and log file for each method
            method_work_dir = each_dir / method
            method_work_dir.mkdir(exist_ok=True)

            # 创建Models符号链接（DL方法需要预训练模型）
            # Create Models symlink (DL method requires pre-trained models)
            models_src = self.config.deepbsa_script.parent / "Models"
            models_link = method_work_dir / "Models"
            if models_src.exists() and not models_link.exists():
                models_link.symlink_to(models_src)
                self.logger.debug(f"创建Models符号链接|Created Models symlink: {models_link}")

            log_file = method_work_dir / f"{method}.log"
            self.logger.info(f"启动|Starting {method} (工作目录|work_dir: {method_work_dir})")

            cmd = [
                "conda", "run", "-n", self.config.conda_env_name, "--no-capture-output",
                "python3",
                str(self.config.deepbsa_script),
                "--i", str(input_file),
                "--m", method,
                "--p", "0"  # 禁用pretreatment，使用已过滤的数据|Disable pretreatment, use pre-filtered data
            ]

            # 启动后台进程，在方法专属目录中运行|Start background process in method-specific directory
            with open(log_file, 'w') as f:
                process = subprocess.Popen(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=str(method_work_dir)  # 设置工作目录为方法专属目录|Set working directory to method-specific directory
                )
                processes.append((method, process))

        # 等待所有进程完成|Wait for all processes to complete
        self.logger.info("")
        self.logger.info("等待所有方法完成|Waiting for all methods to complete...")

        results = {}
        for method, process in processes:
            return_code = process.wait()
            results[method] = (return_code == 0)

        self.logger.info("")
        self.logger.info("所有方法运行完成|All methods completed!")

        # 检查结果|Check results
        check_parallel_results(self.config.output_path, self.config.methods_list, self.logger)

        return results

    def run(self) -> Dict[str, bool]:
        """运行DeepBSA分析|Run DeepBSA analysis

        Returns:
            dict: 每个方法的运行结果|Run result for each method
        """
        # 显示配置信息|Display configuration info
        self.logger.info("=" * 60)
        self.logger.info("DeepBSA 批量分析|DeepBSA Batch Analysis")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input file: {self.config.input_path}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_path}")
        self.logger.info(f"运行方法|Methods: {', '.join(self.config.methods_list)}")
        self.logger.info(f"并行模式|Parallel mode: {'是|Yes' if self.config.parallel else '否|No'}")
        self.logger.info(f"工作目录|Working directories: 每个方法独立目录|Each method in independent directory")
        self.logger.info("=" * 60)
        self.logger.info("")

        # 预处理输入文件（可能创建符号链接）|Preprocess input file (may create symlink)
        actual_input = self.preprocess_input()

        if self.clean_file:
            self.logger.info("")

        # 运行分析|Run analysis
        if self.config.parallel:
            results = self.run_parallel(actual_input)
        else:
            results = self.run_serial(actual_input)

        # 清理临时文件|Clean temporary files
        self.cleanup()

        # 合并结果|Merge results
        self.logger.info("")
        merger = DeepBSAMerger(self.config.output_path, self.logger)
        merger.run(self.config.methods_list)

        return results

    def cleanup(self):
        """清理临时文件|Clean temporary files"""
        if self.clean_file and (self.clean_file.exists() or self.clean_file.is_symlink()):
            self.logger.info("")
            self.logger.info(f"清理临时符号链接|Cleaning temporary symlink: {self.clean_file}")
            self.clean_file.unlink()
            self.logger.info("临时符号链接已删除|Temporary symlink deleted")

    def print_summary(self, results: Dict[str, bool]):
        """打印运行摘要|Print run summary

        Args:
            results: 运行结果|Run results
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("分析完成|Analysis Completed")
        self.logger.info("=" * 60)
        self.logger.info(f"结果保存在|Results saved in: {self.config.output_path}/")

        # 显示输出文件摘要|Display output files summary
        self.logger.info("")
        self.logger.info("输出目录摘要|Output directories summary:")

        # 显示each/目录下的各方法结果|Display individual method results in 'each/' directory
        each_dir = self.config.output_path / "each"
        if each_dir.exists():
            for method_dir in sorted(each_dir.glob("*")):
                if method_dir.is_dir():
                    # 计算目录大小
                    total_size = 0
                    for f in method_dir.rglob('*'):
                        if f.is_file():
                            total_size += f.stat().st_size

                    size_kb = total_size / 1024
                    self.logger.info(f"  each/{method_dir.name}/ ({size_kb:.1f} KB)")

        # 显示合并结果|Display merged results
        merged_dir = self.config.output_path / "merged_results"
        if merged_dir.exists():
            self.logger.info(f"  merged_results/ (合并结果|merged results)")
