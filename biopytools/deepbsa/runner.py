"""
DeepBSA执行器模块|DeepBSA Runner Module
"""

import subprocess
import os
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

    def _is_batch_mode(self) -> bool:
        """检测是否为batch模式|Detect if running in batch mode

        batch模式特征|Batch mode characteristics:
        - 输出路径的最后一部分是方法名|Output path's last part is a method name
        - 只运行一个方法|Only running one method

        Returns:
            bool: 是否为batch模式|Whether batch mode
        """
        if len(self.config.methods_list) != 1:
            return False

        method = self.config.methods_list[0]
        output_dir_name = self.config.output_path.name

        # 检查输出目录名是否是方法名|Check if output dir name is method name
        return output_dir_name == method

    def _calculate_safe_threads(self, num_methods: int) -> int:
        """计算安全的线程数，避免资源冲突|Calculate safe thread count to avoid resource conflicts

        针对64核+300GB集群优化|Optimized for 64-core + 300GB cluster

        Args:
            num_methods: 要并行运行的方法数量|Number of methods to run in parallel

        Returns:
            int: 调整后的线程数|Adjusted thread count
        """
        requested_threads = self.config.threads
        cpu_count = os.cpu_count() or 1

        # 如果用户指定了线程数（非0）
        if requested_threads != 0:
            # 计算总进程数
            total_processes = num_methods * requested_threads

            # 针对64核+300GB集群的激进阈值
            # 警告阈值: CPU核心数的90%（接近满载）
            warning_threshold = int(cpu_count * 0.90)
            # 严重警告阈值: CPU核心数的1.5倍（允许更多超载）
            critical_threshold = int(cpu_count * 1.5)

            if total_processes > critical_threshold:
                self.logger.warning("")
                self.logger.warning("=" * 60)
                self.logger.warning(f"资源警告|Resource Warning: 总进程数|Total processes = {total_processes}")
                self.logger.warning(f"  方法数|Methods: {num_methods}")
                self.logger.warning(f"  每个方法线程|Threads per method: {requested_threads}")
                self.logger.warning(f"  CPU核心数|CPU cores: {cpu_count}")
                self.logger.warning(f"  内存|Memory: 300GB (估算|Estimated)")
                self.logger.warning("")
                self.logger.warning("建议|Recommendation:")
                self.logger.warning(f"  - 并行运行时，建议每个方法使用 4-8 个线程|In parallel mode, recommend 4-8 threads per method")
                self.logger.warning(f"  - 使用 --threads 8 (总进程|total processes {num_methods * 8})")
                self.logger.warning(f"  - 或使用 --no-parallel 串行运行，然后使用 --threads {min(32, cpu_count)}|Or use --no-parallel for serial execution with --threads {min(32, cpu_count)}")
                self.logger.warning("=" * 60)
                self.logger.warning("")

                # 自动限制为 8 个线程（适合64核集群）
                adjusted_threads = 8
                self.logger.warning(f"自动调整为安全值|Auto-adjusted to safe value: {adjusted_threads} threads per method")
                self.logger.warning(f"总进程数|Total processes: {num_methods} × {adjusted_threads} = {num_methods * adjusted_threads}")
                self.logger.warning("如需使用更多线程，请使用 --no-parallel 串行运行|To use more threads, run serially with --no-parallel")
                self.logger.warning("")
                return adjusted_threads

            elif total_processes > warning_threshold:
                self.logger.warning("")
                self.logger.warning(f"注意|Note: 总进程数|Total processes = {total_processes} (接近CPU满载|close to full capacity)")
                self.logger.warning(f"  300GB内存应该足够|Memory should be sufficient")
                self.logger.warning(f"  建议: 并行模式使用 --threads 6-8|Recommendation: use --threads 6-8 in parallel mode")
                self.logger.warning("")

        # 如果 threads=0（自动检测），在并行模式下设置为6（适合64核集群）
        if requested_threads == 0 and num_methods > 1:
            self.logger.info("")
            self.logger.info("并行模式|Parallel mode: 自动设置每个方法为6线程|Auto-set each method to 6 threads")
            self.logger.info(f"总进程数|Total processes: {num_methods} × 6 = {num_methods * 6}")
            self.logger.info(f"针对64核+300GB集群优化|Optimized for 64-core + 300GB cluster")
            self.logger.info("")
            return 6

        return requested_threads

    def check_method_completed(self, method: str) -> bool:
        """检查方法是否已完成|Check if method is already completed

        Args:
            method: 方法名|Method name

        Returns:
            bool: True表示已完成|True indicates completed
        """
        # 根据模式确定结果目录|Determine result directory based on mode
        if self._is_batch_mode():
            # batch模式: 结果在输出目录下的Results/目录
            # Batch mode: results in Results/ directory under output
            method_dir = self.config.output_path
        else:
            # 普通run模式: 结果在methods/{method}/目录
            # Normal run mode: results in methods/{method}/ directory
            method_dir = self.config.output_path / "methods" / method

        if not method_dir.exists():
            return False

        # 检查是否存在values.txt文件（这是方法完成的标志）
        # Check if values.txt file exists (this indicates method completion)
        values_files = list(method_dir.glob("Results/**/* values.txt"))

        if values_files:
            self.logger.info(f"已完成，跳过|Completed, skipping: {method}")
            return True

        return False

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

        # 检测是否为batch模式|Detect if batch mode
        is_batch_mode = self._is_batch_mode()

        if is_batch_mode:
            # batch模式: 直接在输出目录下运行，不创建methods子目录
            # Batch mode: Run directly in output directory, don't create methods subdirectory
            methods_base_dir = self.config.output_path
            self.logger.info(f"Batch模式|Batch mode detected: 输出目录|Output dir = {methods_base_dir}")
        else:
            # 普通run模式: 创建methods子目录存放各方法
            # Normal run mode: Create methods subdirectory for each method
            methods_base_dir = self.config.output_path / "methods"
            methods_base_dir.mkdir(exist_ok=True)
            self.logger.info(f"创建方法目录|Created methods directory: {methods_base_dir}")

        # 检查哪些方法需要运行|Check which methods need to run
        if not self.config.force:
            self.logger.info("")
            self.logger.info("检查已完成的方法|Checking completed methods...")

        # 确定Models源路径|Determine Models source path
        # 对于内置版本，Models在deepbsa_builtin/Models
        # For builtin version, Models is in deepbsa_builtin/Models
        models_src = self.config.deepbsa_script.parent / "Models"

        for method in self.config.methods_list:
            # 检查方法是否已完成|Check if method is already completed
            if not self.config.force and self.check_method_completed(method):
                results[method] = True
                continue

            # 为每个方法创建独立的工作目录|Create independent working directory for each method
            if is_batch_mode:
                # batch模式: 工作目录就是输出目录本身
                # Batch mode: work directory is the output directory itself
                method_work_dir = methods_base_dir
                self.logger.debug(f"Batch模式工作目录|Batch mode work dir: {method_work_dir}")
            else:
                # 普通run模式: 在methods_base_dir下创建子目录
                # Normal run mode: create subdirectory under methods_base_dir
                method_work_dir = methods_base_dir / method
                method_work_dir.mkdir(exist_ok=True)

            # 创建Models符号链接（DL方法需要预训练模型）
            # Create Models symlink (DL method requires pre-trained models)
            models_link = method_work_dir / "Models"
            if models_src.exists() and not models_link.exists():
                models_link.symlink_to(models_src)
                self.logger.debug(f"创建Models符号链接|Created Models symlink: {models_link}")

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info(f"运行方法|Running method: {method}")
            self.logger.info(f"工作目录|Working directory: {method_work_dir}")
            if self.config.threads > 0:
                self.logger.info(f"线程数|Threads: {self.config.threads}")
            else:
                self.logger.info(f"线程数|Threads: auto (自动检测|auto-detect)")
            self.logger.info("=" * 60)

            success = run_deepbsa_method(
                method=method,
                input_file=input_file,
                deepbsa_script=self.config.deepbsa_script,
                output_dir=method_work_dir,
                conda_env_name=self.config.conda_env_name,
                logger=self.logger,
                threads=self.config.threads,
                parallel=False,
                smooth_func=self.config.smooth_func
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

        # 计算安全的线程数|Calculate safe thread count
        num_methods = len(self.config.methods_list)
        safe_threads = self._calculate_safe_threads(num_methods)

        # 创建各方法独立目录|Create individual directories for each method
        # 不使用each中间层，避免路径嵌套|Don't use 'each' intermediate layer to avoid path nesting
        methods_dir = self.config.output_path / "methods"
        methods_dir.mkdir(exist_ok=True)
        self.logger.info(f"创建方法目录|Created methods directory: {methods_dir}")

        # 检查哪些方法需要运行|Check which methods need to run
        methods_to_run = []
        methods_skipped = []

        if not self.config.force:
            self.logger.info("")
            self.logger.info("检查已完成的方法|Checking completed methods...")
            for method in self.config.methods_list:
                if self.check_method_completed(method):
                    methods_skipped.append(method)
                else:
                    methods_to_run.append(method)
        else:
            methods_to_run = self.config.methods_list[:]

        # 显示跳过的方法|Show skipped methods
        if methods_skipped:
            self.logger.info("")
            self.logger.info(f"跳过已完成的方法|Skipping completed methods: {', '.join(methods_skipped)}")
            self.logger.info(f"将运行的方法|Methods to run: {', '.join(methods_to_run) if methods_to_run else '无|None'}")

        # 如果所有方法都已完成|If all methods are completed
        if not methods_to_run:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("所有方法均已完成，无需运行|All methods completed, nothing to run")
            self.logger.info("=" * 60)
            # 返回所有方法为成功|Return all methods as successful
            return {method: True for method in self.config.methods_list}

        # 确定Models源路径|Determine Models source path
        # 对于内置版本，Models在deepbsa_builtin/Models
        # For builtin version, Models is in deepbsa_builtin/Models
        models_src = self.config.deepbsa_script.parent / "Models"

        processes = []

        for method in methods_to_run:
            # 为每个方法创建独立的工作目录和日志文件
            # Create independent working directory and log file for each method
            # 直接在 methods/ 目录下，不嵌套|Directly under methods/, no nesting
            method_work_dir = methods_dir / method
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
            if safe_threads != self.config.threads:
                self.logger.warning(f"  线程数已调整|Threads adjusted: {self.config.threads} → {safe_threads}")
            else:
                if safe_threads > 0:
                    self.logger.info(f"  线程数|Threads: {safe_threads}")
                else:
                    self.logger.info(f"  线程数|Threads: auto (自动检测|auto-detect)")

            cmd = [
                "conda", "run", "-n", self.config.conda_env_name, "--no-capture-output",
                "python3",
                str(self.config.deepbsa_script),
                "--i", str(input_file),
                "--m", method,
                "--p", "0",  # 禁用pretreatment，使用已过滤的数据|Disable pretreatment, use pre-filtered data
                "--s", self.config.smooth_func,  # 平滑函数|Smooth function
                "--threads", str(safe_threads)  # 使用调整后的线程数|Use adjusted thread count
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
        if processes:
            self.logger.info("")
            self.logger.info("等待所有方法完成|Waiting for all methods to complete...")

        results = {}
        for method, process in processes:
            return_code = process.wait()
            results[method] = (return_code == 0)

        # 跳过的方法标记为成功|Mark skipped methods as successful
        for method in methods_skipped:
            results[method] = True

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

        # 显示DeepBSA版本信息|Display DeepBSA version info
        if getattr(self.config, '_using_builtin', False):
            self.logger.info(f"DeepBSA版本|Version: 内置版本|Builtin")
            self.logger.info(f"DeepBSA路径|Path: {self.config.deepbsa_script}")
        else:
            self.logger.info(f"DeepBSA版本|Version: 外部版本|External")
            self.logger.info(f"DeepBSA路径|Path: {self.config.deepbsa_script}")

        self.logger.info(f"线程数|Threads: {self.config.threads if self.config.threads > 0 else 'auto'}")
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
        if not self.config.skip_merge:
            # 检查是否需要跳过合并（如果合并结果已存在）
            merged_csv = self.config.output_path / "merged_results" / "merged_results.csv"

            if not self.config.force and merged_csv.exists():
                self.logger.info("")
                self.logger.info("=" * 60)
                self.logger.info("合并结果已存在，跳过合并步骤|Merged results already exist, skipping merge step")
                self.logger.info(f"如需重新合并，请使用 --force 参数|Use --force to re-merge")
                self.logger.info(f"现有合并结果|Existing merged results: {merged_csv}")
                self.logger.info("=" * 60)
            else:
                self.logger.info("")
                merger = DeepBSAMerger(self.config.output_path, self.logger)
                merger.run(self.config.methods_list)
        else:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("跳过合并步骤|Skipping merge step (--skip-merge)")
            self.logger.info("=" * 60)

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

        # 根据模式显示不同的目录结构|Display different directory structure based on mode
        if self._is_batch_mode():
            # batch模式: 直接显示输出目录下的Results/
            # Batch mode: directly show Results/ under output directory
            results_dir = self.config.output_path / "Results"
            if results_dir.exists():
                # 计算目录大小
                total_size = 0
                for f in results_dir.rglob('*'):
                    if f.is_file():
                        total_size += f.stat().st_size
                size_kb = total_size / 1024
                self.logger.info(f"  Results/ ({size_kb:.1f} KB)")
        else:
            # 普通run模式: 显示methods/目录下的各方法结果
            # Normal run mode: display individual method results in 'methods/' directory
            methods_dir = self.config.output_path / "methods"
            if methods_dir.exists():
                for method_dir in sorted(methods_dir.glob("*")):
                    if method_dir.is_dir():
                        # 计算目录大小
                        total_size = 0
                        for f in method_dir.rglob('*'):
                            if f.is_file():
                                total_size += f.stat().st_size

                        size_kb = total_size / 1024
                        self.logger.info(f"  methods/{method_dir.name}/ ({size_kb:.1f} KB)")

        # 显示合并结果|Display merged results
        merged_dir = self.config.output_path / "merged_results"
        if merged_dir.exists():
            self.logger.info(f"  merged_results/ (合并结果|merged results)")
