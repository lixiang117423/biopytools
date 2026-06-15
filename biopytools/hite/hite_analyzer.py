"""
HiTE 单基因组分析器|HiTE Single-genome Analyzer

提供HiTE单基因组转座子检测与注释的核心功能
Provides core functionality for HiTE single-genome transposon detection and annotation
"""

import sys
from .config import HiteConfig
from .logger import HiteLoggerManager
from .utils import SingularityContainerManager
from .results import HiteResultsProcessor


class HiteAnalyzer:
    """
    HiTE 单基因组分析器|HiTE Single-genome Analyzer

    通过Singularity容器调用HiTE进行转座子检测与注释
    Perform transposon detection and annotation by calling HiTE through Singularity container
    """

    def __init__(self, **kwargs):
        """
        初始化HiTE分析器|Initialize HiTE analyzer

        Args:
            **kwargs: 配置参数|Configuration parameters

        Raises:
            ValueError: 如果配置验证失败|If configuration validation fails
        """
        # 初始化配置|Initialize configuration
        self.config = HiteConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = HiteLoggerManager(
            self.config.output_dir,
            log_prefix="hite",
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化容器管理器|Initialize container manager
        self.container_manager = SingularityContainerManager(
            self.config,
            self.logger
        )

        # 初始化结果处理器|Initialize results processor
        self.results_processor = HiteResultsProcessor(
            self.config,
            self.logger
        )

        self._log_initialization_info()

    def _log_initialization_info(self):
        """
        记录初始化信息|Log initialization information
        """
        self.logger.info("="*80)
        self.logger.info(
            "HiTE 转座子检测与注释分析器初始化 | "
            "HiTE Transposon Detection and Annotation Analyzer Initialized"
        )
        self.logger.info("="*80)
        self.logger.info(f"HiTE版本 | HiTE version: 3.3.3 (Singularity)")
        self.logger.info(f"基因组文件 | Genome file: {self.config.genome}")
        self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数 | Threads: {self.config.threads}")
        self.logger.info(f"植物基因组 | Plant genome: {self.config.plant}")
        self.logger.info(f"基因组注释 | Annotate genome: {self.config.annotate}")
        self.logger.info(f"断点续跑 | Recovery mode: {self.config.recover}")
        self.logger.info(f"结构域预测 | Domain prediction: {self.config.domain}")
        self.logger.info(f"TE类型 | TE type: {self.config.te_type}")
        self.logger.info(f"Singularity | Singularity: {self.config.singularity_cmd}")
        self.logger.info(f"SIF镜像 | SIF image: {self.config.sif_file}")
        self.logger.info("="*80)

    def _build_hite_arguments(self):
        """
        构建HiTE参数列表|Build HiTE argument list

        将配置对象转换为HiTE命令行参数
        Convert configuration object to HiTE command line arguments

        使用容器内输出路径，避免挂载冲突
        Use container-internal output path to avoid mount conflicts

        Returns:
            list: HiTE参数列表|HiTE argument list
        """
        args = []

        # 使用容器内的输出路径，避免挂载冲突
        # Use container-internal output path to avoid mount conflicts
        container_output_dir = '/tmp/hite_output'

        # 必需参数|Required parameters
        args.extend(['--genome', self.config.genome])
        args.extend(['--thread', str(self.config.threads)])
        args.extend(['--out_dir', container_output_dir])

        # 工作目录（使用容器内临时目录，避免挂载冲突）
        # Work directory (use container's temp dir to avoid mount conflicts)
        args.extend(['--work_dir', '/tmp/hite_work'])

        # 布尔参数|Boolean parameters
        args.extend(['--plant', '1' if self.config.plant else '0'])
        args.extend(['--annotate', '1' if self.config.annotate else '0'])
        args.extend(['--recover', '1' if self.config.recover else '0'])
        args.extend(['--domain', '1' if self.config.domain else '0'])

        # 高级参数（仅当非默认值时添加）|Advanced parameters (add only if non-default)
        if self.config.chunk_size != 400:
            args.extend(['--chunk_size', str(self.config.chunk_size)])

        if self.config.miu != 1.3e-8:
            args.extend(['--miu', str(self.config.miu)])

        if self.config.min_te_len != 80:
            args.extend(['--min_TE_len', str(self.config.min_te_len)])

        if self.config.te_type != 'all':
            args.extend(['--te_type', self.config.te_type])

        # 嵌套TE处理|Nested TE handling
        args.extend(['--remove_nested', '1' if self.config.remove_nested else '0'])

        return args

    def _copy_results_from_container(self):
        """
        从容器复制结果文件到主机输出目录
        Copy result files from container to host output directory

        使用容器管理器的复制方法
        Use container manager's copy method
        """
        container_output_dir = '/tmp/hite_output'
        host_output_dir = self.config.output_dir

        return self.container_manager.copy_results_from_container(
            container_output_dir,
            host_output_dir
        )

    def run_analysis(self) -> bool:
        """
        运行HiTE分析|Run HiTE analysis

        执行完整的HiTE转座子检测与注释流程
        Execute complete HiTE transposon detection and annotation pipeline

        Returns:
            bool: 分析成功返回True|Returns True if analysis succeeds

        Raises:
            Exception: 如果分析失败|If analysis fails
        """
        start_time = None

        try:
            import time
            start_time = time.time()

            self.logger.info("")
            self.logger.info("="*80)
            self.logger.info(
                "开始 HiTE 转座子检测与注释分析 | "
                "Starting HiTE Transposon Detection and Annotation Analysis"
            )
            self.logger.info("="*80)
            self.logger.info("")

            # 步骤1: 构建HiTE参数|Step 1: Build HiTE arguments
            self.logger.info("步骤1 | Step 1: 构建HiTE参数 | Building HiTE arguments")
            hite_args = self._build_hite_arguments()
            self.logger.info(f"HiTE参数 | HiTE arguments: {' '.join(hite_args)}")

            # 步骤2: 运行HiTE|Step 2: Run HiTE
            self.logger.info("")
            self.logger.info("步骤2 | Step 2: 运行HiTE | Running HiTE")
            self.logger.info("")

            internal_cmd = ['python', '/HiTE/main.py'] + hite_args
            self.container_manager.run_container_command(
                self.config,
                internal_cmd,
                "HiTE"
            )

            # 步骤3: 从容器复制结果|Step 3: Copy results from container
            self.logger.info("")
            self.logger.info("步骤3 | Step 3: 从容器复制结果 | Copying results from container")
            self.logger.info("")

            self._copy_results_from_container()

            # 步骤4: 处理结果|Step 4: Process results
            self.logger.info("")
            self.logger.info("步骤4 | Step 4: 处理HiTE结果 | Processing HiTE results")
            self.logger.info("")

            self.results_processor.process_results()

            # 步骤5: 生成汇总|Step 5: Generate summary
            self.logger.info("")
            self.logger.info("步骤5 | Step 5: 生成结果汇总 | Generating results summary")
            self.logger.info("")

            self.results_processor.generate_summary()
            self.results_processor.print_summary()

            # 分析完成|Analysis completed
            if start_time:
                elapsed_time = time.time() - start_time
                hours = int(elapsed_time // 3600)
                minutes = int((elapsed_time % 3600) // 60)
                seconds = int(elapsed_time % 60)

                self.logger.info("")
                self.logger.info("="*80)
                self.logger.info(
                    "HiTE 分析成功完成 | HiTE Analysis Completed Successfully"
                )
                self.logger.info(
                    f"总用时 | Total time: "
                    f"{hours:02d}:{minutes:02d}:{seconds:02d} "
                    f"({hours}小时|hours {minutes}分钟|minutes {seconds}秒|seconds)"
                )
                self.logger.info("="*80)
                self.logger.info("")

            return True

        except Exception as e:
            self.logger.error("")
            self.logger.error("="*80)
            self.logger.error(f"HiTE 分析失败 | HiTE Analysis Failed: {str(e)}")
            self.logger.error("="*80)
            self.logger.error("")
            raise


def main(**kwargs) -> HiteAnalyzer:
    """
    主函数入口|Main function entry point

    Args:
        **kwargs: 配置参数|Configuration parameters

    Returns:
        HiteAnalyzer: HiTE分析器实例|HiTE analyzer instance
    """
    return HiteAnalyzer(**kwargs)


if __name__ == '__main__':
    """
    命令行入口|Command line entry point

    支持直接运行此模块进行测试
    Support running this module directly for testing
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="HiTE转座子检测与注释|HiTE Transposon Detection and Annotation"
    )
    parser.add_argument(
        '--genome', '-g',
        required=True,
        help='基因组FASTA文件|Genome FASTA file'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='./hite_output',
        help='输出目录|Output directory'
    )
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=40,
        help='线程数|Number of threads'
    )
    parser.add_argument(
        '--plant/--no-plant',
        default=True,
        help='植物基因组|Plant genome'
    )
    parser.add_argument(
        '--annotate',
        action='store_true',
        help='注释基因组|Annotate genome'
    )
    parser.add_argument(
        '--recover',
        action='store_true',
        help='断点续跑|Resume from interruption'
    )

    args = parser.parse_args()

    try:
        analyzer = HiteAnalyzer(
            genome=args.genome,
            output_dir=args.output_dir,
            threads=args.threads,
            plant=args.plant,
            annotate=args.annotate,
            recover=args.recover
        )
        analyzer.run_analysis()
        sys.exit(0)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)
