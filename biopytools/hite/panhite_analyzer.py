"""
panHiTE 群体基因组分析器|panHiTE Pan-genome Analyzer

提供panHiTE泛基因组转座子检测与比较分析的核心功能
Provides core functionality for panHiTE pan-genome transposon detection and comparative analysis
"""

import sys
from .config import PanHiteConfig
from .logger import HiteLoggerManager
from .utils import SingularityContainerManager
from .results import PanHiteResultsProcessor


class PanHiteAnalyzer:
    """
    panHiTE 群体基因组分析器|panHiTE Pan-genome Analyzer

    通过Singularity容器调用panHiTE进行泛基因组转座子检测与比较分析
    Perform pan-genome transposon detection and comparative analysis by calling
    panHiTE through Singularity container
    """

    def __init__(self, **kwargs):
        """
        初始化panHiTE分析器|Initialize panHiTE analyzer

        Args:
            **kwargs: 配置参数|Configuration parameters

        Raises:
            ValueError: 如果配置验证失败|If configuration validation fails
        """
        # 初始化配置|Initialize configuration
        self.config = PanHiteConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = HiteLoggerManager(
            self.config.output_dir,
            log_prefix="panhite",
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化容器管理器|Initialize container manager
        self.container_manager = SingularityContainerManager(
            self.config,
            self.logger
        )

        # 初始化结果处理器|Initialize results processor
        self.results_processor = PanHiteResultsProcessor(
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
            "panHiTE 群体基因组分析器初始化 | "
            "panHiTE Pan-genome Analyzer Initialized"
        )
        self.logger.info("="*80)
        self.logger.info(f"panHiTE版本 | panHiTE version: 0.0.1 (Singularity)")
        self.logger.info(
            f"泛基因组目录 | Pan-genomes directory: {self.config.pan_genomes_dir}"
        )
        self.logger.info(
            f"基因组列表 | Genome list file: {self.config.genome_list}"
        )
        self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数 | Threads: {self.config.threads}")
        self.logger.info(f"跳过分析 | Skip analysis: {self.config.skip_analyze}")
        self.logger.info(f"TE类型 | TE type: {self.config.te_type}")

        if self.config.genes_dir:
            self.logger.info(
                f"基因注释目录 | Genes directory: {self.config.genes_dir}"
            )
        if self.config.rna_dir:
            self.logger.info(
                f"RNA-seq目录 | RNA-seq directory: {self.config.rna_dir}"
            )

        self.logger.info(f"调试模式 | Debug mode: {self.config.debug}")
        self.logger.info(f"断点续跑 | Recovery mode: {self.config.recover}")
        self.logger.info(f"Singularity | Singularity: {self.config.singularity_cmd}")
        self.logger.info(f"SIF镜像 | SIF image: {self.config.sif_file}")
        self.logger.info("="*80)

    def _build_panhite_arguments(self):
        """
        构建panHiTE参数列表|Build panHiTE argument list

        将配置对象转换为panHiTE命令行参数
        Convert configuration object to panHiTE command line arguments

        Returns:
            list: panHiTE参数列表|panHiTE argument list
        """
        args = []

        # 必需参数|Required parameters
        args.extend(['--pan_genomes_dir', self.config.pan_genomes_dir])
        args.extend(['--genome_list', self.config.genome_list])
        args.extend(['--out_dir', '/tmp/panhite_output'])  # 使用容器内路径
        args.extend(['--thread', str(self.config.threads)])

        # 工作目录（使用容器内临时目录，避免挂载冲突）
        # Work directory (use container's temp dir to avoid mount conflicts)
        args.extend(['--work_dir', '/tmp/panhite_work'])

        # 可选路径参数|Optional path parameters
        if self.config.genes_dir:
            args.extend(['--genes_dir', self.config.genes_dir])
        if self.config.rna_dir:
            args.extend(['--RNA_dir', self.config.rna_dir])

        # 处理参数|Processing parameters
        args.extend(['--miu', str(self.config.miu)])
        args.extend(['--te_type', self.config.te_type])
        args.extend(['--skip_analyze', '1' if self.config.skip_analyze else '0'])
        args.extend(['--debug', '1' if self.config.debug else '0'])
        args.extend(['--recover', '1' if self.config.recover else '0'])

        return args

    def _copy_results_from_container(self):
        """
        从容器复制结果文件到主机输出目录
        Copy result files from container to host output directory

        使用容器管理器的复制方法
        Use container manager's copy method
        """
        container_output_dir = '/tmp/panhite_output'
        host_output_dir = self.config.output_dir

        return self.container_manager.copy_results_from_container(
            container_output_dir,
            host_output_dir
        )

    def run_analysis(self) -> bool:
        """
        运行panHiTE分析|Run panHiTE analysis

        执行完整的panHiTE泛基因组转座子检测与比较分析流程
        Execute complete panHiTE pan-genome transposon detection and
        comparative analysis pipeline

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
                "开始 panHiTE 群体基因组转座子分析 | "
                "Starting panHiTE Pan-genome Transposon Analysis"
            )
            self.logger.info("="*80)
            self.logger.info("")

            # 步骤1: 构建panHiTE参数|Step 1: Build panHiTE arguments
            self.logger.info(
                "步骤1 | Step 1: 构建panHiTE参数 | Building panHiTE arguments"
            )
            panhite_args = self._build_panhite_arguments()
            self.logger.info(
                f"panHiTE参数 | panHiTE arguments: {' '.join(panhite_args)}"
            )

            # 步骤2: 运行panHiTE|Step 2: Run panHiTE
            self.logger.info("")
            self.logger.info("步骤2 | Step 2: 运行panHiTE | Running panHiTE")
            self.logger.info("")

            internal_cmd = ['python', '/HiTE/panHiTE.py'] + panhite_args
            self.container_manager.run_container_command(
                self.config,
                internal_cmd,
                "panHiTE"
            )

            # 步骤3: 从容器复制结果|Step 3: Copy results from container
            self.logger.info("")
            self.logger.info(
                "步骤3 | Step 3: 从容器复制结果 | Copying results from container"
            )
            self.logger.info("")

            self._copy_results_from_container()

            # 步骤4: 处理结果|Step 4: Process results
            self.logger.info("")
            self.logger.info(
                "步骤4 | Step 4: 处理panHiTE结果 | Processing panHiTE results"
            )
            self.logger.info("")

            self.results_processor.process_results()

            # 步骤5: 生成汇总|Step 5: Generate summary
            self.logger.info("")
            self.logger.info(
                "步骤5 | Step 5: 生成结果汇总 | Generating results summary"
            )
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
                    "panHiTE 分析成功完成 | panHiTE Analysis Completed Successfully"
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
            self.logger.error(
                f"panHiTE 分析失败 | panHiTE Analysis Failed: {str(e)}"
            )
            self.logger.error("="*80)
            self.logger.error("")
            raise


def main(**kwargs) -> PanHiteAnalyzer:
    """
    主函数入口|Main function entry point

    Args:
        **kwargs: 配置参数|Configuration parameters

    Returns:
        PanHiteAnalyzer: panHiTE分析器实例|panHiTE analyzer instance
    """
    return PanHiteAnalyzer(**kwargs)


if __name__ == '__main__':
    """
    命令行入口|Command line entry point

    支持直接运行此模块进行测试
    Support running this module directly for testing
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="panHiTE群体基因组转座子分析|panHiTE Pan-genome Transposon Analysis"
    )
    parser.add_argument(
        '--pan-genomes-dir', '-p',
        required=True,
        help='泛基因组目录|Pan-genomes directory'
    )
    parser.add_argument(
        '--genome-list', '-l',
        required=True,
        help='基因组列表文件|Genome list file'
    )
    parser.add_argument(
        '--genes-dir',
        help='基因注释目录|Genes annotation directory'
    )
    parser.add_argument(
        '--rna-dir',
        help='RNA-seq数据目录|RNA-seq data directory'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='./panhite_output',
        help='输出目录|Output directory'
    )
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=80,
        help='线程数|Number of threads'
    )
    parser.add_argument(
        '--skip-analyze',
        action='store_true',
        help='跳过分析，仅生成panTE库|Skip analysis, only generate panTE library'
    )

    args = parser.parse_args()

    try:
        analyzer = PanHiteAnalyzer(
            pan_genomes_dir=args.pan_genomes_dir,
            genome_list=args.genome_list,
            genes_dir=args.genes_dir,
            rna_dir=args.rna_dir,
            output_dir=args.output_dir,
            threads=args.threads,
            skip_analyze=args.skip_analyze
        )
        analyzer.run_analysis()
        sys.exit(0)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)
