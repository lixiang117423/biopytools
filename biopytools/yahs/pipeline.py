"""
YaHS流程管理器|YaHS Pipeline Manager

协调Hi-C scaffolding流程的执行，支持断点续传
Coordinates Hi-C scaffolding pipeline execution with checkpoint resume support
"""

import sys
import time
from pathlib import Path
from datetime import datetime

from .config import YaHSConfig
from .utils import YaHSLogger, CommandRunner, generate_software_versions_yml, calculate_runtime
from .steps import YaHSSteps


class YaHSPipeline:
    """
    YaHS流程管理器|YaHS Pipeline Manager

    管理Hi-C scaffolding的完整流程，支持断点续传
    Manages complete Hi-C scaffolding pipeline with checkpoint resume
    """

    def __init__(self, config: YaHSConfig):
        """
        初始化流程管理器|Initialize pipeline manager

        Args:
            config: 配置对象|Configuration object
        """
        self.config = config

        # 初始化日志|Initialize logging
        self.logger_manager = YaHSLogger(
            log_file=str(config.pipeline_log),
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(
            self.logger,
            working_dir=str(config.output_dir)
        )

        # 初始化步骤执行器|Initialize steps executor
        self.steps = YaHSSteps(config, self.logger, self.cmd_runner)

        # 记录开始时间|Record start time
        self.start_time = time.time()

    def run(self) -> int:
        """
        运行完整流程|Run complete pipeline

        Returns:
            退出码|Exit code (0=成功|Success, 1=失败|Failure)
        """
        try:
            self._print_header()

            # 验证配置|Validate configuration
            self.config.validate()

            # 检查环境|Check environment
            self._check_environment()

            # 生成软件版本信息|Generate software version info
            generate_software_versions_yml(
                self.config,
                str(self.config.output_path / '00_pipeline_info' / 'software_versions.yml'),
                self.logger
            )

            # 执行流程|Execute pipeline
            if self.config.step:
                # 运行指定步骤|Run specified step
                success = self._run_single_step(self.config.step)
            else:
                # 运行完整流程|Run complete pipeline
                success = self._run_complete_pipeline()

            # 输出总结信息|Output summary
            self._print_summary(success)

            return 0 if success else 1

        except KeyboardInterrupt:
            self.logger.warning("流程被用户中断|Pipeline interrupted by user")
            return 130

        except Exception as e:
            self.logger.error(f"流程执行出错|Pipeline execution error: {str(e)}", exc_info=True)
            return 1

    def _print_header(self):
        """打印流程头部信息|Print pipeline header"""
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("  YaHS 高速Hi-C scaffolding流程|YaHS Hi-C Scaffolding Pipeline")
        self.logger.info("=" * 60)
        self.logger.info(f"开始时间|Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info(f"工作目录|Work directory: {self.config.output_dir}")
        self.logger.info(f"参考基因组|Reference genome: {self.config.ref_fa}")
        self.logger.info(f"Hi-C数据|Hi-C data: {self.config.hic_r1}, {self.config.hic_r2}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"酶切位点|Enzyme: {self.config.enzyme_seq}")
        self.logger.info("=" * 60)
        self.logger.info("")

    def _check_environment(self):
        """检查运行环境|Check runtime environment"""
        self.logger.info("检查运行环境|Checking runtime environment")

        # 检查磁盘空间|Check disk space
        import shutil
        total, used, free = shutil.disk_usage(self.config.output_dir)

        free_gb = free // (1024 ** 3)
        self.logger.info(f"可用磁盘空间|Available disk space: {free_gb}G")

        if free_gb < 10:
            self.logger.warning("磁盘空间不足10G，可能影响运行|Disk space < 10G, may affect execution")

        # 检查输入文件大小|Check input file sizes
        for file_desc, file_path in [
            ("参考基因组|Reference genome", self.config.ref_fa),
            ("Hi-C R1|Hi-C R1", self.config.hic_r1),
            ("Hi-C R2|Hi-C R2", self.config.hic_r2)
        ]:
            size_gb = Path(file_path).stat().st_size / (1024 ** 3)
            self.logger.info(f"{file_desc}: {size_gb:.2f}G")

        self.logger.info("环境检查完成|Environment check completed")
        self.logger.info("")

    def _run_complete_pipeline(self) -> bool:
        """
        运行完整流程|Run complete pipeline

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始完整流程|Starting complete pipeline")
        self.logger.info("")

        # 定义步骤|Define steps
        pipeline_steps = [
            ("1", "构建基因组索引|Build genome index", self.steps.step_1_indexing,
             self.config.step1_dir / f"{Path(self.config.ref_fa).stem}.bwt"),
            ("2", "Hi-C数据比对|Hi-C mapping", self.steps.step_2_mapping,
             self.config.final_bam),
            ("3", "YaHS染色体挂载|YaHS scaffolding", self.steps.step_3_scaffolding,
             self.config.final_scaffold_fa),
            ("4", "生成标准Hi-C热图|Generate standard Hi-C heatmap", self.steps.step_4_hic_standard,
             self.config.final_hic),
            ("5", "生成JBAT文件|Generate JBAT files", self.steps.step_5_jbat,
             self.config.jbat_assembly),
            ("6", "组装质量评估|Assembly quality assessment", self.steps.step_6_assessment,
             self.config.assembly_metrics),
        ]

        # 执行各步骤|Execute each step
        for step_num, step_desc, step_func, checkpoint_file in pipeline_steps:
            if self._is_step_completed(checkpoint_file):
                self.logger.info(f"跳过已完成步骤 {step_num}|Skipping completed step {step_num}: {step_desc}")
                self.logger.info("")
                continue

            if not self._run_step(step_num, step_desc, step_func):
                return False

        return True

    def _run_single_step(self, step_num: str) -> bool:
        """
        运行单个步骤|Run single step

        Args:
            step_num: 步骤编号|Step number (1-6)

        Returns:
            是否成功|Whether successful
        """
        step_map = {
            "1": ("构建基因组索引|Build genome index", self.steps.step_1_indexing),
            "2": ("Hi-C数据比对|Hi-C mapping", self.steps.step_2_mapping),
            "3": ("YaHS染色体挂载|YaHS scaffolding", self.steps.step_3_scaffolding),
            "4": ("生成标准Hi-C热图|Generate Hi-C heatmap", self.steps.step_4_hic_standard),
            "5": ("生成JBAT文件|Generate JBAT files", self.steps.step_5_jbat),
            "6": ("组装质量评估|Assembly assessment", self.steps.step_6_assessment),
        }

        if step_num not in step_map:
            self.logger.error(f"无效的步骤编号|Invalid step number: {step_num}")
            return False

        step_desc, step_func = step_map[step_num]
        return self._run_step(step_num, step_desc, step_func)

    def _run_step(self, step_num: str, step_desc: str, step_func) -> bool:
        """
        执行单个步骤|Execute single step

        Args:
            step_num: 步骤编号|Step number
            step_desc: 步骤描述|Step description
            step_func: 步骤函数|Step function

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("")
        self.logger.info(f"{'=' * 60}")
        self.logger.info(f"开始步骤 {step_num}|Starting step {step_num}: {step_desc}")
        self.logger.info(f"{'=' * 60}")
        self.logger.info("")

        step_start = time.time()

        try:
            success = step_func()

            step_runtime = calculate_runtime(step_start)

            if success:
                self.logger.info(f"步骤 {step_num} 完成|Step {step_num} completed [用时|Runtime: {step_runtime}]")
            else:
                self.logger.error(f"步骤 {step_num} 失败|Step {step_num} failed")

            return success

        except Exception as e:
            self.logger.error(f"步骤 {step_num} 执行异常|Step {step_num} execution error: {str(e)}", exc_info=True)
            return False

    def _is_step_completed(self, checkpoint_file) -> bool:
        """
        检查步骤是否已完成|Check if step is completed

        Args:
            checkpoint_file: 检查点文件|Checkpoint file

        Returns:
            是否已完成|Whether completed
        """
        if self.config.force_rerun:
            return False

        return Path(checkpoint_file).exists()

    def _print_summary(self, success: bool):
        """
        打印流程总结|Print pipeline summary

        Args:
            success: 是否成功|Whether successful
        """
        total_runtime = calculate_runtime(self.start_time)

        self.logger.info("")
        self.logger.info("=" * 60)
        if success:
            self.logger.info("流程执行成功|Pipeline executed successfully")
        else:
            self.logger.error("流程执行失败|Pipeline execution failed")
        self.logger.info("=" * 60)
        self.logger.info(f"总运行时间|Total runtime: {total_runtime}")
        self.logger.info(f"工作目录|Work directory: {self.config.output_dir}")
        self.logger.info(f"日志文件|Log file: {self.config.pipeline_log}")
        self.logger.info("")

        if success:
            self.logger.info("主要输出文件|Main output files:")
            self.logger.info(f"  • Scaffold序列|Scaffold sequences: {self.config.final_scaffold_fa}")
            self.logger.info(f"  • Scaffold AGP|Scaffold AGP: {self.config.final_scaffold_agp}")

            if self.config.final_hic.exists():
                self.logger.info(f"  • 标准Hi-C热图|Standard Hi-C heatmap: {self.config.final_hic}")

            if self.config.jbat_hic.exists():
                self.logger.info(f"  • JBAT Hi-C热图|JBAT Hi-C heatmap: {self.config.jbat_hic}")
                self.logger.info(f"  • JBAT assembly文件|JBAT assembly: {self.config.jbat_assembly}")

            if self.config.assembly_metrics.exists():
                self.logger.info(f"  • 质量评估报告|Quality assessment: {self.config.assembly_metrics}")

        self.logger.info("")
        self.logger.info("=" * 60)


def main():
    """
    主函数（用于模块测试）|Main function (for module testing)

    注意：实际入口点在cli/commands/yahs.py
    Note: Actual entry point is at cli/commands/yahs.py
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="YaHS Hi-C Scaffolding Pipeline"
    )

    parser.add_argument('-r', '--ref', required=True,
                        help='参考基因组FASTA|Reference genome FASTA')
    parser.add_argument('-1', '--hic-r1', required=True,
                        help='Hi-C R1文件|Hi-C R1 file')
    parser.add_argument('-2', '--hic-r2', required=True,
                        help='Hi-C R2文件|Hi-C R2 file')
    parser.add_argument('-o', '--output-dir', default='./yahs_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads')
    parser.add_argument('-e', '--enzyme', default='GATC',
                        help='酶切位点|Enzyme sequence')
    parser.add_argument('-s', '--step',
                        choices=['1', '2', '3', '4', '5', '6'],
                        help='运行指定步骤|Run specified step')
    parser.add_argument('--force-rerun', action='store_true',
                        help='强制重新运行|Force rerun all steps')

    args = parser.parse_args()

    # 创建配置|Create configuration
    config = YaHSConfig(
        ref_fa=args.ref,
        hic_r1=args.hic_r1,
        hic_r2=args.hic_r2,
        output_dir=args.output_dir,
        threads=args.threads,
        enzyme_seq=args.enzyme,
        step=args.step,
        force_rerun=args.force_rerun
    )

    # 运行流程|Run pipeline
    pipeline = YaHSPipeline(config)
    sys.exit(pipeline.run())


if __name__ == '__main__':
    main()
