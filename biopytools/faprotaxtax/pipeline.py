"""
FAPROTAX流程编排模块|FAPROTAX Pipeline Orchestration Module
"""

import os
import shutil
import sys
import time
from datetime import datetime
from typing import List

from .config import FaprotaxtaxConfig
from .utils import CommandRunner, detect_input_format, generate_software_versions_yml


class FaprotaxtaxPipeline:
    """FAPROTAX流程编排器|FAPROTAX Pipeline Orchestrator"""

    def __init__(self, config: FaprotaxtaxConfig, logger=None):
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger, config.output_dir) if logger else None

    def _is_completed(self) -> bool:
        """检查流程是否已完成|Check if pipeline is already completed"""
        # 检查关键输出文件|Check key output file
        is_biom = detect_input_format(self.config.input_table) == 'biom'
        ext = '.biom' if (is_biom and self.config.output_format != 'classical') else '.tsv'
        output_file = os.path.join(self.config.collapsed_dir, f"functional_table{ext}")
        return os.path.exists(output_file)

    def _build_collapse_command(self) -> List[str]:
        """
        构建collapse_table.py命令|Build collapse_table.py command

        Returns:
            完整命令列表|Complete command list
        """
        cmd = [self.config.python_interpreter, self.config.collapse_table_path]

        # 必需参数|Required parameters
        cmd.extend(['-i', self.config.input_table])
        cmd.extend(['-g', self.config.groups_file])

        # 输出文件（写入work目录）|Output files (write to work dir)
        is_biom = detect_input_format(self.config.input_table) == 'biom'
        default_ext = '.biom' if (is_biom and self.config.output_format != 'classical') else '.tsv'
        out_ext = '.biom' if self.config.output_format == 'BIOM' else default_ext

        collapsed_out = os.path.join(self.config.work_dir, f"functional_table{out_ext}")
        report_out = os.path.join(self.config.work_dir, "faprotaxtax_report.txt")

        cmd.extend(['-o', collapsed_out])
        cmd.extend(['-r', report_out])

        # 输出格式|Output format
        if self.config.output_format != 'auto':
            cmd.extend(['--output_format_collapsed', self.config.output_format])

        # 核心参数|Core parameters
        if self.config.collapse_by_metadata:
            cmd.extend(['--collapse_by_metadata', self.config.collapse_by_metadata])

        if self.config.group_leftovers_as:
            cmd.extend(['--group_leftovers_as', self.config.group_leftovers_as])

        if self.config.normalize != 'none':
            cmd.extend(['-n', self.config.normalize])

        if self.config.average != 'none':
            cmd.extend(['--average', self.config.average])

        if self.config.row_names_are_in_column:
            cmd.extend(['-d', self.config.row_names_are_in_column])

        # 控制参数|Control parameters
        if self.config.verbose:
            cmd.append('-v')

        if self.config.force:
            cmd.append('-f')

        # 默认使用空格分隔词匹配模式（FAPROTAX标准）|Default to word matching mode
        cmd.extend(['--group_members_defined_as', 'words'])

        return cmd

    def _validate_outputs(self) -> bool:
        """
        验证collapse_table.py输出|Validate collapse_table.py outputs

        Returns:
            bool: 输出是否有效|Whether outputs are valid
        """
        work_dir = self.config.work_dir
        missing = []

        # 检查功能表|Check functional table
        for ext in ['.tsv', '.biom']:
            if os.path.exists(os.path.join(work_dir, f"functional_table{ext}")):
                break
        else:
            missing.append("functional_table (.tsv or .biom)")

        # 检查报告（如果指定了--group_leftovers_as）|Check report (if --group_leftovers_as specified)
        report_file = os.path.join(work_dir, "faprotaxtax_report.txt")
        if not os.path.exists(report_file) and self.config.group_leftovers_as:
            missing.append("faprotaxtax_report.txt")

        if missing:
            self.logger.warning(f"部分输出文件缺失|Some output files missing: {', '.join(missing)}")
            return False

        return True

    def _reorganize_outputs(self) -> None:
        """
        重组输出文件到编号目录|Reorganize output files into numbered directories
        """
        self.logger.info("重组输出文件|Reorganizing output files")

        work_dir = self.config.work_dir

        # 移动功能表到01_collapsed/|Move functional table to 01_collapsed/
        for ext in ['.tsv', '.biom']:
            src_func = os.path.join(work_dir, f"functional_table{ext}")
            if os.path.exists(src_func):
                dst_func = os.path.join(self.config.collapsed_dir, f"functional_table{ext}")
                shutil.copy2(src_func, dst_func)
                self.logger.info(f"功能表已保存|Functional table saved: {dst_func}")
                break

        # 移动报告到02_report/|Move report to 02_report/
        src_report = os.path.join(work_dir, "faprotaxtax_report.txt")
        if os.path.exists(src_report):
            dst_report = os.path.join(self.config.report_dir, "faprotaxtax_report.txt")
            shutil.copy2(src_report, dst_report)
            self.logger.info(f"报告已保存|Report saved: {dst_report}")

        # 清理work目录|Clean up work directory
        if os.path.exists(work_dir):
            try:
                shutil.rmtree(work_dir)
                self.logger.info(f"临时目录已清理|Temporary directory cleaned: {work_dir}")
            except Exception as e:
                self.logger.warning(f"清理临时目录失败|Failed to clean work directory: {e}")

    def _generate_versions(self, start_time: datetime) -> None:
        """生成软件版本信息|Generate software version information"""
        self.logger.info("生成软件版本信息|Generating software version information")
        generate_software_versions_yml(self.config.output_dir, self.config, start_time)

    def run_pipeline(self) -> bool:
        """
        运行FAPROTAX流程|Run FAPROTAX pipeline

        Returns:
            bool: 流程是否成功|Whether pipeline succeeded
        """
        start_time = datetime.now()

        # 记录流程启动信息|Log pipeline startup info
        self.logger.info("=" * 60)
        self.logger.info("FAPROTAX流程开始|FAPROTAX pipeline started")
        self.logger.info(f"输入文件|Input file: {self.config.input_table}")
        self.logger.info(f"功能组数据库|Groups database: {self.config.groups_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"collapse_table.py路径|collapse_table.py path: {self.config.collapse_table_path}")

        # 记录关键参数|Log key parameters
        params = []
        if self.config.collapse_by_metadata:
            params.append(f"collapse_by_metadata={self.config.collapse_by_metadata}")
        if self.config.group_leftovers_as:
            params.append(f"group_leftovers_as={self.config.group_leftovers_as}")
        if self.config.normalize != 'none':
            params.append(f"normalize={self.config.normalize}")
        if self.config.average != 'none':
            params.append(f"average={self.config.average}")
        if params:
            self.logger.info(f"参数|Parameters: {', '.join(params)}")

        # 断点续传检查|Checkpoint resume check
        if self._is_completed():
            self.logger.info("检测到已完成输出，跳过流程|Detected completed outputs, skipping pipeline")
            self._generate_versions(start_time)
            self.logger.info("FAPROTAX流程完成（已跳过）|FAPROTAX pipeline completed (skipped)")
            return True

        # 检测输入格式|Detect input format
        input_format = detect_input_format(self.config.input_table)
        self.logger.info(f"输入格式|Input format: {input_format}")

        # 构建命令|Build command
        cmd = self._build_collapse_command()

        # 执行collapse_table.py|Execute collapse_table.py
        self.logger.info("执行FAPROTAX功能注释|Executing FAPROTAX functional annotation")
        success, stdout, stderr = self.cmd_runner.run_command(
            cmd,
            f"FAPROTAX collapse_table.py - {input_format}输入|{input_format} input"
        )

        if not success:
            self.logger.error("collapse_table.py执行失败|collapse_table.py execution failed")
            self.logger.error(f"终止流程|Pipeline terminated with error")
            return False

        # 验证输出|Validate outputs
        if not self._validate_outputs():
            self.logger.error("输出验证失败，流程终止|Output validation failed, pipeline terminated")
            return False

        # 重组输出|Reorganize outputs
        self._reorganize_outputs()

        # 生成版本信息|Generate version information
        self._generate_versions(start_time)

        # 记录完成信息|Log completion info
        elapsed = (datetime.now() - start_time).total_seconds()
        self.logger.info("=" * 60)
        self.logger.info(f"FAPROTAX流程完成|FAPROTAX pipeline completed")
        self.logger.info(f"总耗时|Total time: {elapsed:.1f}s")
        self.logger.info(f"功能表|Functional table: {self.config.collapsed_dir}/")
        if self.config.group_leftovers_as:
            self.logger.info(f"报告|Report: {self.config.report_dir}/")
        self.logger.info(f"日志|Log: {self.config.log_dir}/")

        return True
