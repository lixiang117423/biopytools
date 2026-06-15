"""
PICRUSt2流程管理模块|PICRUSt2 Pipeline Management Module
"""

import glob
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from .config import Picrust2Config
from .utils import (
    CommandRunner, build_conda_command,
    reorganize_outputs, generate_software_versions_yml, format_number,
    prepare_input_table, preprocess_biom_table, detect_input_format,
    clean_fasta_headers, format_pathway_table
)


class Picrust2Pipeline:
    """PICRUSt2分析流程|PICRUSt2 Analysis Pipeline"""

    def __init__(self, config: Picrust2Config, logger=None):
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger, working_dir=config.output_dir)

    def run_pipeline(self):
        """运行完整PICRUSt2分析流程|Run complete PICRUSt2 analysis pipeline"""
        start_time = datetime.now()
        self.logger.info(f"开始PICRUSt2分析流程|Starting PICRUSt2 analysis pipeline")
        self.logger.info(f"代表序列文件|Study FASTA: {self.config.study_fasta}")
        self.logger.info(f"特征表文件|Input table: {self.config.input_table}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        # 确定使用哪个pipeline脚本
        pipeline_script = self.config._resolve_pipeline_script()
        pipeline_type = self._detect_pipeline_type(pipeline_script)
        self.logger.info(f"使用流程类型|Pipeline type: {pipeline_type}")

        # 断点续传检查|Checkpoint resume check
        if self._is_completed():
            self.logger.info("检测到已完成的分析，跳过|Detected completed analysis, skipping")
            self._generate_versions(start_time)
            return True

        # 自动识别输入特征表格式，Excel则转换为TSV
        # Auto-detect input table format, convert Excel to TSV if needed
        actual_input_table = prepare_input_table(
            self.config.input_table, self.config.work_dir, self.logger
        )

        # BIOM表预处理：自动检测并修复轴方向(行列转置)
        # BIOM preprocess: auto-detect and fix axis orientation (transpose)
        preprocess_dir = self.config.output_dir
        if detect_input_format(actual_input_table) == 'biom':
            actual_input_table = preprocess_biom_table(
                actual_input_table, self.config.study_fasta,
                preprocess_dir, self.logger
            )

        # FASTA header清理：去除描述信息，确保ID与BIOM表一致
        # FASTA header cleanup: strip descriptions to ensure IDs match BIOM table
        actual_fasta = clean_fasta_headers(
            self.config.study_fasta, preprocess_dir, self.logger
        )

        # 构建命令|Build command
        cmd = self._build_picrust2_command(pipeline_script, pipeline_type, actual_input_table, actual_fasta)

        # 执行PICRUSt2 pipeline
        success, stdout, stderr = self.cmd_runner.run_command(
            cmd, description=f"PICRUSt2 pipeline ({pipeline_type})"
        )

        if not success:
            self.logger.error("PICRUSt2流程执行失败|PICRUSt2 pipeline execution failed")
            sys.exit(1)

        # 检查输出文件|Check output files
        if not self._validate_outputs(self.config.work_dir):
            self.logger.error("PICRUSt2输出文件不完整|PICRUSt2 output files incomplete")
            sys.exit(1)

        # 重组输出文件|Reorganize output files
        self.logger.info("重组输出文件到标准目录|Reorganizing outputs to standard directories")
        reorganize_outputs(
            work_dir=self.config.work_dir,
            output_config={
                'placement_dir': self.config.placement_dir,
                'hsp_dir': self.config.hsp_dir,
                'metagenome_dir': self.config.metagenome_dir,
                'pathway_dir': self.config.pathway_dir,
            },
            logger=self.logger
        )

        # 格式化通路丰度表为易读Excel|Format pathway table to readable Excel
        format_pathway_table(self.config.pathway_dir, self.logger)

        # 生成版本信息|Generate version info
        self._generate_versions(start_time)

        elapsed = int((datetime.now() - start_time).total_seconds())
        hours, remainder = divmod(elapsed, 3600)
        minutes, seconds = divmod(remainder, 60)
        self.logger.info(
            f"PICRUSt2分析流程完成|PICRUSt2 analysis pipeline completed "
            f"(耗时|Runtime: {hours}h{minutes}m{seconds}s)"
        )
        return True

    def _detect_pipeline_type(self, pipeline_script: Optional[str]) -> str:
        """检测流程类型|Detect pipeline type"""
        if pipeline_script is None:
            return "unknown"
        basename = os.path.basename(pipeline_script)
        if 'singleRef' in basename:
            return "single_ref"
        return "split_domain"

    def _build_picrust2_command(self, pipeline_script: str, pipeline_type: str,
                                 input_table: str = None, fasta_path: str = None) -> list:
        """构建PICRUSt2 pipeline命令|Build PICRUSt2 pipeline command"""
        args = []

        # 必需参数|Required parameters
        args.extend(['-s', fasta_path or self.config.study_fasta])
        args.extend(['-i', input_table or self.config.input_table])
        args.extend(['-o', self.config.work_dir])

        # 流程参数|Pipeline parameters
        args.extend(['-p', str(self.config.threads)])
        args.extend(['--max_nsti', str(self.config.max_nsti)])
        args.extend(['-t', self.config.placement_tool])
        args.extend(['-m', self.config.hsp_method])
        args.extend(['-e', str(self.config.edge_exponent)])
        args.extend(['--min_align', str(self.config.min_align)])

        if self.config.min_reads > 1:
            args.extend(['--min_reads', str(self.config.min_reads)])

        if self.config.min_samples > 1:
            args.extend(['--min_samples', str(self.config.min_samples)])

        # 功能数据库|Functional databases
        args.extend(['--in_traits', self.config.in_traits])

        # 输出选项|Output options
        if self.config.stratified:
            args.append('--stratified')

        if self.config.coverage:
            args.append('--coverage')

        if self.config.skip_norm:
            args.append('--skip_norm')

        if self.config.per_sequence_contrib:
            args.append('--per_sequence_contrib')

        # 通路控制|Pathway control
        if self.config.skip_pathways:
            args.append('--no_pathways')
        else:
            if self.config.skip_minpath:
                args.append('--skip_minpath')
            if self.config.no_gap_fill:
                args.append('--no_gap_fill')

        # 其他选项|Other options
        if self.config.remove_intermediate:
            args.append('--remove_intermediate')

        if self.config.verbose:
            args.append('--verbose')

        # 单参考流程的额外参数|Single ref pipeline extra params
        if pipeline_type == "single_ref":
            if self.config.in_ref_dir:
                args.extend(['-r', self.config.in_ref_dir])

        # 包装conda run
        cmd = build_conda_command(pipeline_script, args)
        return cmd

    def _is_completed(self) -> bool:
        """检查分析是否已完成（断点续传）|Check if analysis is completed (checkpoint resume)"""
        # 检查关键输出文件|Check key output files
        if self.config.skip_pathways:
            # 跳过通路时只检查metagenome输出|When skipping pathways, only check metagenome output
            marker = os.path.join(self.config.metagenome_dir, "pred_metagenome_unstrat.tsv.gz")
        else:
            marker = os.path.join(self.config.pathway_dir, "path_abun_unstrat.tsv.gz")

        if os.path.exists(marker):
            self.logger.info(f"检测到已完成输出|Detected completed output: {marker}")
            return True

        return False

    def _validate_outputs(self, work_dir: str) -> bool:
        """验证PICRUSt2输出文件完整性|Validate PICRUSt2 output completeness"""
        self.logger.info("验证输出文件完整性|Validating output completeness")

        # 必须存在的输出文件|Required output files
        required_files = ['pred_metagenome_unstrat.tsv.gz']
        if not self.config.skip_pathways:
            required_files.append('path_abun_unstrat.tsv.gz')

        # 在work_dir及其子目录中查找|Search in work_dir and subdirectories
        all_found = True
        for fname in required_files:
            # 直接在work_dir下查找
            direct = os.path.join(work_dir, fname)
            # 在子目录中查找
            matches = glob.glob(os.path.join(work_dir, '**', fname), recursive=True)

            if os.path.exists(direct):
                self.logger.info(f"输出文件验证通过|Output file verified: {fname}")
            elif matches:
                self.logger.info(f"输出文件验证通过|Output file verified: {fname}")
            else:
                self.logger.warning(f"输出文件缺失|Output file missing: {fname}")
                all_found = False

        return all_found

    def _generate_versions(self, start_time: datetime):
        """生成版本信息文件|Generate version info file"""
        try:
            version_file = generate_software_versions_yml(
                self.config.output_dir, self.config, start_time
            )
            self.logger.info(f"版本信息已生成|Version info generated: {version_file}")
        except Exception as e:
            self.logger.warning(f"生成版本信息失败|Failed to generate version info: {e}")
