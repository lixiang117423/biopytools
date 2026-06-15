"""
OrthoFinder执行器模块|OrthoFinder Runner Module
"""

import os
import shutil
from pathlib import Path
from .utils import CommandRunner, build_conda_command


class OrthoFinderRunner:
    """OrthoFinder执行器|OrthoFinder Runner"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.orthofinder_output_dir = None
        self.actual_project_name = None

    def check_existing_results(self) -> bool:
        """检查已有OrthoFinder结果|Check existing OrthoFinder results"""
        locations_to_check = [
            (Path(self.config.input_dir) / "OrthoFinder", "输入目录|Input directory"),
            (Path(self.config.output_dir) / "01_orthofinder_analysis" / "orthofinder_results", "输出目录|Output directory"),
            (Path(self.config.output_dir) / "Results_OrthoFinder", "旧输出目录|Old output directory"),
            (Path(self.config.output_dir) / "OrthoFinder_Results", "旧输出目录备用名称|Old output directory (alt)")
        ]

        results_pattern = f"Results_{self.config.project_name}"

        for check_path, location_name in locations_to_check:
            self.logger.info(f"检查{location_name}|Checking {location_name}: {check_path}")

            if check_path.exists():
                if check_path.name.startswith("Results_") or check_path.name.startswith("OrthoFinder_Results"):
                    results_dir = check_path
                    self.logger.info(f"找到直接结果目录|Found direct results directory: {results_dir}")
                else:
                    possible_dirs = list(check_path.glob(f"{results_pattern}*"))
                    if not possible_dirs:
                        self.logger.info(f"未找到匹配的结果目录|No matching results directory found")
                        continue

                    results_dir = max(possible_dirs, key=os.path.getctime)
                    self.logger.info(f"找到结果目录|Found results directory: {results_dir}")

                orthogroups_file = results_dir / "Orthogroups" / "Orthogroups.tsv"
                gene_count_file = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"

                if not orthogroups_file.exists():
                    orthogroups_file = results_dir / "Orthogroups" / "Orthogroups.txt"

                self.logger.info(f"检查关键文件存在性:")
                self.logger.info(f"  Orthogroups.tsv: {orthogroups_file.exists()}")
                self.logger.info(f"  GeneCount.tsv: {gene_count_file.exists()}")

                if orthogroups_file.exists() and gene_count_file.exists():
                    self.logger.info(f"发现有效的OrthoFinder结果|Found valid OrthoFinder results: {results_dir}")
                    self.orthofinder_output_dir = results_dir
                    return True

        self.logger.info("未发现有效结果，需要运行OrthoFinder|No valid results found, need to run OrthoFinder")
        return False

    def run_orthofinder(self) -> bool:
        """运行OrthoFinder分析|Run OrthoFinder analysis"""

        if self.config.skip_orthofinder:
            self.logger.info("跳过OrthoFinder步骤|Skipping OrthoFinder step")
            return self.check_existing_results()

        if self.config.resume_from_existing and not self.config.force_overwrite:
            if self.check_existing_results():
                if self.config.generate_trees and self.orthofinder_output_dir:
                    species_tree_dir = self.orthofinder_output_dir / "Species_Tree"
                    if not species_tree_dir.exists():
                        self.logger.info("已有结果不包含物种树，需要重新运行|Existing results lack species tree, need to re-run")
                    else:
                        self.logger.info("使用已有结果续跑|Resuming from existing results")
                        return True
                else:
                    self.logger.info("使用已有结果续跑|Resuming from existing results")
                    return True

        if self.config.force_overwrite:
            self._cleanup_existing_results()

        self.logger.info(f"开始OrthoFinder分析|Starting OrthoFinder analysis")
        self.logger.info(f"参数设置|Parameter settings:")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  搜索程序|Search program: {self.config.search_program}")
        self.logger.info(f"  MCL inflation: {self.config.mcl_inflation}")

        # 构建OrthoFinder命令参数列表|Build OrthoFinder command argument list
        cmd_args = [
            '-f', self.config.input_dir,
            '-t', str(self.config.threads),
            '-S', self.config.search_program,
            '-I', str(self.config.mcl_inflation),
            '-n', self.config.project_name
        ]

        self.actual_project_name = self.config.project_name

        if self.config.sequence_type == 'dna':
            cmd_args.append("-d")

        if self.config.basic_analysis_only:
            cmd_args.append("-og")
        elif self.config.generate_trees:
            cmd_args.extend(['-A', self.config.msa_program, '-T', self.config.tree_program])

        cmd = build_conda_command(self.config.orthofinder_path, cmd_args)
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        success = self.cmd_runner.run(cmd, "OrthoFinder同源基因群分析|OrthoFinder orthogroup analysis", long_running=True)

        if not success:
            self._report_log_file_error()

        if success:
            self._locate_output_directory()
            if self.orthofinder_output_dir is None:
                self.logger.error("OrthoFinder命令执行成功但未找到结果目录，可能是已有同名项目的不完整结果|OrthoFinder command succeeded but no results directory found, possibly incomplete results from same project name")
                self._report_log_file_error()
                return False
            self.logger.info(f"OrthoFinder分析完成|OrthoFinder analysis completed")
            self.logger.info(f"结果目录|Results directory: {self.orthofinder_output_dir}")

        return success

    def _report_log_file_error(self):
        """失败时读取OrthoFinder的Log.txt报告详细错误|Read OrthoFinder Log.txt on failure"""
        orthofinder_dir = Path(self.config.input_dir) / "OrthoFinder"
        if not orthofinder_dir.exists():
            return

        results_dirs = sorted(orthofinder_dir.glob("Results_*"), key=os.path.getctime, reverse=True)
        if not results_dirs:
            return

        log_file = results_dirs[0] / "Log.txt"
        if log_file.exists():
            try:
                log_content = log_file.read_text(encoding='utf-8')
                self.logger.error(f"OrthoFinder日志|OrthoFinder Log.txt:\n{log_content}")
            except Exception as e:
                self.logger.error(f"无法读取OrthoFinder日志|Cannot read OrthoFinder log: {e}")

    def _locate_output_directory(self):
        """定位OrthoFinder输出目录并复制到指定位置|Locate OrthoFinder output directory and copy to specified location"""
        results_pattern = f"Results_{self.actual_project_name}"

        input_path = Path(self.config.input_dir)
        orthofinder_path = input_path / "OrthoFinder"
        possible_dirs = list(orthofinder_path.glob(f"{results_pattern}*"))

        if possible_dirs:
            source_dir = max(possible_dirs, key=os.path.getctime)
            self.logger.info(f"找到OrthoFinder结果目录|Found OrthoFinder results directory: {source_dir}")

            analysis_dir = Path(self.config.output_dir) / "01_orthofinder_analysis"
            analysis_dir.mkdir(parents=True, exist_ok=True)

            target_dir = analysis_dir / "orthofinder_results"
            if target_dir.exists():
                shutil.rmtree(target_dir)

            shutil.copytree(str(source_dir), str(target_dir))
            self.orthofinder_output_dir = target_dir
            self.logger.info(f"结果已复制到|Results copied to: {target_dir}")

            if orthofinder_path.exists():
                target_work_dir = analysis_dir / "orthofinder_working_dir"
                if target_work_dir.exists():
                    shutil.rmtree(target_work_dir)
                shutil.copytree(str(orthofinder_path), str(target_work_dir))
                self.logger.info(f"工作目录已复制到|Working directory copied to: {target_work_dir}")
        else:
            self.logger.error("未找到OrthoFinder结果目录|OrthoFinder results directory not found")
            self.orthofinder_output_dir = None

    def get_orthogroups_file(self) -> Path:
        """获取同源基因群文件路径|Get orthogroups file path"""
        if not self.orthofinder_output_dir:
            raise RuntimeError("OrthoFinder结果目录未找到|OrthoFinder results directory not found")

        orthogroups_dir = Path(self.orthofinder_output_dir) / "Orthogroups"

        for filename in ["Orthogroups.tsv", "Orthogroups.txt"]:
            orthogroups_file = orthogroups_dir / filename
            if orthogroups_file.exists():
                return orthogroups_file

        raise FileNotFoundError(f"未找到同源基因群详细文件|Orthogroups detailed file not found in: {orthogroups_dir}")

    def get_gene_count_file(self) -> Path:
        """获取基因计数文件路径|Get gene count file path"""
        if not self.orthofinder_output_dir:
            raise RuntimeError("OrthoFinder结果目录未找到|OrthoFinder results directory not found")

        orthogroups_dir = Path(self.orthofinder_output_dir) / "Orthogroups"
        gene_count_file = orthogroups_dir / "Orthogroups.GeneCount.tsv"

        if gene_count_file.exists():
            return gene_count_file

        raise FileNotFoundError(f"未找到基因计数文件|Gene count file not found: {gene_count_file}")

    def _cleanup_existing_results(self):
        """清理已有结果|Cleanup existing results"""
        input_path = Path(self.config.input_dir)
        orthofinder_path = input_path / "OrthoFinder"

        if orthofinder_path.exists():
            shutil.rmtree(orthofinder_path)
            self.logger.info("已清理旧的OrthoFinder结果|Cleaned up old OrthoFinder results")
