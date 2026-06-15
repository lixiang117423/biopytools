"""
HMMsearch运行器模块|HMMsearch Runner Module
"""

import subprocess
import os
import shlex
import shutil
import re
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name"""
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name
    return None


def build_conda_command_string(command: str, args: str) -> str:
    """构建conda run命令字符串（用于需要shell特性的命令）|Build conda run command string (for commands needing shell features)"""
    conda_env = get_conda_env(command)
    if conda_env:
        return f"conda run -n {conda_env} {command} {args}"
    else:
        return f"{command} {args}"


class HMMsearchRunner:
    """HMMsearch命令运行器|HMMsearch Command Runner"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def run(self) -> str:
        """
        运行hmmsearch命令|Run hmmsearch command

        Returns:
            str: 生成的domtblout文件路径|Path to generated domtblout file
        """
        self.logger.info("=" * 60)
        self.logger.info("开始运行HMMsearch|Starting HMMsearch")
        self.logger.info("=" * 60)

        protein_fastas = self.config.protein_fastas
        total = len(protein_fastas)

        if total > 1:
            self.logger.info(f"检测到{total}个蛋白序列文件|Detected {total} protein FASTA files")

        final_domtblout = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.domtblout"
        )

        # 如果最终domtblout已存在，跳过所有|If final domtblout exists, skip all
        if os.path.exists(final_domtblout):
            self.logger.info(f"domtblout文件已存在，跳过hmmsearch|domtblout file already exists, skipping hmmsearch: {final_domtblout}")
            return final_domtblout

        # 逐文件运行hmmsearch|Run hmmsearch for each file
        individual_domtblouts = []
        for i, protein_fasta in enumerate(protein_fastas, 1):
            self.logger.info(f"处理文件|Processing file [{i}/{total}]: {os.path.basename(protein_fasta)}")

            # 单个文件的domtblout输出|Per-file domtblout output
            suffix = Path(protein_fasta).stem
            per_file_domtblout = os.path.join(
                self.config.output_dir,
                f"{self.config.output_prefix}_{suffix}.domtblout"
            )

            # 断点续传：单个文件结果已存在则跳过|Checkpoint resume: skip if per-file result exists
            if os.path.exists(per_file_domtblout):
                self.logger.info(f"跳过已完成文件|Skipping completed file: {os.path.basename(protein_fasta)}")
                individual_domtblouts.append(per_file_domtblout)
                continue

            cmd_args = self._build_command_args(per_file_domtblout, protein_fasta)
            hmmsearch_cmd = os.path.basename(self.config.hmmsearch_path)
            cmd = build_conda_command_string(hmmsearch_cmd, cmd_args)

            self.logger.info(f"命令|Command: {cmd}")

            try:
                subprocess.run(
                    cmd,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True
                )
                self.logger.info(f"完成|Completed [{i}/{total}]: {os.path.basename(protein_fasta)}")
            except subprocess.CalledProcessError as e:
                self.logger.error(f"HMMsearch执行失败|HMMsearch execution failed: {os.path.basename(protein_fasta)}")
                self.logger.error(f"错误代码|Error code: {e.returncode}")
                self.logger.error(f"错误信息|Error message: {e.stderr}")
                raise

            individual_domtblouts.append(per_file_domtblout)

        # 合并所有domtblout文件|Merge all domtblout files
        if len(individual_domtblouts) > 1:
            self.logger.info(f"合并{len(individual_domtblouts)}个domtblout文件|Merging {len(individual_domtblouts)} domtblout files")
            self._merge_domtblouts(individual_domtblouts, final_domtblout)
        else:
            os.rename(individual_domtblouts[0], final_domtblout)

        self.logger.info("HMMsearch全部执行成功|All HMMsearch executions completed successfully")
        return final_domtblout

    def _merge_domtblouts(self, input_files: list, output_file: str):
        """合并多个domtblout文件（保留第一个文件的注释头）|Merge multiple domtblout files (keep header comments from first file only)"""
        with open(output_file, 'w') as f_out:
            first = True
            for f_in_path in input_files:
                with open(f_in_path, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('#'):
                            if not first:
                                continue
                        f_out.write(line)
                first = False

    def _build_command_args(self, domtblout_file: str, protein_fasta: str) -> str:
        """
        构建hmmsearch命令参数|Build hmmsearch command arguments

        Args:
            domtblout_file: domtblout输出文件路径|domtblout output file path
            protein_fasta: 蛋白序列文件路径|Protein FASTA file path

        Returns:
            str: hmmsearch命令参数字符串|hmmsearch command arguments string
        """
        cmd_parts = []

        # 基础参数|Basic parameters
        cmd_parts.append(f"--cpu {self.config.threads}")

        # 输出参数|Output parameters
        cmd_parts.append(f"--domtblout {shlex.quote(domtblout_file)}")

        # 阈值参数|Threshold parameters
        if self.config.evalue_cutoff:
            cmd_parts.append(f"-E {self.config.evalue_cutoff}")
        if self.config.score_cutoff:
            cmd_parts.append(f"-T {self.config.score_cutoff}")

        # Cut选项|Cut options
        if self.config.use_cut_tc:
            cmd_parts.append("--cut_tc")
        elif self.config.use_cut_ga:
            cmd_parts.append("--cut_ga")
        elif self.config.use_cut_nc:
            cmd_parts.append("--cut_nc")

        # 输入文件|Input files
        cmd_parts.append(shlex.quote(self.config.hmm_file))
        cmd_parts.append(shlex.quote(protein_fasta))

        return " ".join(cmd_parts)
