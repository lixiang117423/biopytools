"""
系统发育树构建模块|Phylogenetic Tree Building Module
"""

import subprocess
from pathlib import Path

from ..common.conda_runner import build_conda_command


class FastTreeBuilder:
    """FastTree系统发育树构建器|FastTree Phylogenetic Tree Builder"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_tree(self, alignment_file: Path, output_file: Path, seq_type: str) -> bool:
        """构建系统发育树|Build phylogenetic tree"""
        self.logger.info("=" * 60)
        self.logger.info("开始构建系统发育树|Starting phylogenetic tree construction")
        self.logger.info("=" * 60)
        
        # 构建FastTree命令|Build FastTree command
        fasttree_args = []
        
        # 根据序列类型添加参数
        if seq_type == 'nucleotide':
            fasttree_args.append("-nt")
            self.logger.info(" 使用核酸模型|Using nucleotide model")
        else:
            self.logger.info(" 使用蛋白质模型|Using protein model")
        
        # 添加用户自定义参数
        if self.config.fasttree_params:
            fasttree_args.extend(self.config.fasttree_params.split())
        
        # 添加输入文件|Add input file
        fasttree_args.append(str(alignment_file.resolve()))

        cmd = build_conda_command(self.config.fasttree_path, fasttree_args)
        self.logger.info(f"命令|Command: {' '.join(cmd)} > {output_file.resolve()}")

        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            success = (result.returncode == 0)
            if not success and result.stderr:
                self.logger.error(f" 系统发育树构建失败|Phylogenetic tree construction failed: {result.stderr}")
        except Exception as e:
            self.logger.error(f" 系统发育树构建失败|Phylogenetic tree construction failed: {e}")
            return False
        
        if success:
            self.logger.info(f" 系统发育树构建完成|Phylogenetic tree completed: {output_file}")
        else:
            self.logger.error(f" 系统发育树构建失败|Phylogenetic tree construction failed")
        
        return success
