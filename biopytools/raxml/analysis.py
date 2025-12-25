"""
🌳 RAxML系统发育分析核心模块 | RAxML Phylogenetic Analysis Core Module
"""

import os
import random
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from .utils import CommandRunner

class RAxMLAnalyzer:
    """RAxML分析器 | RAxML Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def check_sequence_file(self) -> bool:
        """检查序列文件格式 | Check sequence file format"""
        self.logger.info("🔍 检查序列文件格式 | Checking sequence file format")
        
        # 使用RAxML的-f c选项来检查文件 | Use RAxML -f c option to check file
        cmd = (
            f"{self.config.raxml_path} -f c "
            f"-s {self.config.sequence_file} "
            f"-m {self.config.model} "
            f"-n check_alignment"
        )
        
        if self.config.silent:
            cmd += " --silent"
        
        success = self.cmd_runner.run(cmd, "序列文件格式验证 | Sequence file format validation")
        
        # 清理临时文件 | Clean up temporary files
        temp_files = [
            "RAxML_info.check_alignment",
            "RAxML_parsimonyTree.check_alignment",
            "RAxML_result.check_alignment"
        ]
        
        for temp_file in temp_files:
            temp_path = self.config.output_path / temp_file
            if temp_path.exists():
                temp_path.unlink()
        
        return success
    
    def generate_random_seeds(self) -> Tuple[int, int]:
        """生成随机种子 | Generate random seeds"""
        if self.config.parsimony_seed is None:
            parsimony_seed = random.randint(1, 1000000)
            self.logger.info(f"🎲 生成parsimony随机种子 | Generated parsimony random seed: {parsimony_seed}")
        else:
            parsimony_seed = self.config.parsimony_seed
            self.logger.info(f"🎯 使用指定parsimony种子 | Using specified parsimony seed: {parsimony_seed}")
        
        if self.config.bootstrap_seed is None and self.is_bootstrap_analysis():
            bootstrap_seed = random.randint(1, 1000000)
            self.logger.info(f"🎲 生成bootstrap随机种子 | Generated bootstrap random seed: {bootstrap_seed}")
        else:
            bootstrap_seed = self.config.bootstrap_seed
        
        return parsimony_seed, bootstrap_seed
    
    def is_bootstrap_analysis(self) -> bool:
        """检查是否是bootstrap分析 | Check if bootstrap analysis"""
        return (self.config.bootstrap_seed is not None or 
                self.config.rapid_bootstrap_seed is not None or
                self.runs != "1" and self.runs not in ["autoFC", "autoMR", "autoMRE", "autoMRE_IGN"])
    
    def build_raxml_command(self) -> str:
        """构建RAxML命令 | Build RAxML command"""
        self.logger.info("🔧 构建RAxML命令 | Building RAxML command")
        
        # 基础命令 | Basic command
        cmd = [
            self.config.raxml_path,
            f"-s {self.config.sequence_file}",
            f"-n {self.config.output_name}",
            f"-m {self.config.model}",
            f"-w {self.config.output_path.resolve()}"
        ]
        
        # 算法选择 | Algorithm selection
        cmd.append(f"-f {self.config.algorithm}")
        
        # 添加种子 | Add seeds
        parsimony_seed, bootstrap_seed = self.generate_random_seeds()
        cmd.append(f"-p {parsimony_seed}")
        
        # Bootstrap参数 | Bootstrap parameters
        if bootstrap_seed:
            cmd.append(f"-b {bootstrap_seed}")
        
        if self.config.rapid_bootstrap_seed:
            cmd.append(f"-x {self.config.rapid_bootstrap_seed}")
        
        # 运行次数 | Number of runs
        cmd.append(f"-# {self.config.runs}")
        
        # 线程数 | Number of threads
        if "PTHREADS" in self.config.raxml_path.upper():
            cmd.append(f"-T {self.config.threads}")
        
        # 树参数 | Tree parameters
        if self.config.starting_tree:
            cmd.append(f"-t {self.config.starting_tree}")
        
        if self.config.constraint_tree:
            cmd.append(f"-g {self.config.constraint_tree}")
        
        if self.config.outgroup:
            cmd.append(f"-o {self.config.outgroup}")
        
        # 模型参数 | Model parameters
        cmd.append(f"-c {self.config.categories}")
        
        if not self.config.rate_het_model:
            cmd.append("-V")
        
        if self.config.gamma_median:
            cmd.append("-u")
        
        # Bootstrap参数 | Bootstrap parameters
        if self.config.bootstrap_convergence:
            cmd.append(f"-I {self.config.bootstrap_convergence}")
        
        cmd.append(f"-B {self.config.bootstop_threshold}")
        cmd.append(f"--bootstop-perms={self.config.bootstop_perms}")
        
        if self.config.print_bootstrap_trees:
            cmd.append("-k")
        
        # 性能优化参数 | Performance optimization parameters
        cmd.append(f"-e {self.config.likelihood_epsilon}")
        
        if self.config.ml_search_convergence:
            cmd.append("-D")
        
        if self.config.random_starting_tree:
            cmd.append("-d")
        
        if self.config.memory_saving:
            cmd.append("-U")
        
        if not self.config.pattern_compression:
            cmd.append("-H")
        
        # 质量控制参数 | Quality control parameters
        if self.config.no_seq_check:
            cmd.append("--no-seq-check")
        
        if self.config.silent:
            cmd.append("--silent")
        
        return " ".join(cmd)
    
    def run_raxml_analysis(self) -> bool:
        """运行RAxML分析 | Run RAxML analysis"""
        self.logger.info("🌳 开始RAxML系统发育分析 | Starting RAxML phylogenetic analysis")
        
        # 构建命令 | Build command
        cmd = self.build_raxml_command()
        
        # 估算运行时间 | Estimate runtime
        self.estimate_runtime()
        
        # 执行分析 | Execute analysis
        success = self.cmd_runner.run(
            cmd, 
            "RAxML系统发育分析 | RAxML phylogenetic analysis",
            timeout=None  # 不设置超时，让长时间分析可以完成 | No timeout for long analyses
        )
        
        if success:
            self.log_analysis_results()
        
        return success
    
    def estimate_runtime(self):
        """估算运行时间 | Estimate runtime"""
        self.logger.info("⏱️ 运行时间估算 | Runtime estimation:")
        
        # 简单的运行时间估算 | Simple runtime estimation
        try:
            # 尝试获取序列数量 | Try to get sequence count
            with open(self.config.sequence_file, 'r') as f:
                first_line = f.readline().strip()
                if first_line:
                    parts = first_line.split()
                    if len(parts) >= 2:
                        num_taxa = int(parts[0])
                        seq_length = int(parts[1])
                        
                        self.logger.info(f"📊 序列信息 | Sequence info: {num_taxa} 条序列 | sequences, {seq_length} bp")
                        
                        # 粗略估算 | Rough estimation
                        if num_taxa < 50:
                            estimate = "几分钟到几小时 | Several minutes to hours"
                        elif num_taxa < 200:
                            estimate = "几小时到一天 | Several hours to a day"
                        else:
                            estimate = "数天到数周 | Several days to weeks"
                        
                        self.logger.info(f"⏰ 预计运行时间 | Estimated runtime: {estimate}")
                        
                        if self.is_bootstrap_analysis():
                            self.logger.info("🔄 注意：Bootstrap分析会显著增加运行时间 | Note: Bootstrap analysis will significantly increase runtime")
        
        except Exception as e:
            self.logger.warning(f"⚠️ 无法估算运行时间 | Cannot estimate runtime: {e}")
    
    def log_analysis_results(self):
        """记录分析结果 | Log analysis results"""
        self.logger.info("📋 分析结果文件 | Analysis result files:")
        
        # 常见的RAxML输出文件 | Common RAxML output files
        result_patterns = [
            f"RAxML_bestTree.{self.config.output_name}",
            f"RAxML_bipartitions.{self.config.output_name}",
            f"RAxML_bipartitionsBranchLabels.{self.config.output_name}",
            f"RAxML_bootstrap.{self.config.output_name}",
            f"RAxML_info.{self.config.output_name}",
            f"RAxML_log.{self.config.output_name}",
            f"RAxML_parsimonyTree.{self.config.output_name}",
            f"RAxML_result.{self.config.output_name}"
        ]
        
        found_files = []
        for pattern in result_patterns:
            file_path = self.config.output_path / pattern
            if file_path.exists():
                size = file_path.stat().st_size
                found_files.append(f"  📄 {pattern} ({size} bytes)")
        
        if found_files:
            self.logger.info("🎯 生成的文件 | Generated files:")
            for file_info in found_files:
                self.logger.info(file_info)
        else:
            self.logger.warning("⚠️ 未找到预期的输出文件 | No expected output files found")
