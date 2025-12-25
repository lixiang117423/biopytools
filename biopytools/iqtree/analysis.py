"""
🧬 IQ-TREE 核心分析模块 | IQ-TREE Core Analysis Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class TreeBuilder:
    """🌳 系统发育树构建器 | Phylogenetic Tree Builder"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_command(self) -> str:
        """构建IQ-TREE命令 | Build IQ-TREE command"""
        cmd_parts = [self.config.iqtree_path]
        
        # 输入文件 | Input file
        cmd_parts.append(f"-s {self.config.input_file}")
        
        # 输出前缀 | Output prefix
        output_prefix = self.config.output_path / self.config.prefix
        cmd_parts.append(f"--prefix {output_prefix}")
        
        # 模型选择 | Model selection
        if self.config.model:
            cmd_parts.append(f"-m {self.config.model}")
        else:
            # 自动模型选择 | Automatic model selection
            cmd_parts.append("-m MFP")
            self.logger.info("🔬 使用自动模型选择 | Using automatic model selection (MFP)")
        
        # Bootstrap设置 | Bootstrap settings
        if self.config.bootstrap > 0:
            if self.config.bootstrap_type == 'ufboot':
                cmd_parts.append(f"-B {self.config.bootstrap}")
                self.logger.info(f"🔄 使用UFBoot，重复{self.config.bootstrap}次 | Using UFBoot with {self.config.bootstrap} replicates")
                # 如果需要保存bootstrap树 | Save bootstrap trees if requested
                if self.config.save_boot_trees:
                    cmd_parts.append("--wbtl")
                    self.logger.info(f"💾 保存所有bootstrap树 | Saving all bootstrap trees")
            else:
                cmd_parts.append(f"-b {self.config.bootstrap}")
                self.logger.info(f"🔄 使用标准Bootstrap，重复{self.config.bootstrap}次 | Using standard bootstrap with {self.config.bootstrap} replicates")
        
        # 线程数 | Threads
        cmd_parts.append(f"-T {self.config.threads}")
        
        # 外群设置 | Outgroup setting
        if self.config.outgroup:
            cmd_parts.append(f"-o {self.config.outgroup}")
            self.logger.info(f"🎯 设置外群 | Setting outgroup: {self.config.outgroup}")
        
        # 分区分析 | Partition analysis
        if self.config.partition_file:
            if self.config.partition_mode == 'edge-linked':
                cmd_parts.append(f"-p {self.config.partition_file}")
            elif self.config.partition_mode == 'edge-equal':
                cmd_parts.append(f"-q {self.config.partition_file}")
            elif self.config.partition_mode == 'edge-unlinked':
                cmd_parts.append(f"-Q {self.config.partition_file}")
            self.logger.info(f"📊 启用分区分析 | Partition analysis enabled: {self.config.partition_mode}")
        
        # 约束树 | Constraint tree
        if self.config.constraint_tree:
            cmd_parts.append(f"-g {self.config.constraint_tree}")
            self.logger.info(f"🔒 使用约束树 | Using constraint tree")
        
        # 其他参数 | Other parameters
        if self.config.seed:
            cmd_parts.append(f"--seed {self.config.seed}")
        
        if self.config.runs > 1:
            cmd_parts.append(f"--runs {self.config.runs}")
        
        if self.config.redo:
            cmd_parts.append("--redo")
        
        return " ".join(cmd_parts)
    
    def run_tree_inference(self) -> bool:
        """运行系统发育树推断 | Run phylogenetic tree inference"""
        self.logger.info("=" * 60)
        self.logger.info("🌳 开始系统发育树构建 | Starting phylogenetic tree construction")
        self.logger.info("=" * 60)
        
        cmd = self.build_command()
        success = self.cmd_runner.run(
            cmd,
            "系统发育树推断 | Phylogenetic tree inference"
        )
        
        if success:
            self.logger.info("✅ 系统发育树构建完成 | Tree construction completed")
            self._report_output_files()
            self._report_recommended_output()  # ← 添加这一行
                
        return success
    
    def _report_output_files(self):
        """报告输出文件 | Report output files"""
        self.logger.info("📁 输出文件 | Output files:")
        
        output_prefix = self.config.output_path / self.config.prefix
        expected_files = [
            (f"{output_prefix}.treefile", "最佳树文件 | Best tree file"),
            (f"{output_prefix}.iqtree", "主要结果文件 | Main result file"),
            (f"{output_prefix}.log", "运行日志 | Run log"),
        ]
        
        if self.config.bootstrap > 0:
            if self.config.bootstrap_type == 'ufboot' and self.config.save_boot_trees:
                expected_files.append((f"{output_prefix}.ufboot", "UFBoot树文件 | UFBoot trees"))
            expected_files.append((f"{output_prefix}.contree", "共识树 | Consensus tree"))
        
        if not self.config.model or self.config.model == 'MFP':
            expected_files.append((f"{output_prefix}.model.gz", "模型选择结果 | Model selection results"))
        
        for file_path, description in expected_files:
            if os.path.exists(file_path):
                self.logger.info(f"  ✓ {description}: {file_path}")
            else:
                self.logger.warning(f"  ⚠ {description}: {file_path} (未找到 | not found)")
    
    def _report_recommended_output(self):
        """推荐使用的输出文件 | Recommend output files to use"""
        self.logger.info("=" * 60)
        self.logger.info("📌 推荐使用的文件 | Recommended Output Files")
        self.logger.info("=" * 60)
        
        output_prefix = self.config.output_path / self.config.prefix
        
        if self.config.bootstrap > 0:
            self.logger.info(f"⭐ 主要结果（带bootstrap支持值）| Main result (with bootstrap):")
            self.logger.info(f"   👉 {output_prefix}.contree")
            self.logger.info(f"   💡 用途：可视化、发表、展示 | For visualization, publication")
        
        self.logger.info(f"📊 最佳ML树（精确branch length）| Best ML tree:")
        self.logger.info(f"   👉 {output_prefix}.treefile")
        self.logger.info(f"   💡 用途：进一步进化分析 | For further evolutionary analysis")
        
        self.logger.info(f"📄 详细结果文件 | Detailed results:")
        self.logger.info(f"   👉 {output_prefix}.iqtree")
        self.logger.info(f"   💡 包含模型参数、统计信息 | Contains model parameters, statistics")
        
        if self.config.enable_ancestral:
            self.logger.info(f"🧬 祖先状态重建结果 | Ancestral reconstruction:")
            self.logger.info(f"   👉 {output_prefix}_ancestral.treefile")
            self.logger.info(f"   👉 {output_prefix}_ancestral_sequences.fasta")
            self.logger.info(f"   💡 用途：祖先序列分析 | For ancestral sequence analysis")

class ConcordanceAnalyzer:
    """🔬 一致性因子分析器 | Concordance Factor Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_concordance_analysis(self) -> bool:
        """运行一致性因子分析 | Run concordance factor analysis"""
        if not self.config.enable_concordance:
            return True
        
        self.logger.info("=" * 60)
        self.logger.info("🔬 开始一致性因子分析 | Starting concordance factor analysis")
        self.logger.info("=" * 60)
        
        output_prefix = self.config.output_path / self.config.prefix
        tree_file = f"{output_prefix}.treefile"
        
        if not os.path.exists(tree_file):
            self.logger.error(f"❌ 未找到树文件，无法进行一致性因子分析 | Tree file not found: {tree_file}")
            return False
        
        # 构建gCF命令 | Build gCF command
        cf_prefix = self.config.output_path / f"{self.config.prefix}_concordance"
        cmd = f"{self.config.iqtree_path} -t {tree_file} --gcf {self.config.concordance_trees} --prefix {cf_prefix}"
        
        if self.config.threads:
            cmd += f" -T {self.config.threads}"
        
        success = self.cmd_runner.run(
            cmd,
            "基因一致性因子分析 | Gene concordance factor analysis"
        )
        
        if success:
            self.logger.info("✅ 一致性因子分析完成 | Concordance factor analysis completed")
        
        return success

class AncestralReconstructor:
    """🧬 祖先状态重建器 | Ancestral State Reconstructor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_ancestral_reconstruction(self) -> bool:
        """运行祖先状态重建 | Run ancestral state reconstruction"""
        if not self.config.enable_ancestral:
            return True
        
        self.logger.info("=" * 60)
        self.logger.info("🧬 开始祖先状态重建 | Starting ancestral state reconstruction")
        self.logger.info("=" * 60)
        
        output_prefix = self.config.output_path / self.config.prefix
        tree_file = f"{output_prefix}.treefile"
        
        if not os.path.exists(tree_file):
            self.logger.error(f"❌ 未找到树文件，无法进行祖先状态重建 | Tree file not found: {tree_file}")
            return False
        
        # 构建祖先重建命令 | Build ancestral reconstruction command
        asr_prefix = self.config.output_path / f"{self.config.prefix}_ancestral"
        cmd = f"{self.config.iqtree_path} -s {self.config.input_file} -t {tree_file} --ancestral --prefix {asr_prefix}"
        
        if self.config.model:
            cmd += f" -m {self.config.model}"
        
        if self.config.threads:
            cmd += f" -T {self.config.threads}"
        
        success = self.cmd_runner.run(
            cmd,
            "祖先状态重建 | Ancestral state reconstruction"
        )
        
        if success:
            self.logger.info("✅ 祖先状态重建完成 | Ancestral state reconstruction completed")
        
        return success
