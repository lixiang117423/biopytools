"""
IQ-TREE 核心分析模块|IQ-TREE Core Analysis Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, build_conda_command


class TreeBuilder:
    """系统发育树构建器|Phylogenetic Tree Builder"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_command(self) -> list:
        """构建IQ-TREE命令|Build IQ-TREE command"""
        args = []

        # 输入文件|Input file
        args.extend(['-s', self.config.input_file])

        # 输出前缀|Output prefix
        args.extend(['--prefix', self.config.prefix])

        # 模型选择|Model selection
        if self.config.model:
            args.extend(['-m', self.config.model])
        else:
            args.extend(['-m', 'MFP'])
            self.logger.info("使用自动模型选择|Using automatic model selection (MFP)")

        # Bootstrap设置|Bootstrap settings
        if self.config.bootstrap > 0:
            if self.config.bootstrap_type == 'ufboot':
                args.extend([f'-B', str(self.config.bootstrap)])
                self.logger.info(f"使用UFBoot，重复{self.config.bootstrap}次|Using UFBoot with {self.config.bootstrap} replicates")
                if self.config.save_boot_trees:
                    args.append('--wbtl')
                    self.logger.info("保存所有bootstrap树|Saving all bootstrap trees")
            else:
                args.extend([f'-b', str(self.config.bootstrap)])
                self.logger.info(f"使用标准Bootstrap，重复{self.config.bootstrap}次|Using standard bootstrap with {self.config.bootstrap} replicates")

        # 线程数|Threads
        args.extend(['-T', str(self.config.threads)])

        # 外群设置|Outgroup setting
        if self.config.outgroup:
            args.extend(['-o', self.config.outgroup])
            self.logger.info(f"设置外群|Setting outgroup: {self.config.outgroup}")

        # 分区分析|Partition analysis
        if self.config.partition_file:
            if self.config.partition_mode == 'edge-linked':
                args.extend(['-p', self.config.partition_file])
            elif self.config.partition_mode == 'edge-equal':
                args.extend(['-q', self.config.partition_file])
            elif self.config.partition_mode == 'edge-unlinked':
                args.extend(['-Q', self.config.partition_file])
            self.logger.info(f"启用分区分析|Partition analysis enabled: {self.config.partition_mode}")

        # 约束树|Constraint tree
        if self.config.constraint_tree:
            args.extend(['-g', self.config.constraint_tree])
            self.logger.info("使用约束树|Using constraint tree")

        # 其他参数|Other parameters
        if self.config.seed:
            args.extend(['--seed', str(self.config.seed)])

        if self.config.runs > 1:
            args.extend(['--runs', str(self.config.runs)])

        if self.config.redo:
            args.append('--redo')

        return build_conda_command(self.config.iqtree_path, args)

    def run_tree_inference(self) -> bool:
        """运行系统发育树推断|Run phylogenetic tree inference"""
        self.logger.info("=" * 60)
        self.logger.info("开始系统发育树构建|Starting phylogenetic tree construction")
        self.logger.info("=" * 60)

        output_prefix = self.config.output_path / self.config.prefix
        tree_file = f"{output_prefix}.treefile"

        if self._is_step_completed(tree_file) and not self.config.redo:
            self.logger.info("跳过已完成步骤|Skipping completed step: 系统发育树推断|phylogenetic tree inference")
            self._report_output_files()
            self._report_recommended_output()
            return True

        cmd = self.build_command()
        success = self.cmd_runner.run(
            cmd,
            "系统发育树推断|Phylogenetic tree inference"
        )

        if success:
            self.logger.info("系统发育树构建完成|Tree construction completed")
            self._report_output_files()
            self._report_recommended_output()

        return success

    def _is_step_completed(self, output_file: str) -> bool:
        """检查步骤是否已完成|Check if step is completed"""
        return Path(output_file).exists()

    def _report_output_files(self):
        """报告输出文件|Report output files"""
        self.logger.info("输出文件|Output files:")

        output_prefix = self.config.output_path / self.config.prefix
        expected_files = [
            (f"{output_prefix}.treefile", "最佳树文件|Best tree file"),
            (f"{output_prefix}.iqtree", "主要结果文件|Main result file"),
            (f"{output_prefix}.log", "运行日志|Run log"),
        ]

        if self.config.bootstrap > 0:
            if self.config.bootstrap_type == 'ufboot' and self.config.save_boot_trees:
                expected_files.append((f"{output_prefix}.ufboot", "UFBoot树文件|UFBoot trees"))
            expected_files.append((f"{output_prefix}.contree", "共识树|Consensus tree"))

        if not self.config.model or self.config.model == 'MFP':
            expected_files.append((f"{output_prefix}.model.gz", "模型选择结果|Model selection results"))

        for file_path, description in expected_files:
            if os.path.exists(file_path):
                self.logger.info(f"  {description}: {file_path}")
            else:
                self.logger.warning(f"  {description}: {file_path} (未找到|not found)")

    def _report_recommended_output(self):
        """推荐使用的输出文件|Recommend output files to use"""
        self.logger.info("=" * 60)
        self.logger.info("推荐使用的文件|Recommended Output Files")
        self.logger.info("=" * 60)

        output_prefix = self.config.output_path / self.config.prefix

        if self.config.bootstrap > 0:
            self.logger.info("主要结果（带bootstrap支持值）|Main result (with bootstrap):")
            self.logger.info(f"   {output_prefix}.contree")
            self.logger.info("   用途：可视化、发表、展示|For visualization, publication")

        self.logger.info("最佳ML树（精确branch length）|Best ML tree:")
        self.logger.info(f"   {output_prefix}.treefile")
        self.logger.info("   用途：进一步进化分析|For further evolutionary analysis")

        self.logger.info("详细结果文件|Detailed results:")
        self.logger.info(f"   {output_prefix}.iqtree")
        self.logger.info("   包含模型参数、统计信息|Contains model parameters, statistics")

        if self.config.enable_ancestral:
            self.logger.info("祖先状态重建结果|Ancestral reconstruction:")
            self.logger.info(f"   {output_prefix}_ancestral.treefile")
            self.logger.info(f"   {output_prefix}_ancestral_sequences.fasta")
            self.logger.info("   用途：祖先序列分析|For ancestral sequence analysis")


class ConcordanceAnalyzer:
    """一致性因子分析器|Concordance Factor Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_concordance_analysis(self) -> bool:
        """运行一致性因子分析|Run concordance factor analysis"""
        if not self.config.enable_concordance:
            return True

        self.logger.info("=" * 60)
        self.logger.info("开始一致性因子分析|Starting concordance factor analysis")
        self.logger.info("=" * 60)

        output_prefix = self.config.output_path / self.config.prefix
        tree_file = f"{output_prefix}.treefile"

        if not os.path.exists(tree_file):
            self.logger.error(f"未找到树文件，无法进行一致性因子分析|Tree file not found: {tree_file}")
            return False

        cf_prefix = f"{self.config.prefix}_concordance"
        cf_output = f"{self.config.output_path / cf_prefix}.cf.tree"

        if self._is_step_completed(cf_output) and not self.config.redo:
            self.logger.info("跳过已完成步骤|Skipping completed step: 一致性因子分析|concordance factor analysis")
            return True

        args = ['-t', tree_file, '--gcf', str(self.config.concordance_trees),
                '--prefix', cf_prefix]
        if self.config.threads:
            args.extend(['-T', str(self.config.threads)])

        cmd = build_conda_command(self.config.iqtree_path, args)

        success = self.cmd_runner.run(
            cmd,
            "基因一致性因子分析|Gene concordance factor analysis"
        )

        if success:
            self.logger.info("一致性因子分析完成|Concordance factor analysis completed")

        return success

    def _is_step_completed(self, output_file: str) -> bool:
        return Path(output_file).exists()


class AncestralReconstructor:
    """祖先状态重建器|Ancestral State Reconstructor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_ancestral_reconstruction(self) -> bool:
        """运行祖先状态重建|Run ancestral state reconstruction"""
        if not self.config.enable_ancestral:
            return True

        self.logger.info("=" * 60)
        self.logger.info("开始祖先状态重建|Starting ancestral state reconstruction")
        self.logger.info("=" * 60)

        output_prefix = self.config.output_path / self.config.prefix
        tree_file = f"{output_prefix}.treefile"

        if not os.path.exists(tree_file):
            self.logger.error(f"未找到树文件，无法进行祖先状态重建|Tree file not found: {tree_file}")
            return False

        asr_prefix = f"{self.config.prefix}_ancestral"
        asr_output = f"{self.config.output_path / asr_prefix}.state"

        if self._is_step_completed(asr_output) and not self.config.redo:
            self.logger.info("跳过已完成步骤|Skipping completed step: 祖先状态重建|ancestral state reconstruction")
            return True

        args = ['-s', self.config.input_file, '-t', tree_file,
                '--ancestral', '--prefix', asr_prefix]
        if self.config.model:
            args.extend(['-m', self.config.model])
        if self.config.threads:
            args.extend(['-T', str(self.config.threads)])

        cmd = build_conda_command(self.config.iqtree_path, args)

        success = self.cmd_runner.run(
            cmd,
            "祖先状态重建|Ancestral state reconstruction"
        )

        if success:
            self.logger.info("祖先状态重建完成|Ancestral state reconstruction completed")

        return success

    def _is_step_completed(self, output_file: str) -> bool:
        return Path(output_file).exists()
