"""
HapHiC流程管理模块 - Pipeline模式 | HapHiC Pipeline Management Module - Pipeline Mode
"""

import os
import sys
import time
import glob
from pathlib import Path
from typing import Dict, List, Optional
import shutil
import subprocess

from .config import HapHiCConfig
from .utils import HapHiCLogger, CommandRunner, FileManager


class HapHiCPipeline:
    """HapHiC流程管理器 - Pipeline模式 | HapHiC Pipeline Manager - Pipeline Mode"""

    def __init__(self, config: HapHiCConfig, logger):
        """初始化流程管理器 | Initialize pipeline manager"""
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger, config.dry_run)
        self.step_dirs = {}

    def run_complete_pipeline(self) -> bool:
        """运行完整的HapHiC pipeline | Run complete HapHiC pipeline"""
        try:
            self.logger.info("🚀 开始HapHiC Pipeline流程 | Starting HapHiC Pipeline")
            start_time = time.time()

            # 创建输出目录结构
            self._create_directory_structure()

            # 如果强制重新运行，清理所有目录
            if self.config.force_rerun:
                self.logger.info("🔄 强制重新运行模式，清理所有目录 | Force rerun mode, cleaning all directories")
                self._force_cleanup_all_directories()

            # 检查HapHiC pipeline是否已完成
            if not self.config.force_rerun and self._check_haphic_pipeline_completed():
                self.logger.info("✅ HapHiC Pipeline已完成，跳过执行 | HapHiC Pipeline already completed, skipping execution")
            else:
                # 运行HapHiC pipeline（一步完成）
                if not self._run_haphic_pipeline():
                    return False

            # 生成可视化结果（默认执行）
            if not self.config.force_rerun and self._check_visualization_completed():
                self.logger.info("✅ 可视化结果已存在，跳过生成 | Visualization results already exist, skipping generation")
            else:
                if self.config.force_rerun:
                    self.logger.info("🔄 强制重新运行模式，重新生成可视化 | Force rerun mode, regenerating visualization")
                # 如果禁用可视化，直接跳过
                if not self.config.generate_plots:
                    self.logger.info("📊 可视化生成已禁用 | Visualization generation disabled")
                else:
                    # 尝试生成可视化，如果失败则警告但继续执行
                    if not self._generate_visualization():
                        self.logger.warning("⚠️ 可视化生成失败，但继续执行其他步骤 | Visualization failed, but continuing with other steps")

            # 生成Juicebox文件（如果需要）
            if self.config.output_juicebox:
                if not self.config.force_rerun and self._check_juicebox_files_completed():
                    self.logger.info("✅ Juicebox文件已存在，跳过生成 | Juicebox files already exist, skipping generation")
                else:
                    if self.config.force_rerun:
                        self.logger.info("🔄 强制重新运行模式，重新生成Juicebox文件 | Force rerun mode, regenerating Juicebox files")
                    # 尝试生成Juicebox文件，如果失败则警告但继续执行
                    if not self._generate_juicebox_files():
                        self.logger.warning("⚠️ Juicebox文件生成失败，但继续执行其他步骤 | Juicebox file generation failed, but continuing with other steps")

            # 统计信息
            elapsed_time = time.time() - start_time
            self.logger.info(f"✅ HapHiC Pipeline完成 | HapHiC Pipeline completed in {elapsed_time:.2f}秒")

            return True

        except Exception as e:
            self.logger.error(f"💥 HapHiC Pipeline执行失败 | HapHiC Pipeline execution failed: {e}")
            return False

    def run_single_step(self, step_name: str) -> bool:
        """运行单个步骤 - 保持兼容性 | Run single step - compatibility mode"""
        self.logger.warning(f"⚠️ 单步运行模式已弃用，建议使用完整的pipeline模式 | Single-step mode deprecated, recommend using pipeline mode")
        return self.run_complete_pipeline()

    def continue_from_step(self, step_name: str) -> bool:
        """从指定步骤继续 - 保持兼容性 | Continue from specified step - compatibility mode"""
        self.logger.warning(f"⚠️ 断点续传模式已弃用，建议使用完整的pipeline模式 | Continue-from mode deprecated, recommend using pipeline mode")
        return self.run_complete_pipeline()

    def _create_directory_structure(self):
        """创建目录结构 | Create directory structure"""
        self.logger.info("📁 设置目录结构 | Setting up directory structure")

        # 创建主输出目录（如果不存在）
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 设置步骤目录路径（不创建，让HapHiC pipeline自己处理）
        self.step_dirs = {
            "cluster": os.path.join(self.config.output_dir, "01.cluster"),
            "reassign": os.path.join(self.config.output_dir, "02.reassign"),
            "sort": os.path.join(self.config.output_dir, "03.sort"),
            "build": os.path.join(self.config.output_dir, "04.build"),
            "plots": os.path.join(self.config.output_dir, "05.plots"),
            "juicebox": os.path.join(self.config.output_dir, "06.juicebox")
        }

        # 只创建juicebox目录，plots目录会在需要时创建
        os.makedirs(self.step_dirs["juicebox"], exist_ok=True)

        # plots目录延迟到实际生成可视化时创建

        self.logger.info("✅ 目录结构设置完成 | Directory structure setup complete")

    def _run_haphic_pipeline(self) -> bool:
        """运行HapHiC pipeline | Run HapHiC pipeline"""
        self.logger.info("🔄 运行HapHiC Pipeline | Running HapHiC Pipeline")

        # 注意：不再清理已存在的目录，HapHiC可以处理已存在的目录
        self.logger.info("ℹ️ HapHiC将处理已存在的目录 | HapHiC will handle existing directories")

        # 构建pipeline命令
        cmd = [self.config.haphic_bin, "pipeline"]
        cmd.extend([self.config.asm_file, self.config.bam_file, str(self.config.nchrs)])

        # 添加所有pipeline参数
        cmd.extend(self._get_pipeline_options())

        # 设置工作目录为输出目录（避免目录冲突）
        original_cwd = os.getcwd()
        try:
            os.chdir(self.config.output_dir)
            self.logger.info(f"📂 工作目录: {self.config.output_dir}")

            success = self.cmd_runner.run_command(cmd, "HapHiC Pipeline")

            # Pipeline完成后，检查输出文件
            if success:
                self._verify_pipeline_output()

            return success

        finally:
            os.chdir(original_cwd)

    def _check_haphic_pipeline_completed(self) -> bool:
        """检查HapHiC pipeline是否已完成 | Check if HapHiC pipeline is completed"""
        self.logger.info("🔍 检查HapHiC Pipeline完成状态 | Checking HapHiC Pipeline completion status")

        # 检查各步骤的特征文件
        checks = []

        # 1. Cluster步骤检查
        cluster_dir = self.step_dirs["cluster"]
        corrected_asm = os.path.join(cluster_dir, "corrected_asm.fa")
        paired_links_clm = os.path.join(cluster_dir, "paired_links.clm")
        cluster_completed = os.path.exists(corrected_asm) and os.path.exists(paired_links_clm)
        checks.append(("Cluster", cluster_completed, [corrected_asm, paired_links_clm]))

        # 2. Reassign步骤检查
        reassign_dir = self.step_dirs["reassign"]
        final_clusters = os.path.join(reassign_dir, "final_groups", "final_clusters.txt")
        reassign_completed = os.path.exists(final_clusters)
        checks.append(("Reassign", reassign_completed, [final_clusters]))

        # 3. Sort步骤检查
        sort_dir = self.step_dirs["sort"]
        # 检查是否有至少12个tour文件（对应12个染色体）
        tour_files = []
        if os.path.exists(sort_dir):
            tour_files = [f for f in os.listdir(sort_dir) if f.endswith('.tour')]
        sort_completed = len(tour_files) >= self.config.nchrs
        checks.append(("Sort", sort_completed, [f"{sort_dir}/*.tour ({len(tour_files)} files)"]))

        # 4. Build步骤检查
        build_dir = self.step_dirs["build"]
        scaffolds_fa = os.path.join(build_dir, f"{self.config.prefix}.fa")
        scaffolds_agp = os.path.join(build_dir, f"{self.config.prefix}.agp")
        build_completed = os.path.exists(scaffolds_fa) and os.path.exists(scaffolds_agp)
        checks.append(("Build", build_completed, [scaffolds_fa, scaffolds_agp]))

        # 记录检查结果
        all_completed = True
        for step_name, completed, files in checks:
            if completed:
                self.logger.info(f"✅ {step_name}步骤已完成 | {step_name} step completed")
            else:
                self.logger.info(f"❌ {step_name}步骤未完成 | {step_name} step not completed")
                self.logger.info(f"   期待文件: {', '.join(files)}")
                all_completed = False

        return all_completed

    def _force_cleanup_all_directories(self):
        """强制清理所有步骤目录 | Force clean all step directories"""
        all_dirs = ["01.cluster", "02.reassign", "03.sort", "04.build", "05.plots", "06.juicebox"]

        for dir_name in all_dirs:
            dir_path = os.path.join(self.config.output_dir, dir_name)
            if os.path.exists(dir_path):
                try:
                    shutil.rmtree(dir_path)
                    self.logger.info(f"🗑️ 强制清理目录: {dir_name}")
                except Exception as e:
                    self.logger.warning(f"⚠️ 无法清理目录 {dir_name}: {e}")

        # 重新创建必要的目录
        os.makedirs(self.step_dirs["juicebox"], exist_ok=True)

    def _cleanup_existing_directories(self):
        """清理可能存在的HapHiC子目录 - 仅在需要时执行 | Clean up existing HapHiC subdirectories - only execute when needed"""
        # 注意：由于实现了断点续传，默认不再清理已存在的目录
        # 只清理非HapHiC pipeline管理的目录（如05.plots）

        # 只清理plots目录，因为它不由HapHiC pipeline管理
        plots_dir = self.step_dirs["plots"]
        if os.path.exists(plots_dir):
            try:
                shutil.rmtree(plots_dir)
                if self.logger:
                    self.logger.info(f"🗑️ 清理plots目录以重新生成可视化 | Cleaning plots directory to regenerate visualization")
                # 重新创建空目录
                os.makedirs(plots_dir, exist_ok=True)
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"⚠️ 无法清理plots目录: {e}")

        # HapHiC pipeline管理的目录（01.cluster到04.build）不再清理，支持断点续传
        if self.logger:
            self.logger.info("ℹ️ HapHiC pipeline目录保留以支持断点续传 | HapHiC pipeline directories preserved for resume support")

    def _get_pipeline_options(self) -> List[str]:
        """获取pipeline命令选项 | Get pipeline command options"""
        options = []

        # 输入文件选项
        options.extend([
            "--aln_format", "bam",
            "--outdir", self.config.output_dir
        ])

        # 基本pipeline控制选项
        options.extend([
            "--RE", self.config.re_sites,
            "--steps", "1,2,3,4",  # 运行所有步骤
            "--threads", str(self.config.threads),
            "--processes", str(self.config.processes)
        ])

        # 组装校正选项
        if self.config.correct_nrounds > 0:
            # 使用更高的分辨率来减少内存使用
            if self.config.correct_resolution is not None:
                # 使用用户指定的分辨率
                correct_resolution = self.config.correct_resolution
                if self.logger:
                    self.logger.info(f"🔧 使用用户指定的分辨率: {correct_resolution}bp")
            else:
                # 自动优化：根据覆盖度选择分辨率
                correct_resolution = 1000 if self.config.correct_min_coverage <= 5 else 500
                if self.logger:
                    self.logger.info(f"🔧 自动优化分辨率: {correct_resolution}bp (基于覆盖率{self.config.correct_min_coverage})")
            options.extend([
                "--correct_nrounds", str(self.config.correct_nrounds),
                "--correct_resolution", str(correct_resolution),
                "--median_cov_ratio", "0.2",
                "--region_len_ratio", "0.1",
                "--min_region_cutoff", "5000"
            ])
            if self.logger:
                self.logger.info(f"🔧 启用组装校正，使用分辨率: {correct_resolution}bp (内存优化)")

        # 过滤参数
        options.extend([
            "--Nx", str(self.config.nx),
            "--RE_site_cutoff", str(self.config.min_re_sites),
            "--density_lower", "0.2X",
            "--density_upper", "1.9X",
            "--read_depth_upper", "1.5X",
            "--topN", "10"
        ])

        # MCL聚类参数
        options.extend([
            "--min_inflation", str(self.config.min_inflation),
            "--max_inflation", str(self.config.max_inflation),
            "--inflation_step", str(self.config.inflation_step),
            "--expansion", "2",
            "--max_iter", "200",
            "--pruning", "0.0001"
        ])

        # 重新分配参数
        options.extend([
            "--min_group_len", str(self.config.min_group_len),
            "--max_ctg_len", "10000",  # 10Mbp = 10000Kbp
            "--min_RE_sites", str(self.config.min_re_sites),
            "--min_links", "25",  # HapHiC pipeline默认值
            "--min_link_density", "0.0001",  # HapHiC pipeline默认值
            "--min_density_ratio", "4",
            "--ambiguous_cutoff", "0.6",
            "--reassign_nrounds", "5"
        ])

        # 排序参数
        if not self.config.fast_sorting:
            options.append("--skip_fast_sort")

        if not self.config.allhic_optimization:
            options.append("--skip_allhic")

        # ALLHiC优化参数
        options.extend([
            "--mutprob", "0.2",
            "--ngen", "5000",
            "--npop", "100",
            "--seed", "42"
        ])

        # 构建参数
        options.extend([
            "--Ns", "100",
            "--max_width", "60",
            "--prefix", self.config.prefix
        ])

        if self.config.keep_letter_case:
            options.append("--keep_letter_case")

        # 等位基因处理选项
        if self.config.remove_allelic_links:
            options.extend([
                "--remove_allelic_links", str(self.config.remove_allelic_links),
                "--concordance_ratio_cutoff", "0.2",
                "--nwindows", "50",
                "--max_read_pairs", "200",
                "--min_read_pairs", "20"
            ])

        if self.config.gfa_files:
            options.extend(["--gfa", self.config.gfa_files])
            options.extend(["--phasing_weight", str(self.config.phasing_weight)])

        # 性能参数已在前面添加，这里不需要重复

        # 日志选项
        if self.config.verbose:
            options.append("--verbose")

        return options

    def _verify_pipeline_output(self):
        """验证pipeline输出 | Verify pipeline output"""
        self.logger.info("🔍 验证Pipeline输出 | Verifying pipeline output")

        output_files = self.config.get_output_files()

        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"✅ {file_type}: {os.path.basename(file_path)} ({size:,} bytes)")
            else:
                self.logger.warning(f"⚠️ {file_type}不存在: {file_path}")

    def _generate_juicebox_files(self) -> bool:
        """生成Juicebox兼容文件 | Generate Juicebox compatible files"""
        self.logger.info("🥤 生成Juicebox文件 | Generating Juicebox files")

        try:
            # 确保juicebox目录存在
            juicebox_dir = self.step_dirs["juicebox"]
            os.makedirs(juicebox_dir, exist_ok=True)
            self.logger.info(f"📂 Juicebox目录: {juicebox_dir}")

            build_dir = self.step_dirs["build"]

            # 查找必要的文件 - 修正为HapHiC实际生成的位置
            # 优先使用raw.agp文件，格式更标准
            agp_file = os.path.join(build_dir, f"{self.config.prefix}.raw.agp")
            if not os.path.exists(agp_file):
                agp_file = os.path.join(build_dir, f"{self.config.prefix}.agp")
                self.logger.info(f"🥤 使用标准AGP文件: {agp_file}")
            else:
                self.logger.info(f"🥤 使用Raw AGP文件（格式更标准）: {agp_file}")

            bam_file = self.config.bam_file

            if not FileManager.check_file_exists(agp_file):
                self.logger.error(f"❌ AGP文件不存在: {agp_file}")
                return False

            if not FileManager.check_file_exists(bam_file):
                self.logger.error(f"❌ BAM文件不存在: {bam_file}")
                return False

            # 1. 使用matlock生成mnd文件
            self.logger.info("🔧 步骤1: 使用matlock生成.mnd文件 | Step 1: Generate .mnd file with matlock")
            mnd_file = os.path.join(juicebox_dir, "out.links.mnd")
            cmd1 = ["matlock", "bam2", "juicer", bam_file, mnd_file]

            if not self.cmd_runner.run_command(cmd1, "生成.mnd文件"):
                return False

            # 2. 排序mnd文件
            self.logger.info("🔧 步骤2: 排序.mnd文件 | Step 2: Sorting .mnd file")
            sorted_mnd_file = os.path.join(juicebox_dir, "out.sorted.links.mnd")
            cmd2 = ["sort", "-k2,2", "-k6,6", mnd_file]

            try:
                with open(sorted_mnd_file, 'w') as outfile:
                    result = subprocess.run(cmd2, stdout=outfile, stderr=subprocess.PIPE, text=True)
                    if result.returncode != 0:
                        self.logger.error(f"❌ 排序.mnd文件失败: {result.stderr}")
                        return False
            except Exception as e:
                self.logger.error(f"❌ 排序.mnd文件异常: {e}")
                return False

            # 3. 生成assembly文件
            self.logger.info("🔧 步骤3: 生成.assembly文件 | Step 3: Generate .assembly file")
            assembly_file = os.path.join(juicebox_dir, f"{self.config.prefix}.assembly")
            cmd3 = ["python3", "~/software/3d-dna/utils/agp2assembly.py", agp_file, assembly_file]

            # 展开波浪号
            cmd3[1] = os.path.expanduser(cmd3[1])

            if not self.cmd_runner.run_command(cmd3, "生成assembly文件"):
                return False

            # 4. 使用asm-visualizer生成.hic文件
            self.logger.info("🔧 步骤4: 生成.hic文件 | Step 4: Generate .hic file")
            hic_file = os.path.join(juicebox_dir, f"{self.config.prefix}.hic")
            cmd4 = ["bash", "~/software/3d-dna/visualize/run-asm-visualizer.sh", "-p", "false", assembly_file, sorted_mnd_file]

            # 展开波浪号
            cmd4[0] = os.path.expanduser(cmd4[0])

            # 设置工作目录为juicebox目录
            original_cwd = os.getcwd()
            try:
                os.chdir(juicebox_dir)

                if not self.cmd_runner.run_command(cmd4, "生成.hic文件"):
                    return False

                # asm-visualizer可能会生成不同名称的.hic文件，需要重命名
                expected_hic = assembly_file.replace(".assembly", ".hic")
                if os.path.exists(expected_hic) and not os.path.exists(hic_file):
                    os.rename(expected_hic, hic_file)
                elif os.path.exists(f"{self.config.prefix}.hic") and not os.path.exists(hic_file):
                    os.rename(f"{self.config.prefix}.hic", hic_file)
                elif os.path.exists("assembly.hic") and not os.path.exists(hic_file):
                    if os.path.exists(hic_file):
                        os.remove(hic_file)
                    os.rename("assembly.hic", hic_file)

            finally:
                os.chdir(original_cwd)

            # 检查最终文件
            if FileManager.check_file_exists(hic_file):
                size = os.path.getsize(hic_file)
                self.logger.info(f"✅ .hic文件生成成功: {hic_file} ({size:,} bytes)")
            else:
                self.logger.error(f"❌ .hic文件生成失败: {hic_file}")
                return False

            if FileManager.check_file_exists(assembly_file):
                size = os.path.getsize(assembly_file)
                self.logger.info(f"✅ .assembly文件生成成功: {assembly_file} ({size:,} bytes)")
            else:
                self.logger.error(f"❌ .assembly文件生成失败: {assembly_file}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"💥 Juicebox文件生成失败: {e}")
            return False

    def get_step_status(self) -> Dict[str, bool]:
        """获取各步骤状态 | Get status of each step - 兼容性方法"""
        # Pipeline模式下检查最终输出
        output_files = self.config.get_output_files()

        # 简化的状态检查 - 主要看scaffolds文件是否存在
        has_scaffolds = FileManager.check_file_exists(output_files["scaffolds_fa"])

        return {
            "cluster": has_scaffolds,
            "reassign": has_scaffolds,
            "sort": has_scaffolds,
            "build": has_scaffolds
        }

    def _check_visualization_completed(self) -> bool:
        """检查可视化是否已完成 | Check if visualization is completed"""
        plots_dir = self.step_dirs["plots"]

        # 检查是否有PDF或PNG文件
        if not os.path.exists(plots_dir):
            return False

        plot_files = []
        for ext in ['*.pdf', '*.png', '*.jpg', '*.jpeg']:
            plot_files.extend(glob.glob(os.path.join(plots_dir, ext)))

        if plot_files:
            self.logger.info(f"📊 找到 {len(plot_files)} 个可视化文件 | Found {len(plot_files)} visualization files")
            return True
        return False

    def _check_juicebox_files_completed(self) -> bool:
        """检查Juicebox文件是否已完成 | Check if Juicebox files are completed"""
        juicebox_dir = self.step_dirs["juicebox"]

        # 检查关键文件
        hic_file = os.path.join(juicebox_dir, f"{self.config.prefix}.hic")
        assembly_file = os.path.join(juicebox_dir, f"{self.config.prefix}.assembly")

        if os.path.exists(hic_file) and os.path.exists(assembly_file):
            return True
        return False

    def _generate_visualization(self) -> bool:
        """生成可视化结果 | Generate visualization results"""
        self.logger.info("📊 生成Hi-C接触图可视化 | Generating Hi-C contact map visualization")

        try:
            plots_dir = self.step_dirs["plots"]

            # 创建plots目录（如果不存在）
            os.makedirs(plots_dir, exist_ok=True)

            # 使用HapHiC实际生成的AGP文件（优先使用raw.agp，格式更标准）
            build_dir = self.step_dirs["build"]
            agp_file = os.path.join(build_dir, f"{self.config.prefix}.raw.agp")

            # 如果raw.agp不存在，使用标准agp文件
            if not os.path.exists(agp_file):
                agp_file = os.path.join(build_dir, f"{self.config.prefix}.agp")
                self.logger.info(f"📊 使用标准AGP文件: {agp_file}")
            else:
                self.logger.info(f"📊 使用Raw AGP文件（格式更标准）: {agp_file}")

            if not FileManager.check_file_exists(agp_file):
                self.logger.warning(f"⚠️ 未找到AGP文件: {agp_file}，跳过可视化")
                return True

            # 预处理AGP文件，确保格式正确
            processed_agp_file = self._preprocess_agp_file(agp_file)
            if not processed_agp_file:
                self.logger.error("❌ AGP文件预处理失败")
                return False

            cmd = [self.config.haphic_bin, "plot"]
            cmd.extend([self.config.asm_file, processed_agp_file])
            cmd.extend(["--bin_size", str(self.config.bin_size)])
            cmd.extend(["--min_len", str(int(self.config.min_len))])

            if self.config.separate_plots:
                cmd.append("--separate_plots")

            # 设置工作目录
            original_cwd = os.getcwd()
            try:
                os.chdir(plots_dir)
                return self.cmd_runner.run_command(cmd, "生成可视化")
            finally:
                os.chdir(original_cwd)

        except Exception as e:
            self.logger.error(f"💥 可视化生成失败: {e}")
            return False

    def _preprocess_agp_file(self, original_agp_file: str) -> Optional[str]:
        """预处理AGP文件，确保格式正确 | Preprocess AGP file to ensure correct format"""
        try:
            processed_file = original_agp_file.replace('.agp', '.processed.agp')

            with open(original_agp_file, 'r') as infile, open(processed_file, 'w') as outfile:
                line_count = 0
                skipped_count = 0

                for line in infile:
                    line = line.strip()
                    line_count += 1

                    # 跳过空行和注释行
                    if not line or line.startswith('#'):
                        skipped_count += 1
                        continue

                    # 检查列数
                    cols = line.split('\t')
                    if len(cols) < 5:
                        self.logger.warning(f"⚠️ 跳过格式不正确的行 {line_count}: {line}")
                        skipped_count += 1
                        continue

                    # 写入格式正确的行
                    outfile.write(line + '\n')

            self.logger.info(f"📊 AGP文件预处理完成: 共 {line_count} 行，跳过 {skipped_count} 行，输出到 {processed_file}")
            return processed_file

        except Exception as e:
            self.logger.error(f"💥 AGP文件预处理失败: {e}")
            return None