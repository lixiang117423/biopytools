"""
HapHiC流程管理模块 - Pipeline模式|HapHiC Pipeline Management Module - Pipeline Mode
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
from ..common.paths import resolve_legacy_path


class HapHiCPipeline:
    """HapHiC流程管理器 - Pipeline模式|HapHiC Pipeline Manager - Pipeline Mode"""

    def __init__(self, config: HapHiCConfig, logger):
        """初始化流程管理器|Initialize pipeline manager"""
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger, config.dry_run)
        self.step_dirs = {}

    def run_complete_pipeline(self) -> bool:
        """运行完整的HapHiC pipeline|Run complete HapHiC pipeline"""
        try:
            self.logger.info("开始HapHiC Pipeline流程|Starting HapHiC Pipeline")
            start_time = time.time()

            # 创建输出目录结构
            self._create_directory_structure()

            # 如果强制重新运行，清理所有目录
            if self.config.force_rerun:
                self.logger.info("强制重新运行模式，清理所有目录|Force rerun mode, cleaning all directories")
                self._force_cleanup_all_directories()

            # 检查HapHiC pipeline是否已完成
            if not self.config.force_rerun and self._check_haphic_pipeline_completed():
                self.logger.info("HapHiC Pipeline已完成，跳过执行|HapHiC Pipeline already completed, skipping execution")
            else:
                # 运行HapHiC pipeline（一步完成）
                if not self._run_haphic_pipeline():
                    return False

            # 生成可视化结果（默认执行）
            if not self.config.force_rerun and self._check_visualization_completed():
                self.logger.info("可视化结果已存在，跳过生成|Visualization results already exist, skipping generation")
            else:
                if self.config.force_rerun:
                    self.logger.info("强制重新运行模式，重新生成可视化|Force rerun mode, regenerating visualization")
                # 如果禁用可视化，直接跳过
                if not self.config.generate_plots:
                    self.logger.info("可视化生成已禁用|Visualization generation disabled")
                else:
                    # 尝试生成可视化，如果失败则警告但继续执行
                    if not self._generate_visualization():
                        self.logger.warning("可视化生成失败，但继续执行其他步骤|Visualization failed, but continuing with other steps")

            # 生成Juicebox文件（如果需要）
            if self.config.output_juicebox:
                if not self.config.force_rerun and self._check_juicebox_files_completed():
                    self.logger.info("Juicebox文件已存在，跳过生成|Juicebox files already exist, skipping generation")
                else:
                    if self.config.force_rerun:
                        self.logger.info("强制重新运行模式，重新生成Juicebox文件|Force rerun mode, regenerating Juicebox files")
                    # 尝试生成Juicebox文件，如果失败则警告但继续执行
                    if not self._generate_juicebox_files():
                        self.logger.warning("Juicebox文件生成失败，但继续执行其他步骤|Juicebox file generation failed, but continuing with other steps")

            # 统计信息
            elapsed_time = time.time() - start_time
            self.logger.info(f"HapHiC Pipeline完成|HapHiC Pipeline completed in {elapsed_time:.2f}秒")

            return True

        except Exception as e:
            self.logger.error(f"HapHiC Pipeline执行失败|HapHiC Pipeline execution failed: {e}")
            return False

    def run_single_step(self, step_name: str) -> bool:
        """运行单个步骤 - 保持兼容性|Run single step - compatibility mode"""
        self.logger.warning(f"单步运行模式已弃用，建议使用完整的pipeline模式|Single-step mode deprecated, recommend using pipeline mode")
        return self.run_complete_pipeline()

    def continue_from_step(self, step_name: str) -> bool:
        """从指定步骤继续 - 保持兼容性|Continue from specified step - compatibility mode"""
        self.logger.warning(f"断点续传模式已弃用，建议使用完整的pipeline模式|Continue-from mode deprecated, recommend using pipeline mode")
        return self.run_complete_pipeline()

    def _create_directory_structure(self):
        """创建目录结构|Create directory structure"""
        self.logger.info("设置目录结构|Setting up directory structure")

        # 创建主输出目录（如果不存在）
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 设置步骤目录路径（不创建，让HapHiC pipeline自己处理）
        self.step_dirs = {
            "cluster": resolve_legacy_path(self.config.output_dir, "01_cluster"),
            "reassign": resolve_legacy_path(self.config.output_dir, "02_reassign"),
            "sort": resolve_legacy_path(self.config.output_dir, "03_sort"),
            "build": resolve_legacy_path(self.config.output_dir, "04_build"),
            "plots": resolve_legacy_path(self.config.output_dir, "05_plots"),
            "juicebox": resolve_legacy_path(self.config.output_dir, "06_juicebox")
        }

        # 只创建juicebox目录，plots目录会在需要时创建
        os.makedirs(self.step_dirs["juicebox"], exist_ok=True)

        # plots目录延迟到实际生成可视化时创建

        self.logger.info("目录结构设置完成|Directory structure setup complete")

    def _run_haphic_pipeline(self) -> bool:
        """运行HapHiC pipeline|Run HapHiC pipeline"""
        self.logger.info("运行HapHiC Pipeline|Running HapHiC Pipeline")

        # 需要重跑时清理所有步骤目录，HapHiC上游使用os.mkdir不支持已存在目录
        # Clean all step directories when re-running, HapHiC upstream uses os.mkdir
        if not self._check_haphic_pipeline_completed():
            step_dirs = ["01_cluster", "02_reassign", "03_sort", "04_build"]
            for step_dir in step_dirs:
                dir_path = os.path.join(self.config.output_dir, step_dir)
                if os.path.exists(dir_path):
                    shutil.rmtree(dir_path)
                    self.logger.info(f"清理步骤目录|Cleaned step directory: {step_dir}")

        self.logger.info("HapHiC将处理已存在的目录|HapHiC will handle existing directories")

        # 构建pipeline命令
        cmd = [self.config.haphic_bin, "pipeline"]
        cmd.extend([self.config.asm_file, self.config.bam_file, str(self.config.nchrs)])

        # 添加所有pipeline参数
        cmd.extend(self._get_pipeline_options())

        # 设置工作目录为输出目录（避免目录冲突）
        original_cwd = os.getcwd()
        try:
            os.chdir(self.config.output_dir)
            self.logger.info(f"工作目录|Working directory: {self.config.output_dir}")

            success = self.cmd_runner.run_command(cmd, "HapHiC Pipeline")

            # Pipeline完成后，检查输出文件
            if success:
                self._verify_pipeline_output()

            return success

        finally:
            os.chdir(original_cwd)

    def _check_haphic_pipeline_completed(self) -> bool:
        """检查HapHiC pipeline是否已完成|Check if HapHiC pipeline is completed"""
        self.logger.info("检查HapHiC Pipeline完成状态|Checking HapHiC Pipeline completion status")

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
                self.logger.info(f"{step_name}步骤已完成|{step_name} step completed")
            else:
                self.logger.info(f"{step_name}步骤未完成|{step_name} step not completed")
                self.logger.info(f"  期待文件|Expected files: {', '.join(files)}")
                all_completed = False

        return all_completed

    def _force_cleanup_all_directories(self):
        """强制清理所有步骤目录|Force clean all step directories"""
        all_dirs = ["01_cluster", "02_reassign", "03_sort", "04_build", "05_plots", "06_juicebox"]

        for dir_name in all_dirs:
            dir_path = os.path.join(self.config.output_dir, dir_name)
            if os.path.exists(dir_path):
                try:
                    shutil.rmtree(dir_path)
                    self.logger.info(f"强制清理目录|Force cleaning directory: {dir_name}")
                except Exception as e:
                    self.logger.warning(f"无法清理目录|Cannot clean directory {dir_name}: {e}")

        # 重新创建必要的目录
        os.makedirs(self.step_dirs["juicebox"], exist_ok=True)

    def _cleanup_existing_directories(self):
        """清理可能存在的HapHiC子目录 - 仅在需要时执行|Clean up existing HapHiC subdirectories - only execute when needed"""
        # 注意：由于实现了断点续传，默认不再清理已存在的目录
        # 只清理非HapHiC pipeline管理的目录（如05_plots）

        # 只清理plots目录，因为它不由HapHiC pipeline管理
        plots_dir = self.step_dirs["plots"]
        if os.path.exists(plots_dir):
            try:
                shutil.rmtree(plots_dir)
                if self.logger:
                    self.logger.info(f"清理plots目录以重新生成可视化|Cleaning plots directory to regenerate visualization")
                # 重新创建空目录
                os.makedirs(plots_dir, exist_ok=True)
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"无法清理plots目录|Cannot clean plots directory: {e}")

        # HapHiC pipeline管理的目录（01_cluster到04_build）不再清理，支持断点续传
        if self.logger:
            self.logger.info("HapHiC pipeline目录保留以支持断点续传|HapHiC pipeline directories preserved for resume support")

    def _get_pipeline_options(self) -> List[str]:
        """获取pipeline命令选项|Get pipeline command options"""
        options = []

        # 输入文件选项|Input file options
        options.extend([
            "--aln_format", self.config.aln_format,
            "--outdir", self.config.output_dir
        ])

        # 基本pipeline控制|Basic pipeline control
        options.extend([
            "--RE", self.config.re_sites,
            "--steps", "1,2,3,4",
            "--threads", str(self.config.threads),
            "--processes", str(self.config.processes)
        ])

        if self.config.quick_view:
            options.append("--quick_view")
        if self.config.normalize_by_nlinks:
            options.append("--normalize_by_nlinks")

        # 组装校正|Assembly correction
        if self.config.correct_nrounds > 0:
            if self.config.correct_resolution is not None:
                correct_resolution = self.config.correct_resolution
                self.logger.info(f"使用用户指定的分辨率|Using user-specified resolution: {correct_resolution}bp")
            else:
                correct_resolution = 1000 if self.config.correct_min_coverage <= 5 else 500
                self.logger.info(f"自动优化分辨率|Auto-optimized resolution: {correct_resolution}bp (基于覆盖率|based on coverage {self.config.correct_min_coverage})")
            options.extend([
                "--correct_nrounds", str(self.config.correct_nrounds),
                "--correct_resolution", str(correct_resolution),
                "--median_cov_ratio", str(self.config.median_cov_ratio),
                "--region_len_ratio", str(self.config.region_len_ratio),
                "--min_region_cutoff", str(self.config.min_region_cutoff)
            ])
            self.logger.info(f"启用组装校正|Assembly correction enabled, resolution: {correct_resolution}bp")

        # 过滤参数|Filtering parameters
        options.extend([
            "--Nx", str(self.config.nx),
            "--RE_site_cutoff", str(self.config.re_site_cutoff),
            "--density_lower", self.config.density_lower,
            "--density_upper", self.config.density_upper,
            "--read_depth_upper", self.config.read_depth_upper,
            "--topN", str(self.config.topN),
            "--rank_sum_hard_cutoff", str(self.config.rank_sum_hard_cutoff),
            "--rank_sum_upper", self.config.rank_sum_upper
        ])

        # 等位基因处理|Allelic link handling
        if self.config.remove_allelic_links:
            options.extend([
                "--remove_allelic_links", str(self.config.remove_allelic_links),
                "--concordance_ratio_cutoff", str(self.config.concordance_ratio_cutoff),
                "--nwindows", str(self.config.nwindows),
                "--max_read_pairs", str(self.config.max_read_pairs),
                "--min_read_pairs", str(self.config.min_read_pairs)
            ])
            if self.config.remove_concentrated_links:
                options.append("--remove_concentrated_links")

        # 分相信息|Phasing information
        if self.config.gfa_files:
            options.extend([
                "--gfa", self.config.gfa_files,
                "--phasing_weight", str(self.config.phasing_weight)
            ])

        # 超长读长|Ultra-long reads
        if self.config.ul_file:
            options.extend([
                "--ul", self.config.ul_file,
                "--min_ul_mapq", str(self.config.min_ul_mapq),
                "--min_ul_alignment_length", str(self.config.min_ul_alignment_length),
                "--max_distance_to_end", str(self.config.max_distance_to_end),
                "--max_overlap_ratio", str(self.config.max_overlap_ratio),
                "--max_gap_len", str(self.config.max_gap_len),
                "--min_ul_support", str(self.config.min_ul_support)
            ])

        # MCL聚类|Markov clustering
        if self.config.skip_clustering:
            options.append("--skip_clustering")

        options.extend([
            "--bin_size", str(self.config.bin_size_kbp),
            "--flank", str(self.config.flank),
            "--min_inflation", str(self.config.min_inflation),
            "--max_inflation", str(self.config.max_inflation),
            "--inflation_step", str(self.config.inflation_step),
            "--expansion", str(self.config.expansion),
            "--max_iter", str(self.config.max_iter),
            "--pruning", str(self.config.pruning)
        ])

        # 重分配|Reassignment
        options.extend([
            "--min_group_len", str(self.config.min_group_len),
            "--max_ctg_len", str(self.config.max_ctg_len),
            "--min_RE_sites", str(self.config.min_re_sites),
            "--min_links", str(self.config.min_links),
            "--min_link_density", str(self.config.min_link_density),
            "--min_density_ratio", str(self.config.min_density_ratio),
            "--ambiguous_cutoff", str(self.config.ambiguous_cutoff),
            "--reassign_nrounds", str(self.config.reassign_nrounds)
        ])

        if self.config.no_additional_rescue:
            options.append("--no_additional_rescue")

        # 排序|Sorting
        if self.config.skip_fast_sort:
            options.append("--skip_fast_sort")

        options.extend([
            "--flanking_region", str(self.config.flanking_region),
            "--density_cal_method", self.config.density_cal_method,
            "--confidence_cutoff", str(self.config.confidence_cutoff)
        ])

        # ALLHiC优化|ALLHiC optimization
        if self.config.skip_allhic:
            options.append("--skip_allhic")

        options.extend([
            "--mutprob", str(self.config.mutprob),
            "--ngen", str(self.config.ngen),
            "--npop", str(self.config.npop),
            "--seed", str(self.config.seed)
        ])

        if self.config.skip_ga:
            options.append("--skipGA")

        # 构建|Building
        options.extend([
            "--Ns", str(self.config.Ns),
            "--max_width", str(self.config.max_width),
            "--prefix", self.config.prefix
        ])

        if self.config.sort_by_input:
            options.append("--sort_by_input")
        if self.config.keep_letter_case:
            options.append("--keep_letter_case")

        # 性能|Performance
        if self.config.dense_matrix:
            options.append("--dense_matrix")

        # 日志|Logging
        if self.config.verbose:
            options.append("--verbose")

        return options

    def _verify_pipeline_output(self):
        """验证pipeline输出|Verify pipeline output"""
        self.logger.info("验证Pipeline输出|Verifying pipeline output")

        output_files = self.config.get_output_files()

        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"{file_type}: {os.path.basename(file_path)} ({size:,} bytes)")
            else:
                self.logger.warning(f"{file_type}不存在|{file_type} not found: {file_path}")

    def _generate_juicebox_files(self) -> bool:
        """生成Juicebox兼容文件 (.hic 和 .assembly)|Generate Juicebox compatible files (.hic and .assembly)"""
        self.logger.info("生成Juicebox文件|Generating Juicebox files")

        try:
            juicebox_dir = self.step_dirs["juicebox"]
            build_dir = self.step_dirs["build"]

            # 查找必要的文件
            output_files = self.config.get_output_files()
            agp_file = output_files["scaffolds_raw_agp"]

            if not FileManager.check_file_exists(agp_file):
                # 尝试使用另一个AGP文件
                agp_file = output_files["scaffolds_agp"]
                if not FileManager.check_file_exists(agp_file):
                    self.logger.error("未找到AGP文件，无法生成assembly文件|AGP file not found, cannot generate assembly file")
                    return False

            # 查找filtered BAM文件 - 可能在多个位置
            possible_bam_locations = [
                os.path.join(self.step_dirs["cluster"], "hic.filtered.bam"),
                self.config.bam_file  # 使用原始BAM文件作为备选
            ]

            hic_filtered_bam = None
            for bam_file in possible_bam_locations:
                if FileManager.check_file_exists(bam_file):
                    hic_filtered_bam = bam_file
                    break

            if not hic_filtered_bam:
                self.logger.warning("未找到hic.filtered.bam文件，使用原始BAM文件|hic.filtered.bam not found, using original BAM file")
                hic_filtered_bam = self.config.bam_file

            # 1. 使用matlock生成.mnd文件
            self.logger.info("步骤1: 生成.mnd文件|Step 1: Generating .mnd file")

            mnd_file = os.path.join(juicebox_dir, "out.links.mnd")

            # 检查mnd文件是否已存在
            if os.path.exists(mnd_file):
                try:
                    size = os.path.getsize(mnd_file)
                    if size > 1000000:  # 至少1MB
                        self.logger.info(f"mnd文件已存在，跳过生成|mnd file already exists: {mnd_file} ({size:,} bytes)")
                    else:
                        self.logger.info(f"mnd文件太小，重新生成|mnd file too small, regenerating: {size:,} bytes")
                        cmd1 = [self.config.matlock_bin, "bam2", "juicer", hic_filtered_bam, mnd_file]
                        if not self.cmd_runner.run_command(cmd1, "生成.mnd文件"):
                            return False
                except Exception as e:
                    self.logger.warning(f"无法检查mnd文件大小|Cannot check mnd file size: {e}")
                    cmd1 = [self.config.matlock_bin, "bam2", "juicer", hic_filtered_bam, mnd_file]
                    if not self.cmd_runner.run_command(cmd1, "生成.mnd文件"):
                        return False
            else:
                cmd1 = [self.config.matlock_bin, "bam2", "juicer", hic_filtered_bam, mnd_file]
                if not self.cmd_runner.run_command(cmd1, "生成.mnd文件"):
                    return False

            # 2. 排序.mnd文件
            self.logger.info("步骤2: 排序.mnd文件|Step 2: Sorting .mnd file")

            sorted_mnd_file = os.path.join(juicebox_dir, "out.sorted.links.mnd")

            # 检查sorted_mnd文件是否已存在
            if os.path.exists(sorted_mnd_file):
                try:
                    size = os.path.getsize(sorted_mnd_file)
                    if size > 1000000:  # 至少1MB
                        self.logger.info(f"排序后的mnd文件已存在，跳过排序|sorted mnd file already exists: {sorted_mnd_file} ({size:,} bytes)")
                        dry_run_success = True
                    else:
                        self.logger.info(f"排序后的mnd文件太小，重新排序|sorted mnd file too small, resorting: {size:,} bytes")
                        cmd2 = ["sort", "-k2,2", "-k6,6", mnd_file]
                        try:
                            with open(sorted_mnd_file, 'w') as outfile:
                                result = subprocess.run(cmd2, stdout=outfile, stderr=subprocess.PIPE, text=True)
                                if result.returncode != 0:
                                    self.logger.error(f"排序.mnd文件失败|Sorting .mnd file failed: {result.stderr}")
                                    return False
                            dry_run_success = True
                        except Exception as e:
                            self.logger.error(f"排序.mnd文件异常|Error sorting .mnd file: {e}")
                            return False
                except Exception as e:
                    self.logger.warning(f"无法检查sorted_mnd文件大小|Cannot check sorted_mnd file size: {e}")
                    cmd2 = ["sort", "-k2,2", "-k6,6", mnd_file]
                    try:
                        with open(sorted_mnd_file, 'w') as outfile:
                            result = subprocess.run(cmd2, stdout=outfile, stderr=subprocess.PIPE, text=True)
                            if result.returncode != 0:
                                self.logger.error(f"排序.mnd文件失败|Sorting .mnd file failed: {result.stderr}")
                                return False
                        dry_run_success = True
                    except Exception as e:
                        self.logger.error(f"排序.mnd文件异常|Error sorting .mnd file: {e}")
                        return False
            else:
                cmd2 = ["sort", "-k2,2", "-k6,6", mnd_file]
                # 使用shell重定向
                if self.config.dry_run:
                    self.logger.info(f"[DRY RUN] {cmd2[0]} {cmd2[1]} > {sorted_mnd_file}")
                    dry_run_success = True
                else:
                    try:
                        with open(sorted_mnd_file, 'w') as outfile:
                            result = subprocess.run(cmd2, stdout=outfile, stderr=subprocess.PIPE, text=True)
                            if result.returncode != 0:
                                self.logger.error(f"排序.mnd文件失败|Sorting .mnd file failed: {result.stderr}")
                                return False
                        dry_run_success = True
                    except Exception as e:
                        self.logger.error(f"排序.mnd文件异常|Error sorting .mnd file: {e}")
                        return False

            if not dry_run_success:
                return False

            # 3. 生成.assembly文件
            self.logger.info("步骤3: 生成.assembly文件|Step 3: Generating .assembly file")

            assembly_file = output_files["assembly_file"]
            cmd3 = ["python3", self.config.agp2assembly_script, agp_file, assembly_file]

            if not self.cmd_runner.run_command(cmd3, "生成.assembly文件"):
                return False

            # 4. 运行asm-visualizer生成.hic文件
            self.logger.info("步骤4: 生成.hic文件|Step 4: Generating .hic file")

            hic_file = output_files["hic_file"]
            cmd4 = ["bash", self.config.asm_visualizer_script, "-p", "false", assembly_file, sorted_mnd_file]

            # 添加调试信息：显示命令和参数
            self.logger.info(f"工作目录|Working directory: {juicebox_dir}")
            self.logger.info(f"Assembly文件|Assembly file: {assembly_file}")
            self.logger.info(f"排序后的mnd文件|Sorted mnd file: {sorted_mnd_file}")
            self.logger.info(f"预期输出.hic文件|Expected .hic file: {hic_file}")
            self.logger.info(f"执行命令|Command: {' '.join(cmd4)}")

            # 设置工作目录为juicebox目录
            original_cwd = os.getcwd()
            try:
                os.chdir(juicebox_dir)
                self.logger.info(f"切换到目录|Changed to directory: {os.getcwd()}")

                if self.config.dry_run:
                    self.logger.info(f"[DRY RUN] cd {juicebox_dir} && {' '.join(cmd4)}")
                    success = True
                else:
                    # 设置环境变量，将Java临时目录指向juicebox目录（避免NFS和权限问题）
                    # 创建一个tmp子目录用于临时文件
                    tmp_dir = os.path.join(juicebox_dir, "tmp_java")
                    os.makedirs(tmp_dir, exist_ok=True)
                    env = os.environ.copy()
                    env["TMPDIR"] = tmp_dir
                    env["TMP"] = tmp_dir
                    env["TEMP"] = tmp_dir
                    env["JAVA_TMPDIR"] = tmp_dir
                    java_opts = env.get("_JAVA_OPTIONS", "")
                    if java_opts:
                        env["_JAVA_OPTIONS"] = f"{java_opts} -Djava.io.tmpdir={tmp_dir}"
                    else:
                        env["_JAVA_OPTIONS"] = f"-Djava.io.tmpdir={tmp_dir}"
                    self.logger.info(f"设置Java临时目录|Set Java temp directory: {tmp_dir}")

                    # 运行asm-visualizer，它会在当前目录生成.hic文件
                    self.logger.info("开始运行run-assembly-visualizer.sh脚本|Starting run-assembly-visualizer.sh script")
                    result = subprocess.run(cmd4, capture_output=True, text=True, cwd=juicebox_dir, env=env)
                    self.logger.info(f"脚本执行完成|Script execution completed with return code: {result.returncode}")

                    # 即使返回码非0，也检查.hic文件是否生成成功
                    # 因为脚本有时会有警告但仍然成功生成文件
                    if result.returncode != 0:
                        self.logger.info(f"检查.hic文件生成状态|Checking .hic file generation status")

                        # 列出当前目录的所有文件，方便调试
                        self.logger.info(f"当前目录内容|Current directory contents:")
                        try:
                            import subprocess as sp
                            ls_result = sp.run(['ls', '-lh', '.'], capture_output=True, text=True)
                            if ls_result.returncode == 0:
                                # 只显示.hic文件
                                for line in ls_result.stdout.split('\n'):
                                    if '.hic' in line or line.startswith('total'):
                                        self.logger.info(f"  {line}")
                        except Exception:
                            pass

                        # 检查可能的.hic文件位置
                        # asm-visualizer脚本生成的文件名是assembly文件名的basename
                        expected_hic_basename = os.path.basename(assembly_file).replace(".assembly", ".hic")
                        expected_hic_relative = expected_hic_basename  # 相对于juicebox_dir
                        expected_hic_absolute = os.path.join(juicebox_dir, expected_hic_basename)

                        self.logger.info(f"检查可能的.hic文件名|Checking possible .hic filenames:")
                        self.logger.info(f"  1. 相对路径(相对名)|Relative filename: {expected_hic_relative}")
                        self.logger.info(f"  2. 绝对路径|Absolute path: {expected_hic_absolute}")
                        self.logger.info(f"  3. assembly.hic|assembly.hic")
                        self.logger.info(f"  4. 最终目标|Target: {hic_file}")

                        # 检查文件是否存在
                        exists_relative = os.path.exists(expected_hic_relative)
                        exists_absolute = os.path.exists(expected_hic_absolute)
                        exists_assembly_hic = os.path.exists("assembly.hic")
                        exists_target = os.path.exists(hic_file)

                        self.logger.info(f"文件存在状态|File existence status:")
                        self.logger.info(f"  {expected_hic_relative}: {exists_relative}")
                        self.logger.info(f"  {expected_hic_absolute}: {exists_absolute}")
                        self.logger.info(f"  assembly.hic: {exists_assembly_hic}")
                        self.logger.info(f"  {hic_file}: {exists_target}")

                        if exists_relative or exists_absolute or exists_assembly_hic or exists_target:
                            self.logger.warning(f"asm-visualizer返回码非0，但.hic文件已生成|asm-visualizer returned non-zero code, but .hic file exists")
                            self.logger.warning(f"返回码|Return code: {result.returncode}")

                            # 输出完整的stdout和stderr，方便调试
                            if result.stdout:
                                self.logger.info(f"标准输出完整内容|Full Stdout:")
                                for line in result.stdout.split('\n')[-50:]:  # 最后50行
                                    self.logger.info(f"  {line}")
                            if result.stderr:
                                self.logger.info(f"标准错误完整内容|Full Stderr:")
                                for line in result.stderr.split('\n')[-50:]:  # 最后50行
                                    self.logger.info(f"  {line}")
                            # 继续处理，不返回False
                        else:
                            self.logger.error(f"生成.hic文件失败|Generating .hic file failed")
                            self.logger.error(f"返回码|Return code: {result.returncode}")
                            if result.stdout:
                                self.logger.error(f"标准输出|Stdout: {result.stdout[:2000]}")
                            if result.stderr:
                                self.logger.error(f"标准错误|Stderr: {result.stderr[:2000]}")
                            return False

                    # asm-visualizer会生成以assembly文件名命名的.hic文件（相对路径）
                    expected_hic_basename = os.path.basename(assembly_file).replace(".assembly", ".hic")
                    self.logger.info(f"查找生成的.hic文件|Looking for generated .hic file: {expected_hic_basename}")

                    if os.path.exists(expected_hic_basename):
                        # 检查是否需要重命名（如果目标文件已存在且不同）
                        if os.path.exists(hic_file) and os.path.abspath(expected_hic_basename) != os.path.abspath(hic_file):
                            self.logger.info(f"目标文件已存在，将覆盖|Target file exists, will overwrite: {hic_file}")
                            os.remove(hic_file)
                        elif os.path.abspath(expected_hic_basename) == os.path.abspath(hic_file):
                            self.logger.info(f".hic文件已在正确位置|.hic file already at correct location: {hic_file}")
                        else:
                            # 重命名文件
                            self.logger.info(f"重命名.hic文件|Renaming .hic file: {expected_hic_basename} -> {hic_file}")
                            os.rename(expected_hic_basename, hic_file)
                            self.logger.info(f"重命名成功|Rename successful")
                    elif os.path.exists("assembly.hic"):
                        # 备选文件名
                        self.logger.info(f"使用备选文件名|Using alternate filename: assembly.hic")
                        if os.path.exists(hic_file):
                            os.remove(hic_file)
                        os.rename("assembly.hic", hic_file)
                        self.logger.info(f"重命名.hic文件（备选名称）|Renamed .hic file (alternate name): assembly.hic -> {hic_file}")
                    else:
                        self.logger.error("未找到生成的.hic文件|Generated .hic file not found")
                        self.logger.error(f"预期文件名|Expected filename: {expected_hic_basename}")
                        return False

                    success = True
                    self.logger.info(f".hic文件处理成功|.hic file processing successful: {hic_file}")

            finally:
                os.chdir(original_cwd)

            if not success:
                return False

            # 5. 验证输出文件
            self.logger.info("验证Juicebox文件|Verifying Juicebox files")

            if not FileManager.check_file_exists(hic_file):
                self.logger.error(f".hic文件未生成|.hic file not generated: {hic_file}")
                return False

            if not FileManager.check_file_exists(assembly_file):
                self.logger.error(f".assembly文件未生成|.assembly file not generated: {assembly_file}")
                return False

            # 获取文件大小信息
            hic_size = FileManager.get_file_size(hic_file) / (1024**2)  # MB
            assembly_size = FileManager.get_file_size(assembly_file) / 1024  # KB

            self.logger.info(f"Juicebox文件生成成功|Juicebox files generated successfully")
            self.logger.info(f" .hic文件大小|.hic file size: {hic_size:.2f} MB")
            self.logger.info(f" .assembly文件大小|.assembly file size: {assembly_size:.2f} KB")
            self.logger.info(f" 文件位置|File location: {juicebox_dir}")

            return True

        except Exception as e:
            self.logger.error(f"Juicebox文件生成失败|Juicebox file generation failed: {e}")
            return False

    def get_step_status(self) -> Dict[str, bool]:
        """获取各步骤状态|Get status of each step - 兼容性方法"""
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
        """检查可视化是否已完成|Check if visualization is completed"""
        plots_dir = self.step_dirs["plots"]

        # 检查是否有PDF或PNG文件
        if not os.path.exists(plots_dir):
            return False

        plot_files = []
        for ext in ['*.pdf', '*.png', '*.jpg', '*.jpeg']:
            plot_files.extend(glob.glob(os.path.join(plots_dir, ext)))

        if plot_files:
            self.logger.info(f"找到 {len(plot_files)} 个可视化文件|Found {len(plot_files)} visualization files")
            return True
        return False

    def _check_juicebox_files_completed(self) -> bool:
        """检查Juicebox文件是否已完成|Check if Juicebox files are completed"""
        juicebox_dir = self.step_dirs["juicebox"]

        # 检查关键文件
        hic_file = os.path.join(juicebox_dir, f"{self.config.prefix}.hic")
        assembly_file = os.path.join(juicebox_dir, f"{self.config.prefix}.assembly")

        if os.path.exists(hic_file) and os.path.exists(assembly_file):
            return True
        return False

    def _generate_visualization(self) -> bool:
        """生成可视化结果|Generate visualization results"""
        self.logger.info("生成Hi-C接触图可视化|Generating Hi-C contact map visualization")

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
                self.logger.info(f"使用标准AGP文件|Using standard AGP file: {agp_file}")
            else:
                self.logger.info(f"使用Raw AGP文件（格式更标准）|Using Raw AGP file (more standard format): {agp_file}")

            if not FileManager.check_file_exists(agp_file):
                self.logger.warning(f"未找到AGP文件|AGP file not found: {agp_file}，跳过可视化|skipping visualization")
                return True

            # 预处理AGP文件，确保格式正确
            processed_agp_file = self._preprocess_agp_file(agp_file)
            if not processed_agp_file:
                self.logger.error("AGP文件预处理失败|AGP file preprocessing failed")
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
            self.logger.error(f"可视化生成失败|Visualization generation failed: {e}")
            return False

    def _preprocess_agp_file(self, original_agp_file: str) -> Optional[str]:
        """预处理AGP文件，确保格式正确|Preprocess AGP file to ensure correct format"""
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
                        self.logger.warning(f"跳过格式不正确的行|Skipping incorrectly formatted line {line_count}: {line}")
                        skipped_count += 1
                        continue

                    # 写入格式正确的行
                    outfile.write(line + '\n')

            self.logger.info(f"AGP文件预处理完成|AGP file preprocessing complete: 共 {line_count} 行|total lines，跳过 {skipped_count} 行|skipped lines，输出到 {processed_file}")
            return processed_file

        except Exception as e:
            self.logger.error(f"AGP文件预处理失败|AGP file preprocessing failed: {e}")
            return None