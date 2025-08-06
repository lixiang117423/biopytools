"""
Bismark流程核心分析模块 | Bismark Pipeline Core Analysis Module
"""
# (此文件无改动 | No changes in this file)
import os
import glob
import shutil
from pathlib import Path
from .utils import CommandRunner
from .results import ResultsManager

class CorePipeline:
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.results_manager = ResultsManager(config, logger)

    def build_bismark_index(self) -> bool:
        self.logger.info(f"🏗️ 步骤 1: 检查并构建Bismark索引于基因组目录 | Step 1: Check and build Bismark index in genome directory: {self.config.genome_dir}")
        
        index_complete_file = Path(self.config.genome_dir) / "Bisulfite_Genome" / "CT_conversion" / "genome_mfa.CT_conversion.fa"
        if index_complete_file.exists():
            self.logger.info("✅ Bismark索引已存在，跳过构建 | Bismark index already exists, skipping build.")
            return True

        cmd_runner = CommandRunner(self.logger, Path(self.config.genome_dir))

        bowtie2_path = shutil.which(self.config.bowtie2_path)
        if not bowtie2_path:
            self.logger.error("❌ 找不到bowtie2，请检查安装 | bowtie2 not found, please check installation.")
            return False
        bowtie2_dir = os.path.dirname(bowtie2_path)
        self.logger.info(f"  -> 使用bowtie2路径 | Using bowtie2 path: {bowtie2_dir}")

        cmd = (
            f"\"{self.config.bismark_genome_preparation_path}\" "
            f"--path_to_aligner \"{bowtie2_dir}\" "
            f"--verbose "
            f"\"{self.config.genome_dir}\""
        )
        
        if cmd_runner.run(cmd, "构建Bismark索引 | Building Bismark index"):
            self.logger.info("✅ Bismark索引构建完成 | Bismark index build complete.")
            return True
        else:
            self.logger.warning("  ⚠️ 主索引构建命令失败，尝试备用方案 | Main index build command failed, trying fallback.")
            cmd_auto = (
                f"\"{self.config.bismark_genome_preparation_path}\" "
                f"--verbose \"{self.config.genome_dir}\""
            )
            if cmd_runner.run(cmd_auto, "构建Bismark索引 (自动查找bowtie2) | Building Bismark index (auto-detect bowtie2)"):
                self.logger.info("✅ Bismark索引构建完成 (备用方案成功) | Bismark index build complete (fallback succeeded).")
                return True
            else:
                self.logger.error("❌ Bismark索引构建失败 (主方案和备用方案均失败) | Bismark index build failed (main and fallback attempts failed).")
                return False

    def _get_r2_pattern(self, r1_pattern: str) -> str:
        if '1' in r1_pattern:
            return r1_pattern.replace('1', '2', 1)
        if 'R1' in r1_pattern:
            return r1_pattern.replace('R1', 'R2', 1)
        raise ValueError(
            f"无法从R1模式 '{r1_pattern}' 推断出R2模式。模式中应包含 '1' 或 'R1'。"
            f"Cannot infer R2 pattern from R1 pattern '{r1_pattern}'. It should contain '1' or 'R1'."
        )

    def run_mapping_and_extraction(self) -> bool:
        self.logger.info("🧬 步骤 2: 开始甲基化比对和提取 | Step 2: Starting methylation mapping and extraction...")
        
        r1_pattern = self.config.pattern
        try: r2_pattern = self._get_r2_pattern(r1_pattern)
        except ValueError as e: self.logger.error(e); return False

        self.logger.info(f"  使用文件模式 | Using file patterns: R1: *{r1_pattern}, R2: *{r2_pattern}")
        
        samples = {}
        for f in os.listdir(self.config.raw_dir):
            if f.endswith(r1_pattern):
                sample_name = f[:-len(r1_pattern)]
                r1_path = os.path.join(self.config.raw_dir, f)
                r2_path = os.path.join(self.config.raw_dir, f"{sample_name}{r2_pattern}")
                if os.path.exists(r2_path):
                    samples[sample_name] = {'R1': r1_path, 'R2': r2_path}
                else:
                    self.logger.warning(f"  ⚠️ 发现R1文件但未找到对应的R2文件 | Found R1 but missing R2: {f} (期望 | expected: {os.path.basename(r2_path)})")
        
        if not samples:
            self.logger.error(f"❌ 未根据模式 '*{r1_pattern}' 找到配对的 FASTQ 文件 | No paired FASTQ files found based on pattern '*{r1_pattern}'.")
            return False

        self.logger.info(f"🎯 发现 {len(samples)} 个样品 | Found {len(samples)} samples: {', '.join(samples.keys())}")
        
        for sample, files in samples.items():
            self.logger.info(f"--- 开始处理样品 | Starting to process sample: {sample} ---")
            
            sample_result_dir = Path(self.config.result_dir) / sample
            if list(sample_result_dir.glob("*CX_report.txt")):
                self.logger.info(f"🚀 跳过样品 {sample} - 最终结果文件已存在 | Skipping sample {sample} - final result files already exist.")
                continue

            sample_tmp_dir = Path(self.config.tmp_dir) / sample
            sample_tmp_dir.mkdir(parents=True, exist_ok=True)
            cmd_runner = CommandRunner(self.logger, sample_tmp_dir)
            
            align_cmd = (
                f"\"{self.config.bismark_path}\" --genome \"{self.config.genome_dir}\" "
                f"--multicore {max(1, self.config.threads // 4)} "
                f"-1 \"{files['R1']}\" -2 \"{files['R2']}\" "
                f"--output_dir \"{sample_tmp_dir}\" "
                f"--temp_dir \"{sample_tmp_dir}\" "
                f"--bowtie2 --non_directional"
            )

            if not cmd_runner.run(align_cmd, f"Bismark比对 | Bismark alignment: {sample}"):
                self.logger.error(f"❌ 样品 {sample} 的比对步骤失败，跳过 | Alignment step failed for sample {sample}, skipping.")
                continue
            
            temp_bam_files = list(sample_tmp_dir.glob("*_pe.bam"))
            if not temp_bam_files:
                self.logger.error(f"❌ 未在临时目录中找到BAM文件 | BAM file not found in temp directory: {sample}")
                continue
            temp_bam_path = temp_bam_files[0]

            extract_cmd_parts = [
                f"\"{self.config.bismark_methylation_extractor_path}\"",
                "--paired-end",
                "--report", "--bedGraph", "--cytosine_report", "--CX_context",
                f"--genome_folder \"{self.config.genome_dir}\"",
                f"--buffer_size {self.config.sort_buffer}",
                f"--multicore {self.config.threads}",
                f"\"{temp_bam_path}\""
            ]
            if self.config.no_overlap: extract_cmd_parts.insert(1, "--no_overlap")
            
            if not cmd_runner.run(" ".join(extract_cmd_parts), f"甲基化提取 | Methylation extraction: {sample}"):
                self.logger.error(f"❌ 样品 {sample} 的甲基化提取步骤失败，跳过 | Methylation extraction step failed for sample {sample}, skipping.")
                continue

            self._organize_and_split_results(sample, temp_bam_path, sample_tmp_dir, sample_result_dir)
            self.logger.info(f"✅ 样品 {sample} 处理完成 | Sample {sample} processing complete.")
        
        self.logger.info("✅ 所有样品处理流程已完成 | All sample processing pipelines are complete.")
        return True

    def _organize_and_split_results(self, sample: str, temp_bam_path: Path, sample_tmp_dir: Path, sample_result_dir: Path):
        self.logger.info(f"  整理并拆分最终文件 | Organizing and splitting final files...")
        
        final_bam_dest = Path(self.config.mapping_dir) / temp_bam_path.name
        self.logger.info(f"    -> 移动BAM文件到 | Moving BAM file to: {final_bam_dest}")
        shutil.move(str(temp_bam_path), str(final_bam_dest))
        
        sample_result_dir.mkdir(parents=True, exist_ok=True)
        result_patterns = ["*CX_report.txt", "*.bedGraph", "*_splitting_report.txt", "*M-bias.txt"]
        final_cx_report_path = None
        for pattern in result_patterns:
            for f in sample_tmp_dir.glob(pattern):
                dest_path = sample_result_dir / f.name
                shutil.move(str(f), str(dest_path))
                self.logger.info(f"    -> 移动结果文件到 | Moving result file to: {dest_path}")
                if dest_path.name.endswith("CX_report.txt"):
                    final_cx_report_path = dest_path

        if final_cx_report_path:
            self.logger.info(f"  拆分CX报告文件 | Splitting CX_report file: {final_cx_report_path.name}")
            self.results_manager.split_cx_report(final_cx_report_path)
        else:
            self.logger.warning(f"  ⚠️ 未找到CX报告文件进行拆分 | CX_report file not found for splitting.")
