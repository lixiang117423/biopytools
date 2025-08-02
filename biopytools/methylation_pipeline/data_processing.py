"""
数据处理模块 | Data Processing Module
包含文件重命名、质控、比对等基础数据处理功能
"""

import os
import shutil
from pathlib import Path

from .utils import CommandRunner


class FileManager:
    """文件管理器 | File Manager"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def rename_files(self) -> bool:
        """重命名和清理原始文件 | Rename and clean raw files"""
        self.logger.info("步骤1: 开始重命名和清理原始文件...")

        try:
            # 重命名包含URL参数的文件
            for filename in os.listdir(self.config.raw_dir):
                if "?" in filename:
                    old_path = os.path.join(self.config.raw_dir, filename)
                    new_name = filename.split("?")[0]
                    new_path = os.path.join(self.config.raw_dir, new_name)

                    self.logger.info(f"重命名: {filename} -> {new_name}")
                    shutil.move(old_path, new_path)

            # 清理文件名前缀
            for filename in os.listdir(self.config.raw_dir):
                if filename.startswith("FZYM412_"):
                    old_path = os.path.join(self.config.raw_dir, filename)
                    new_name = filename.replace("FZYM412_", "")
                    new_path = os.path.join(self.config.raw_dir, new_name)

                    self.logger.info(f"去前缀: {filename} -> {new_name}")
                    shutil.move(old_path, new_path)

            self.logger.info("步骤1: 文件重命名和清理完成")
            return True

        except Exception as e:
            self.logger.error(f"文件重命名失败: {e}")
            return False


class QualityController:
    """质量控制器 | Quality Controller"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_fastp(self) -> bool:
        """使用fastp进行质控 | Run fastp for quality control"""
        self.logger.info("步骤2: 开始使用fastp进行质控...")

        # 创建输出目录
        os.makedirs(self.config.clean_dir, exist_ok=True)

        # 获取所有样品名称
        samples = []
        for filename in os.listdir(self.config.raw_dir):
            if filename.endswith(".fq.gz") and (
                "_1.fq.gz" in filename or "_2.fq.gz" in filename
            ):
                sample = filename.replace("_1.fq.gz", "").replace("_2.fq.gz", "")
                if sample not in samples:
                    samples.append(sample)

        if not samples:
            self.logger.error(f"在 {self.config.raw_dir} 中没有找到 *.fq.gz 文件")
            return False

        self.logger.info(f"发现样品: {samples}")

        # 处理每个样品
        for sample in samples:
            clean_sample = sample.replace("FZYM412_", "")
            group = self.config.sample_groups.get(clean_sample, "Unknown")
            self.logger.info(
                f"处理样品: {sample} (清理后: {clean_sample}, 分组: {group})"
            )

            r1 = os.path.join(self.config.raw_dir, f"{sample}_1.fq.gz")
            r2 = os.path.join(self.config.raw_dir, f"{sample}_2.fq.gz")
            clean_r1 = os.path.join(
                self.config.clean_dir, f"{clean_sample}_1_clean.fq.gz"
            )
            clean_r2 = os.path.join(
                self.config.clean_dir, f"{clean_sample}_2_clean.fq.gz"
            )

            # 检查是否已经存在clean文件
            if os.path.exists(clean_r1) and os.path.exists(clean_r2):
                self.logger.info(f"  ✅ 跳过 {clean_sample} - clean文件已存在")
                continue

            if os.path.exists(r1) and os.path.exists(r2):
                cmd = (
                    f"{self.config.fastp_path} "
                    f"-i {r1} -I {r2} "
                    f"-o {clean_r1} -O {clean_r2} "
                    f"-h {os.path.join(self.config.clean_dir, f'{clean_sample}_fastp.html')} "
                    f"-j {os.path.join(self.config.clean_dir, f'{clean_sample}_fastp.json')} "
                    f"--thread {self.config.threads} "
                    f"--detect_adapter_for_pe --correction --cut_front --cut_tail "
                    f"--cut_mean_quality 20 --qualified_quality_phred 20 "
                    f"--unqualified_percent_limit 40 --n_base_limit 10 --length_required 36"
                )

                if not self.cmd_runner.run(cmd, f"fastp处理: {clean_sample}"):
                    return False

                self.logger.info(f"  fastp完成: {clean_sample}")
            else:
                self.logger.warning(f"  警告: 找不到配对文件 - {sample}")

        self.logger.info("步骤2: fastp质控完成")
        return True


class IndexBuilder:
    """索引构建器 | Index Builder"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_bismark_index(self) -> bool:
        """构建Bismark索引 | Build Bismark index"""
        self.logger.info("步骤3: 检查并构建Bismark索引...")

        bismark_index_dir = os.path.join(self.config.mapping_dir, "bismark_index")
        os.makedirs(bismark_index_dir, exist_ok=True)

        # 复制基因组文件到索引目录
        genome_dest = os.path.join(
            bismark_index_dir, os.path.basename(self.config.genome_fa)
        )
        if not os.path.exists(genome_dest):
            shutil.copy2(self.config.genome_fa, genome_dest)

        # 获取bowtie2路径
        bowtie2_path = shutil.which(self.config.bowtie2_path)
        if not bowtie2_path:
            self.logger.error("找不到bowtie2，请检查安装")
            return False

        bowtie2_dir = os.path.dirname(bowtie2_path)
        self.logger.info(f"使用bowtie2路径: {bowtie2_dir}")

        # 构建索引
        cmd = (
            f"{self.config.bismark_genome_preparation_path} "
            f"--path_to_aligner {bowtie2_dir} "
            f"--verbose {bismark_index_dir}"
        )

        if self.cmd_runner.run(cmd, "构建Bismark索引"):
            self.logger.info("Bismark索引构建完成")
            return True
        else:
            # 尝试让bismark自动查找
            self.logger.info("尝试让bismark自动查找bowtie2...")
            cmd_auto = (
                f"{self.config.bismark_genome_preparation_path} "
                f"--verbose {bismark_index_dir}"
            )

            if self.cmd_runner.run(cmd_auto, "构建Bismark索引（自动查找）"):
                self.logger.info("Bismark索引构建完成")
                return True
            else:
                self.logger.error("Bismark索引构建失败")
                return False
