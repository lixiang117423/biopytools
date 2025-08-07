"""
PacBio HiFi结构变异数据处理模块 | PacBio HiFi Structural Variant Data Processing Module
"""

import os
from pathlib import Path
from typing import Dict, List
from .utils import SVCommandRunner


class InputValidator:
    """输入文件验证器 | Input File Validator"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def validate_files(self) -> bool:
        """验证输入文件 | Validate input files"""
        self.runner.logger.info("步骤1: 检查输入文件 | Step 1: Validating input files...")
        
        # 检查BAM索引 | Check BAM index
        bam_index = f"{self.config.bam_file}.bai"
        if not Path(bam_index).exists():
            self.runner.logger.info("创建BAM索引 | Creating BAM index...")
            if not self.runner.run(
                f"{self.config.samtools_path} index {self.config.bam_file}",
                "创建BAM索引"
            ):
                return False
        
        # 检查参考基因组索引 | Check reference genome index
        ref_index = f"{self.config.ref_genome}.fai"
        if not Path(ref_index).exists():
            self.runner.logger.info("创建参考基因组索引 | Creating reference genome index...")
            if not self.runner.run(
                f"{self.config.samtools_path} faidx {self.config.ref_genome}",
                "创建参考基因组索引"
            ):
                return False
        
        # BAM基本统计 | BAM basic statistics
        self.runner.logger.info("生成BAM文件统计 | Generating BAM statistics...")
        stats_file = self.config.logs_dir / f"{self.config.sample_name}_bam_stats.txt"
        flagstat_file = self.config.logs_dir / f"{self.config.sample_name}_flagstat.txt"
        
        self.runner.run(
            f"{self.config.samtools_path} stats {self.config.bam_file} > {stats_file}",
            "BAM统计"
        )
        self.runner.run(
            f"{self.config.samtools_path} flagstat {self.config.bam_file} > {flagstat_file}",
            "BAM flagstat"
        )
        
        self.runner.logger.info("步骤1完成 | Step 1 completed.")
        return True


class SVCaller:
    """结构变异检测器 | Structural Variant Caller"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def run_pbsv(self) -> bool:
        """运行pbsv检测 | Run pbsv detection"""
        self.runner.logger.info("步骤2: 使用pbsv检测结构变异 | Step 2: Running pbsv SV detection...")
        
        # pbsv discover
        svsig_file = self.config.temp_dir / f"{self.config.sample_name}.svsig.gz"
        if not self.runner.run(
            f"{self.config.pbsv_path} discover --hifi {self.config.bam_file} {svsig_file}",
            "pbsv discover"
        ):
            return False
        
        # pbsv call
        vcf_file = self.config.results_dir / f"{self.config.sample_name}_pbsv.vcf"
        if not self.runner.run(
            f"{self.config.pbsv_path} call -j {self.config.threads} -m {self.config.min_sv_length} "
            f"--ccs {self.config.ref_genome} {svsig_file} {vcf_file}",
            "pbsv call"
        ):
            return False
        
        self.runner.logger.info("pbsv检测完成 | pbsv detection completed.")
        return True
    
    def run_sniffles2(self) -> bool:
        """运行Sniffles2检测 | Run Sniffles2 detection"""
        self.runner.logger.info("步骤3: 使用Sniffles2检测结构变异 | Step 3: Running Sniffles2 SV detection...")
        
        vcf_file = self.config.results_dir / f"{self.config.sample_name}_sniffles2.vcf"
        if not self.runner.run(
            f"{self.config.sniffles_path} --input {self.config.bam_file} --vcf {vcf_file} "
            f"--reference {self.config.ref_genome} --threads {self.config.threads} "
            f"--minsvlen {self.config.min_sv_length} --minsupport {self.config.min_support} "
            f"--sample-id {self.config.sample_name} --output-rnames",
            "Sniffles2检测"
        ):
            return False
        
        self.runner.logger.info("Sniffles2检测完成 | Sniffles2 detection completed.")
        return True
    
    def run_cutesv(self) -> bool:
        """运行cuteSV检测 | Run cuteSV detection"""
        self.runner.logger.info("步骤4: 使用cuteSV检测结构变异 | Step 4: Running cuteSV SV detection...")
        
        vcf_file = self.config.results_dir / f"{self.config.sample_name}_cutesv.vcf"
        temp_dir = self.config.temp_dir / "cutesv_temp"
        temp_dir.mkdir(exist_ok=True)
        
        if not self.runner.run(
            f"{self.config.cutesv_path} {self.config.bam_file} {self.config.ref_genome} "
            f"{vcf_file} {temp_dir} --threads {self.config.threads} "
            f"--min_size {self.config.min_sv_length} --min_support {self.config.min_support} "
            f"--genotype --sample {self.config.sample_name}",
            "cuteSV检测"
        ):
            return False
        
        self.runner.logger.info("cuteSV检测完成 | cuteSV detection completed.")
        return True
