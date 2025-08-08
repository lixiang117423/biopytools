"""
🔬 parabricks WGS数据处理模块 | parabricks WGS Data Processing Module 🔬
"""

from pathlib import Path
from .utils import CommandRunner, FileProcessor

class parabricksProcessor:
    """🔬 parabricks数据处理器 | parabricks Data Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_processor = FileProcessor(config, logger)
    
    def process_sample(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
        """处理单个样品 🧬 | Process single sample"""
        self.logger.info(f"🚀 开始处理样品 | Starting to process sample: {sample_name}")
        
        # 检查输出文件是否已存在 🧐 | Check if output files already exist
        if self.file_processor.check_output_exists(sample_name):
            self.logger.info(f"⏭️ 样品 {sample_name} 已处理完成，跳过 | Sample {sample_name} already processed, skipping")
            return True
        
        # 定义输出文件路径 📂 | Define output file paths
        output_vcf = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
        output_bam = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        # 构建Read Group信息 🏷️ | Build Read Group information
        read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:{sample_name}"
        
        self.logger.info(f"  📄 R1文件 | R1 file: {Path(r1_file).name}")
        self.logger.info(f"  📄 R2文件 | R2 file: {Path(r2_file).name}")
        self.logger.info(f"  📜 输出VCF | Output VCF: vcf/{output_vcf.name}")
        self.logger.info(f"  💾 输出BAM | Output BAM: bam/{output_bam.name}")
        
        # 构建parabricks WGS命令 💻 | Build parabricks WGS command
        cmd = (
            f"{self.config.parabricks_path} wgs "
            f"-R \"{read_group}\" "
            f"-o {output_vcf} "
            f"-b {output_bam} "
            f"-t {self.config.threads} "
            f"--tmp-dir {self.config.tmp_dir} "
            f"--pcr-indel-model {self.config.pcr_indel_model} "
            f"--standard-min-confidence-threshold-for-calling {self.config.min_confidence_threshold} "
            f"--min-base-quality-score {self.config.min_base_quality} "
            f"--ploidy {self.config.ploidy} "
            f"{self.config.reference} "
            f"{r1_file} "
            f"{r2_file}"
        )
        
        # 运行parabricks命令 ▶️ | Run parabricks command
        success = self.cmd_runner.run(cmd, f"🧬 parabricks WGS分析 (样品: {sample_name}) | parabricks WGS analysis (sample: {sample_name})")
        
        if success:
            self.logger.info(f"  ✅ ✓ 样品 {sample_name} 处理完成 | Sample {sample_name} processing completed")
            
            # 检查并报告输出文件大小 📊 | Check and report output file sizes
            if output_vcf.exists():
                vcf_size = self.file_processor.get_file_size(str(output_vcf))
                self.logger.info(f"  📏 VCF文件大小 | VCF file size: {vcf_size}")
            
            if output_bam.exists():
                bam_size = self.file_processor.get_file_size(str(output_bam))
                self.logger.info(f"  📏 BAM文件大小 | BAM file size: {bam_size}")
            
            return True
        else:
            self.logger.error(f"  ❌ ✗ 样品 {sample_name} 处理失败 | Sample {sample_name} processing failed")
            return False
