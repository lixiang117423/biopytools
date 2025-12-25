# """
# 🔬 GTX WGS数据处理模块 | GTX WGS Data Processing Module 🔬
# """

# from pathlib import Path
# from .utils import CommandRunner, FileProcessor

# class GTXProcessor:
#     """🔬 GTX数据处理器 | GTX Data Processor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.file_processor = FileProcessor(config, logger)
    
#     def process_sample(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
#         """处理单个样品 🧬 | Process single sample"""
#         self.logger.info(f"🚀 开始处理样品 | Starting to process sample: {sample_name}")
        
#         # 检查输出文件是否已存在 🧐 | Check if output files already exist
#         if self.file_processor.check_output_exists(sample_name):
#             self.logger.info(f"⏭️ 样品 {sample_name} 已处理完成，跳过 | Sample {sample_name} already processed, skipping")
#             return True
        
#         # 定义输出文件路径 📂 | Define output file paths
#         output_vcf = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
#         output_bam = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
#         # 构建Read Group信息 🏷️ | Build Read Group information
#         read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:{sample_name}"
        
#         self.logger.info(f"  📄 R1文件 | R1 file: {Path(r1_file).name}")
#         self.logger.info(f"  📄 R2文件 | R2 file: {Path(r2_file).name}")
#         self.logger.info(f"  📜 输出VCF | Output VCF: vcf/{output_vcf.name}")
#         self.logger.info(f"  💾 输出BAM | Output BAM: bam/{output_bam.name}")
        
#         # 构建GTX WGS命令 💻 | Build GTX WGS command
#         cmd = (
#             f"{self.config.gtx_path} wgs "
#             f"-R \"{read_group}\" "
#             f"-o {output_vcf} "
#             f"-b {output_bam} "
#             f"-t {self.config.threads} "
#             f"--tmp-dir {self.config.tmp_dir} "
#             f"--pcr-indel-model {self.config.pcr_indel_model} "
#             f"--standard-min-confidence-threshold-for-calling {self.config.min_confidence_threshold} "
#             f"--min-base-quality-score {self.config.min_base_quality} "
#             f"--ploidy {self.config.ploidy} "
#             f"{self.config.reference} "
#             f"{r1_file} "
#             f"{r2_file}"
#         )
        
#         # 运行GTX命令 ▶️ | Run GTX command
#         success = self.cmd_runner.run(cmd, f"🧬 GTX WGS分析 (样品: {sample_name}) | GTX WGS analysis (sample: {sample_name})")
        
#         if success:
#             self.logger.info(f"  ✅ ✓ 样品 {sample_name} 处理完成 | Sample {sample_name} processing completed")
            
#             # 检查并报告输出文件大小 📊 | Check and report output file sizes
#             if output_vcf.exists():
#                 vcf_size = self.file_processor.get_file_size(str(output_vcf))
#                 self.logger.info(f"  📏 VCF文件大小 | VCF file size: {vcf_size}")
            
#             if output_bam.exists():
#                 bam_size = self.file_processor.get_file_size(str(output_bam))
#                 self.logger.info(f"  📏 BAM文件大小 | BAM file size: {bam_size}")
            
#             return True
#         else:
#             self.logger.error(f"  ❌ ✗ 样品 {sample_name} 处理失败 | Sample {sample_name} processing failed")
#             return False


"""
🔬 GTX WGS数据处理模块 | GTX WGS Data Processing Module 🔬
"""

from pathlib import Path
from typing import List, Tuple
from .utils import CommandRunner, FileProcessor

class GTXProcessor:
    """🔬 GTX数据处理器 | GTX Data Processor"""
    
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
        
        # 构建GTX WGS命令 💻 | Build GTX WGS command
        cmd = (
            f"faketime '2020-10-20 00:00:00' {self.config.gtx_path} wgs -g "
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
        
        # 运行GTX命令 ▶️ | Run GTX command
        success = self.cmd_runner.run(cmd, f"🧬 GTX WGS分析 (样品: {sample_name}) | GTX WGS analysis (sample: {sample_name})")
        
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


class JointProcessor:
    """🤝 GTX Joint Calling处理器 | GTX Joint Calling Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_processor = FileProcessor(config, logger)
    
    def generate_sample_map(self, processed_samples: List[str]) -> bool:
        """生成sample-name-map文件 📋 | Generate sample-name-map file"""
        self.logger.info("📋 生成样品映射文件 | Generating sample mapping file")
        
        try:
            with open(self.config.sample_map_file, 'w', encoding='utf-8') as f:
                for sample_name in processed_samples:
                    vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
                    if vcf_file.exists():
                        f.write(f"{sample_name}\t{vcf_file}\n")
                        self.logger.debug(f"  添加样品映射 | Added sample mapping: {sample_name} -> {vcf_file}")
                    else:
                        self.logger.warning(f"⚠️ 样品VCF文件不存在 | Sample VCF file does not exist: {vcf_file}")
            
            # 检查生成的文件 ✅ | Check generated file
            if self.config.sample_map_file.exists():
                sample_count = sum(1 for line in open(self.config.sample_map_file) if line.strip())
                self.logger.info(f"✅ 样品映射文件已生成 | Sample mapping file generated: {self.config.sample_map_file}")
                self.logger.info(f"  📊 包含样品数 | Contains samples: {sample_count}")
                
                # 显示前几行内容作为示例
                with open(self.config.sample_map_file, 'r') as f:
                    lines = f.readlines()[:3]
                    self.logger.info("  📝 文件内容示例 | File content example:")
                    for i, line in enumerate(lines, 1):
                        self.logger.info(f"    {i}. {line.strip()}")
                    if len(lines) < sample_count:
                        self.logger.info(f"    ... 还有 {sample_count - len(lines)} 行 | ... and {sample_count - len(lines)} more lines")
                
                return True
            else:
                self.logger.error("❌ 样品映射文件生成失败 | Sample mapping file generation failed")
                return False
                
        except Exception as e:
            self.logger.error(f"❌ 生成样品映射文件时出错 | Error generating sample mapping file: {e}")
            return False
    
    def check_joint_prerequisites(self) -> bool:
        """检查joint calling先决条件 🔍 | Check joint calling prerequisites"""
        self.logger.info("🔍 检查Joint Calling先决条件 | Checking joint calling prerequisites")
        
        # 检查是否有可用的VCF文件
        vcf_files = list(self.config.vcf_output_dir.glob("*.vcf.gz"))
        vcf_count = len(vcf_files)
        
        self.logger.info(f"  📜 可用VCF文件数 | Available VCF files: {vcf_count}")
        
        if vcf_count < 2:
            self.logger.warning("⚠️ Joint calling需要至少2个样品的VCF文件 | Joint calling requires at least 2 sample VCF files")
            return False
        
        # 检查样品映射文件是否存在
        if not self.config.sample_map_file.exists():
            self.logger.warning("⚠️ 样品映射文件不存在 | Sample mapping file does not exist")
            return False
        
        self.logger.info("✅ Joint calling先决条件检查通过 | Joint calling prerequisites check passed")
        return True
    
    def run_joint_calling(self) -> bool:
        """执行joint calling 🤝 | Run joint calling"""
        self.logger.info("🤝 开始GTX Joint Calling | Starting GTX Joint Calling")
        
        # 检查先决条件
        if not self.check_joint_prerequisites():
            self.logger.error("❌ Joint calling先决条件不满足 | Joint calling prerequisites not met")
            return False
        
        # 检查输出文件是否已存在
        if self.config.joint_output_file.exists():
            self.logger.info(f"⏭️ Joint calling结果已存在，跳过 | Joint calling result already exists, skipping: {self.config.joint_output_file}")
            return True
        
        # 构建GTX joint命令
        cmd = (
            f"{self.config.gtx_path} joint "
            f"-r {self.config.reference} "
            f"--sample-name-map {self.config.sample_map_file} "
            f"-o {self.config.joint_output_file}"
            f"--ploidy {self.config.ploidy} "

        )
        
        # 如果配置了joint calling专用的线程数，添加线程参数
        if hasattr(self.config, 'joint_threads') and self.config.joint_threads != self.config.threads:
            cmd += f" -t {self.config.joint_threads}"
        
        self.logger.info(f"  📜 输出文件 | Output file: {self.config.joint_output_file}")
        self.logger.info(f"  📋 样品映射文件 | Sample mapping file: {self.config.sample_map_file}")
        
        # 运行joint calling命令
        success = self.cmd_runner.run(cmd, "🤝 GTX Joint Calling")
        
        if success:
            self.logger.info("✅ GTX Joint Calling完成 | GTX Joint Calling completed")
            
            # 检查并报告输出文件大小
            if self.config.joint_output_file.exists():
                joint_size = self.file_processor.get_file_size(str(self.config.joint_output_file))
                self.logger.info(f"  📏 Joint VCF文件大小 | Joint VCF file size: {joint_size}")
            
            return True
        else:
            self.logger.error("❌ GTX Joint Calling失败 | GTX Joint Calling failed")
            return False
    
    def get_processed_samples(self) -> List[str]:
        """获取已处理的样品列表 📝 | Get list of processed samples"""
        processed_samples = []
        
        # 从VCF输出目录获取已处理的样品
        vcf_files = list(self.config.vcf_output_dir.glob("*.vcf.gz"))
        
        for vcf_file in vcf_files:
            # 从文件名提取样品名（去掉.vcf.gz后缀）
            sample_name = vcf_file.stem.replace('.vcf', '')
            processed_samples.append(sample_name)
        
        processed_samples.sort()  # 排序以确保一致性
        
        self.logger.info(f"📊 检测到已处理样品数 | Detected processed samples: {len(processed_samples)}")
        if processed_samples:
            # 显示前几个样品名作为示例
            sample_examples = processed_samples[:3]
            self.logger.info(f"  样品示例 | Sample examples: {', '.join(sample_examples)}")
            if len(processed_samples) > 3:
                self.logger.info(f"  ... 还有 {len(processed_samples) - 3} 个样品 | ... and {len(processed_samples) - 3} more samples")
        
        return processed_samples