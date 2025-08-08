"""
📊 parabricks WGS分析结果处理模块 | parabricks WGS Analysis Results Processing Module 📊
"""

import os
from pathlib import Path

class SummaryGenerator:
    """📋 结果汇总生成器 | Results Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self):
        """生成总结报告 📄 | Generate summary report"""
        self.logger.info("📄 生成处理报告 | Generating processing report")
        
        # 统计输出文件 🔢 | Count output files
        vcf_files = list(self.config.vcf_output_dir.glob("*.vcf.gz"))
        bam_files = list(self.config.bam_output_dir.glob("*.sorted.bam"))
        
        vcf_count = len(vcf_files)
        bam_count = len(bam_files)
        
        self.logger.info(f"📊 生成的文件统计 | Generated files statistics:")
        self.logger.info(f"  📜 VCF文件 | VCF files: {vcf_count} 个")
        self.logger.info(f"  💾 BAM文件 | BAM files: {bam_count} 个")
        
        # 计算目录大小 📏 | Calculate directory sizes
        vcf_dir_size = self._get_directory_size(self.config.vcf_output_dir) if vcf_count > 0 else "0 B"
        bam_dir_size = self._get_directory_size(self.config.bam_output_dir) if bam_count > 0 else "0 B"
        total_size = self._get_directory_size(self.config.output_path)
        
        self.logger.info(f"  📏 VCF目录大小 | VCF directory size: {vcf_dir_size}")
        self.logger.info(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}")
        self.logger.info(f"  📏 总输出大小 | Total output size: {total_size}")
        
        # 生成详细报告文件 ✍️ | Generate detailed report file
        report_file = self.config.output_path / "parabricks_analysis_summary.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("📊 parabricks WGS批处理分析总结报告 | parabricks WGS Batch Processing Analysis Summary Report 📊\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"⚙️ 分析配置 | Analysis Configuration:\n")
            f.write(f"  📂 输入目录 | Input directory: {self.config.input_dir}\n")
            f.write(f"  📤 输出目录 | Output directory: {self.config.output_dir}\n")
            f.write(f"  🧬 参考基因组 | Reference genome: {self.config.reference}\n")
            f.write(f"  🧵 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  💻 parabricks程序路径 | parabricks path: {self.config.parabricks_path}\n\n")
            
            f.write(f"🔬 质控参数 | Quality Control Parameters:\n")
            f.write(f"  🎯 最小置信度阈值 | Min confidence threshold: {self.config.min_confidence_threshold}\n")
            f.write(f"  ✨ 最小碱基质量 | Min base quality: {self.config.min_base_quality}\n")
            f.write(f"  🧬 倍性 | Ploidy: {self.config.ploidy}\n")
            f.write(f"  🔬 PCR indel模型 | PCR indel model: {self.config.pcr_indel_model}\n\n")
            
            f.write(f"📊 文件统计 | File Statistics:\n")
            f.write(f"  📜 VCF文件数量 | VCF file count: {vcf_count}\n")
            f.write(f"  💾 BAM文件数量 | BAM file count: {bam_count}\n")
            f.write(f"  📏 VCF目录大小 | VCF directory size: {vcf_dir_size}\n")
            f.write(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}\n")
            f.write(f"  📏 总输出大小 | Total output size: {total_size}\n\n")
            
            f.write(f"🏗️ 输出目录结构 | Output Directory Structure:\n")
            f.write(f"  {self.config.output_dir}/\n")
            f.write(f"  ├── bam/          # 💾 BAM文件存放目录 | BAM files directory\n")
            f.write(f"  ├── vcf/          # 📜 VCF文件存放目录 | VCF files directory\n")
            f.write(f"  ├── tmp/          # 🗑️ 临时文件目录 | Temporary files directory\n")
            f.write(f"  └── parabricks_analysis_summary.txt  # 📄 分析总结报告 | Analysis summary report\n")
        
        self.logger.info(f"✅ 总结报告已生成 | Summary report generated: {report_file}")
    
    def _get_directory_size(self, path: Path) -> str:
        """获取目录大小 📏 | Get directory size"""
        try:
            total_size = 0
            for dirpath, dirnames, filenames in os.walk(path):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    if os.path.exists(filepath):
                        total_size += os.path.getsize(filepath)
            
            # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
            for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
                if total_size < 1024.0:
                    return f"{total_size:.1f} {unit}"
                total_size /= 1024.0
            return f"{total_size:.1f} PB"
        except:
            return "Unknown"
