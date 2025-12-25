# """
# 📊 GTX WGS分析结果处理模块 | GTX WGS Analysis Results Processing Module 📊
# """

# import os
# from pathlib import Path

# class SummaryGenerator:
#     """📋 结果汇总生成器 | Results Summary Generator"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_summary_report(self):
#         """生成总结报告 📄 | Generate summary report"""
#         self.logger.info("📄 生成处理报告 | Generating processing report")
        
#         # 统计输出文件 🔢 | Count output files
#         vcf_files = list(self.config.vcf_output_dir.glob("*.vcf.gz"))
#         bam_files = list(self.config.bam_output_dir.glob("*.sorted.bam"))
        
#         vcf_count = len(vcf_files)
#         bam_count = len(bam_files)
        
#         self.logger.info(f"📊 生成的文件统计 | Generated files statistics:")
#         self.logger.info(f"  📜 VCF文件 | VCF files: {vcf_count} 个")
#         self.logger.info(f"  💾 BAM文件 | BAM files: {bam_count} 个")
        
#         # 计算目录大小 📏 | Calculate directory sizes
#         vcf_dir_size = self._get_directory_size(self.config.vcf_output_dir) if vcf_count > 0 else "0 B"
#         bam_dir_size = self._get_directory_size(self.config.bam_output_dir) if bam_count > 0 else "0 B"
#         total_size = self._get_directory_size(self.config.output_path)
        
#         self.logger.info(f"  📏 VCF目录大小 | VCF directory size: {vcf_dir_size}")
#         self.logger.info(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}")
#         self.logger.info(f"  📏 总输出大小 | Total output size: {total_size}")
        
#         # 生成详细报告文件 ✍️ | Generate detailed report file
#         report_file = self.config.output_path / "gtx_analysis_summary.txt"
#         with open(report_file, 'w', encoding='utf-8') as f:
#             f.write("📊 GTX WGS批处理分析总结报告 | GTX WGS Batch Processing Analysis Summary Report 📊\n")
#             f.write("=" * 80 + "\n\n")
            
#             f.write(f"⚙️ 分析配置 | Analysis Configuration:\n")
#             f.write(f"  📂 输入目录 | Input directory: {self.config.input_dir}\n")
#             f.write(f"  📤 输出目录 | Output directory: {self.config.output_dir}\n")
#             f.write(f"  🧬 参考基因组 | Reference genome: {self.config.reference}\n")
#             f.write(f"  🧵 线程数 | Threads: {self.config.threads}\n")
#             f.write(f"  💻 GTX程序路径 | GTX path: {self.config.gtx_path}\n\n")
            
#             f.write(f"🔬 质控参数 | Quality Control Parameters:\n")
#             f.write(f"  🎯 最小置信度阈值 | Min confidence threshold: {self.config.min_confidence_threshold}\n")
#             f.write(f"  ✨ 最小碱基质量 | Min base quality: {self.config.min_base_quality}\n")
#             f.write(f"  🧬 倍性 | Ploidy: {self.config.ploidy}\n")
#             f.write(f"  🔬 PCR indel模型 | PCR indel model: {self.config.pcr_indel_model}\n\n")
            
#             f.write(f"📊 文件统计 | File Statistics:\n")
#             f.write(f"  📜 VCF文件数量 | VCF file count: {vcf_count}\n")
#             f.write(f"  💾 BAM文件数量 | BAM file count: {bam_count}\n")
#             f.write(f"  📏 VCF目录大小 | VCF directory size: {vcf_dir_size}\n")
#             f.write(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}\n")
#             f.write(f"  📏 总输出大小 | Total output size: {total_size}\n\n")
            
#             f.write(f"🏗️ 输出目录结构 | Output Directory Structure:\n")
#             f.write(f"  {self.config.output_dir}/\n")
#             f.write(f"  ├── bam/          # 💾 BAM文件存放目录 | BAM files directory\n")
#             f.write(f"  ├── vcf/          # 📜 VCF文件存放目录 | VCF files directory\n")
#             f.write(f"  ├── tmp/          # 🗑️ 临时文件目录 | Temporary files directory\n")
#             f.write(f"  └── gtx_analysis_summary.txt  # 📄 分析总结报告 | Analysis summary report\n")
        
#         self.logger.info(f"✅ 总结报告已生成 | Summary report generated: {report_file}")
    
#     def _get_directory_size(self, path: Path) -> str:
#         """获取目录大小 📏 | Get directory size"""
#         try:
#             total_size = 0
#             for dirpath, dirnames, filenames in os.walk(path):
#                 for filename in filenames:
#                     filepath = os.path.join(dirpath, filename)
#                     if os.path.exists(filepath):
#                         total_size += os.path.getsize(filepath)
            
#             # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
#             for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
#                 if total_size < 1024.0:
#                     return f"{total_size:.1f} {unit}"
#                 total_size /= 1024.0
#             return f"{total_size:.1f} PB"
#         except:
#             return "Unknown"


"""
📊 GTX WGS分析结果处理模块 | GTX WGS Analysis Results Processing Module 📊
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
        
        # 统计Joint calling结果
        joint_files = []
        joint_vcf_exists = False
        sample_map_exists = False
        
        if self.config.enable_joint_calling:
            joint_files = list(self.config.joint_output_dir.glob("*"))
            joint_vcf_exists = self.config.joint_output_file.exists()
            sample_map_exists = self.config.sample_map_file.exists()
        
        self.logger.info(f"📊 生成的文件统计 | Generated files statistics:")
        self.logger.info(f"  📜 单样品VCF文件 | Individual VCF files: {vcf_count} 个")
        self.logger.info(f"  💾 BAM文件 | BAM files: {bam_count} 个")
        
        if self.config.enable_joint_calling:
            self.logger.info(f"  🤝 Joint calling文件 | Joint calling files: {len(joint_files)} 个")
            self.logger.info(f"    ✓ Joint VCF文件: {'存在' if joint_vcf_exists else '不存在'}")
            self.logger.info(f"    ✓ 样品映射文件: {'存在' if sample_map_exists else '不存在'}")
        
        # 计算目录大小 📏 | Calculate directory sizes
        vcf_dir_size = self._get_directory_size(self.config.vcf_output_dir) if vcf_count > 0 else "0 B"
        bam_dir_size = self._get_directory_size(self.config.bam_output_dir) if bam_count > 0 else "0 B"
        joint_dir_size = "0 B"
        
        if self.config.enable_joint_calling and joint_files:
            joint_dir_size = self._get_directory_size(self.config.joint_output_dir)
        
        total_size = self._get_directory_size(self.config.output_path)
        
        self.logger.info(f"  📏 单样品VCF目录大小 | Individual VCF directory size: {vcf_dir_size}")
        self.logger.info(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}")
        
        if self.config.enable_joint_calling:
            self.logger.info(f"  📏 Joint calling目录大小 | Joint calling directory size: {joint_dir_size}")
        
        self.logger.info(f"  📏 总输出大小 | Total output size: {total_size}")
        
        # 生成详细报告文件 ✍️ | Generate detailed report file
        report_file = self.config.output_path / "gtx_analysis_summary.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("📊 GTX WGS批处理分析总结报告 | GTX WGS Batch Processing Analysis Summary Report 📊\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"⚙️ 分析配置 | Analysis Configuration:\n")
            f.write(f"  📂 输入目录 | Input directory: {self.config.input_dir}\n")
            f.write(f"  📤 输出目录 | Output directory: {self.config.output_dir}\n")
            f.write(f"  🧬 参考基因组 | Reference genome: {self.config.reference}\n")
            f.write(f"  🧵 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  💻 GTX程序路径 | GTX path: {self.config.gtx_path}\n")
            f.write(f"  🤝 Joint Calling: {'启用' if self.config.enable_joint_calling else '未启用'} | {'Enabled' if self.config.enable_joint_calling else 'Disabled'}\n")
            
            if self.config.enable_joint_calling:
                f.write(f"  📜 Joint输出文件 | Joint output file: {self.config.joint_output_name}\n")
                f.write(f"  🧵 Joint线程数 | Joint threads: {self.config.joint_threads}\n")
            
            f.write("\n")
            
            f.write(f"🔬 质控参数 | Quality Control Parameters:\n")
            f.write(f"  🎯 最小置信度阈值 | Min confidence threshold: {self.config.min_confidence_threshold}\n")
            f.write(f"  ✨ 最小碱基质量 | Min base quality: {self.config.min_base_quality}\n")
            f.write(f"  🧬 倍性 | Ploidy: {self.config.ploidy}\n")
            f.write(f"  🔬 PCR indel模型 | PCR indel model: {self.config.pcr_indel_model}\n\n")
            
            f.write(f"📊 文件统计 | File Statistics:\n")
            f.write(f"  📜 单样品VCF文件数量 | Individual VCF file count: {vcf_count}\n")
            f.write(f"  💾 BAM文件数量 | BAM file count: {bam_count}\n")
            
            if self.config.enable_joint_calling:
                f.write(f"  🤝 Joint calling文件数量 | Joint calling file count: {len(joint_files)}\n")
                f.write(f"    - Joint VCF: {'✓' if joint_vcf_exists else '✗'}\n")
                f.write(f"    - 样品映射文件 | Sample mapping: {'✓' if sample_map_exists else '✗'}\n")
            
            f.write(f"  📏 单样品VCF目录大小 | Individual VCF directory size: {vcf_dir_size}\n")
            f.write(f"  📏 BAM目录大小 | BAM directory size: {bam_dir_size}\n")
            
            if self.config.enable_joint_calling:
                f.write(f"  📏 Joint calling目录大小 | Joint calling directory size: {joint_dir_size}\n")
            
            f.write(f"  📏 总输出大小 | Total output size: {total_size}\n\n")
            
            # 输出目录结构
            f.write(f"🏗️ 输出目录结构 | Output Directory Structure:\n")
            f.write(f"  {self.config.output_dir}/\n")
            f.write(f"  ├── bam/          # 💾 BAM文件存放目录 | BAM files directory ({bam_count} 个文件)\n")
            f.write(f"  ├── vcf/          # 📜 单样品VCF文件存放目录 | Individual VCF files directory ({vcf_count} 个文件)\n")
            
            if self.config.enable_joint_calling:
                f.write(f"  ├── joint/        # 🤝 Joint calling结果目录 | Joint calling results directory\n")
                f.write(f"  │   ├── {self.config.joint_output_name}  # 📜 Joint VCF文件 | Joint VCF file\n")
                f.write(f"  │   └── sample_map.txt     # 📋 样品映射文件 | Sample mapping file\n")
            
            f.write(f"  ├── tmp/          # 🗑️ 临时文件目录 | Temporary files directory\n")
            f.write(f"  └── gtx_analysis_summary.txt  # 📄 分析总结报告 | Analysis summary report\n\n")
            
            # 添加样品详细信息
            if vcf_files:
                f.write(f"📝 处理样品详细信息 | Processed Samples Details:\n")
                for i, vcf_file in enumerate(vcf_files, 1):
                    sample_name = vcf_file.stem.replace('.vcf', '')
                    vcf_size = self._get_file_size(str(vcf_file))
                    
                    # 查找对应的BAM文件
                    bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
                    bam_size = self._get_file_size(str(bam_file)) if bam_file.exists() else "N/A"
                    
                    f.write(f"  {i:2d}. {sample_name}\n")
                    f.write(f"      VCF: {vcf_size:>10s}  BAM: {bam_size:>10s}\n")
            
            # Joint calling详细信息
            if self.config.enable_joint_calling and joint_vcf_exists:
                f.write(f"\n🤝 Joint Calling结果详情 | Joint Calling Results Details:\n")
                joint_vcf_size = self._get_file_size(str(self.config.joint_output_file))
                f.write(f"  📜 Joint VCF文件 | Joint VCF file: {self.config.joint_output_file.name}\n")
                f.write(f"  📏 文件大小 | File size: {joint_vcf_size}\n")
                
                # 读取样品映射文件信息
                if sample_map_exists:
                    sample_count = sum(1 for line in open(self.config.sample_map_file) if line.strip())
                    f.write(f"  📋 包含样品数 | Included samples: {sample_count}\n")
                    f.write(f"  📄 样品映射文件 | Sample mapping file: {self.config.sample_map_file.name}\n")
            
            f.write(f"\n🎯 分析完成状态 | Analysis Completion Status:\n")
            f.write(f"  ✅ 单样品分析: 完成 {vcf_count} 个样品 | Individual analysis: {vcf_count} samples completed\n")
            
            if self.config.enable_joint_calling:
                joint_status = "完成" if joint_vcf_exists else "失败"
                f.write(f"  🤝 Joint calling: {joint_status} | {joint_status}\n")
        
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
            
            return self._format_size(total_size)
        except:
            return "Unknown"
    
    def _get_file_size(self, file_path: str) -> str:
        """获取单个文件大小 📏 | Get single file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            return self._format_size(size_bytes)
        except:
            return "Unknown"
    
    def _format_size(self, size_bytes: int) -> str:
        """格式化文件大小 📏 | Format file size"""
        # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f} PB"