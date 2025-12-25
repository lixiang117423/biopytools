"""
📊 ENA下载结果汇总模块 | ENA Download Results Summary Module
"""

from pathlib import Path
from datetime import datetime

class ResultsSummary:
    """📋 结果汇总生成器 | Results Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, metadata_file: Path, script_file: Path = None, download_links_count: int = 0):
        """📊 生成汇总报告 | Generate summary report"""
        summary_file = self.config.output_path / 'download_summary.txt'
        
        try:
            with summary_file.open('w', encoding='utf-8') as f:
                f.write("="*60 + "\n")
                f.write("📊 ENA数据下载汇总报告 | ENA Data Download Summary Report\n")
                f.write("="*60 + "\n\n")
                
                # 📋 基本信息 | Basic information
                f.write("📋 基本信息 | Basic Information:\n")
                f.write(f"  - 🎯 项目编号 | Accession: {self.config.accession}\n")
                f.write(f"  - 🌐 下载协议 | Protocol: {self.config.protocol}\n")
                f.write(f"  - ⚙️ 执行方式 | Method: {self.config.method}\n")
                f.write(f"  - 📁 输出目录 | Output Directory: {self.config.output_dir}\n")
                f.write(f"  - 🕐 生成时间 | Generated Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("\n")
                
                # 📊 元数据信息 | Metadata information
                f.write("📊 元数据信息 | Metadata Information:\n")
                if metadata_file and metadata_file.exists():
                    f.write(f"  - 📄 元数据文件 | Metadata File: {metadata_file.name}\n")
                    f.write(f"  - 📏 文件大小 | File Size: {metadata_file.stat().st_size:,} bytes\n")
                    f.write(f"  - 📋 格式 | Format: {self.config.metadata_format.upper()}\n")
                else:
                    f.write("  - ❌ 状态 | Status: 元数据下载失败 | Metadata download failed\n")
                f.write("\n")
                
                # 📥 下载信息 | Download information
                f.write("📥 下载信息 | Download Information:\n")
                f.write(f"  - 📁 发现的FASTQ文件数量 | FASTQ Files Found: {download_links_count}\n")
                
                if script_file and script_file.exists():
                    f.write(f"  - 📜 下载脚本 | Download Script: {script_file.name}\n")
                    f.write(f"  - 📏 脚本大小 | Script Size: {script_file.stat().st_size:,} bytes\n")
                    
                    if self.config.method == "save":
                        f.write("\n")
                        f.write("🚀 下一步操作 | Next Steps:\n")
                        f.write(f"  💡 执行以下命令开始下载 | Run the following command to start download:\n")
                        f.write(f"  bash {script_file.name}\n")
                else:
                    if self.config.method == "save":
                        f.write("  - ❌ 状态 | Status: 下载脚本生成失败 | Download script generation failed\n")
                    else:
                        f.write("  - ✅ 状态 | Status: 直接下载已执行 | Direct download executed\n")
                
                f.write("\n")
                
                # 📁 文件列表 | File list
                f.write("📁 输出文件 | Output Files:\n")
                for file_path in sorted(self.config.output_path.iterdir()):
                    if file_path.is_file():
                        f.write(f"  - 📄 {file_path.name}\n")
                
                f.write("\n")
                f.write("="*60 + "\n")
            
            self.logger.info(f"✅ 汇总报告已生成 | Summary report generated: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 生成汇总报告失败 | Failed to generate summary report: {str(e)}")