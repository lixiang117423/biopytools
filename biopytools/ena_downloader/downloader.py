"""
📥 ENA FASTQ文件下载模块 | ENA FASTQ File Download Module
"""

import os
from pathlib import Path
from typing import List

class FastqDownloader:
    """📁 FASTQ文件下载器 | FASTQ File Downloader"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_download_command(self, link: str) -> str:
        """🔨 构建下载命令 | Build download command"""
        if self.config.protocol == 'ftp':
            return f'wget -c {link} -P {self.config.output_dir}'
        elif self.config.protocol == 'aspera':
            return f'ascp -v -k 1 -T -l 1000m -P 33001 -i {self.config.aspera_key} era-fasp@{link} {self.config.output_dir}/'
        else:
            raise ValueError(f"❌ Unsupported protocol: {self.config.protocol}")
    
    def generate_download_script(self, download_links: List[str]) -> Path:
        """📜 生成下载脚本 | Generate download script"""
        if not download_links:
            self.logger.warning("⚠️ 没有找到下载链接 | No download links found")
            return None
        
        # 🔨 构建下载命令 | Build download commands
        commands = []
        for link in download_links:
            if link.strip():
                cmd = self.build_download_command(link)
                if cmd:
                    commands.append(cmd)
        
        if not commands:
            self.logger.warning("⚠️ 没有生成有效的下载命令 | No valid download commands generated")
            return None
        
        # 📝 确定脚本文件名 | Determine script filename
        if self.config.protocol == 'aspera':
            script_name = f'download_{self.config.accession}_fastq_by_aspera.sh'
        else:
            script_name = f'download_{self.config.accession}_fastq_by_wget.sh'
        
        script_path = self.config.output_path / script_name
        
        # ✍️ 写入脚本文件 | Write script file
        try:
            with script_path.open('w', newline='\n', encoding='utf-8') as f:
                f.write('#!/bin/bash\n')
                f.write(f'# 📥 ENA FASTQ下载脚本 | ENA FASTQ Download Script\n')
                f.write(f'# 🎯 项目: {self.config.accession}\n')
                f.write(f'# 🌐 协议: {self.config.protocol}\n')
                f.write(f'# 🕐 生成时间: $(date)\n')
                f.write('\n')
                f.write('set -e  # ❌ 遇到错误时退出 | Exit on error\n')
                f.write('\n')
                f.write('echo "🚀 开始下载FASTQ文件..."\n')
                f.write('echo "🚀 Starting FASTQ file download..."\n')
                f.write('\n')
                
                for i, cmd in enumerate(commands, 1):
                    f.write(f'echo "📥 下载文件 {i}/{len(commands)} | Downloading file {i}/{len(commands)}"\n')
                    f.write(f'{cmd}\n')
                    f.write('\n')
                
                f.write('echo "✅ 所有文件下载完成！"\n')
                f.write('echo "✅ All files downloaded successfully!"\n')
            
            # 🔐 设置执行权限 | Set executable permission
            script_path.chmod(0o755)
            
            self.logger.info(f"✅ 下载脚本已生成 | Download script generated: {script_path}")
            return script_path
            
        except IOError as e:
            self.logger.error(f"❌ 生成下载脚本失败 | Failed to generate download script: {str(e)}")
            return None
    
    def execute_downloads(self, download_links: List[str]) -> bool:
        """🚀 直接执行下载 | Execute downloads directly"""
        if not download_links:
            self.logger.warning("⚠️ 没有找到下载链接 | No download links found")
            return False
        
        self.logger.info(f"🚀 开始直接下载 {len(download_links)} 个文件 | Starting direct download of {len(download_links)} files")
        
        success_count = 0
        for i, link in enumerate(download_links, 1):
            if not link.strip():
                continue
                
            self.logger.info(f"📥 下载文件 {i}/{len(download_links)} | Downloading file {i}/{len(download_links)}")
            cmd = self.build_download_command(link)
            
            if cmd and self.cmd_runner.run_command(cmd):
                success_count += 1
            else:
                self.logger.error(f"❌ 下载失败 | Download failed: {link}")
        
        self.logger.info(f"📊 下载完成: {success_count}/{len(download_links)} 个文件成功 | Download completed: {success_count}/{len(download_links)} files successful")
        return success_count > 0