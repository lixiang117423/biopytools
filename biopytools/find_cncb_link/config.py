"""
CNCB FTP配置管理模块 | CNCB FTP Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

@dataclass
class Config:
    """CNCB FTP配置类 | CNCB FTP Configuration Class"""
    
    # FTP服务器配置 | FTP Server Configuration
    ftp_host: str = "download2.cncb.ac.cn"
    ftp_timeout: int = 60
    
    # 基础目录配置 | Base Directory Configuration
    base_dirs: List[str] = None
    
    # 存档类型映射 | Archive Type Mapping
    archive_map: Dict[str, str] = None
    
    # 文件名模板 | File Name Templates
    filename_templates: List[str] = None
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        if self.base_dirs is None:
            self.base_dirs = ['INSDC', 'INSDC2', 'INSDC3', 'INSDC4', 'INSDC5']
        
        if self.archive_map is None:
            self.archive_map = {
                "SRR": "SRA",
                "ERR": "ERA", 
                "DRR": "DRA"
            }
        
        if self.filename_templates is None:
            self.filename_templates = [
                "{run_id}.sra",
                "{run_id}_1.fastq.gz", 
                "{run_id}_2.fastq.gz",
                "{run_id}.fastq.gz",
                "{run_id}"
            ]
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        if not self.ftp_host:
            errors.append("❌ FTP主机地址不能为空 | FTP host cannot be empty")
        
        if self.ftp_timeout <= 0:
            errors.append("❌ FTP超时时间必须为正数 | FTP timeout must be positive")
        
        if not self.base_dirs:
            errors.append("❌ 基础目录列表不能为空 | Base directories list cannot be empty")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
