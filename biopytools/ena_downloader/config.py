"""
ENA下载配置管理模块 | ENA Download Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class DownloadConfig:
    """ENA下载配置类 | ENA Download Configuration Class"""
    
    # 必需参数 | Required parameters
    accession: str
    
    # 输出设置 | Output settings
    output_dir: Optional[str] = None
    create_dir: bool = False  # 是否创建目录 | Whether to create directory
    
    # 下载协议设置 | Download protocol settings
    protocol: str = "ftp"  # ftp 或 aspera
    method: str = "save"   # save 或 run
    aspera_key: Optional[str] = None
    
    # 元数据设置 | Metadata settings
    fields: Optional[List[str]] = None
    metadata_format: str = "tsv"  # tsv, csv, xlsx
    
    # API设置 | API settings
    api_url: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    max_retries: int = 3
    
    # 内部属性 | Internal attributes
    base_name: str = ""
    metadata_file: str = ""
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 设置输出目录 | Set output directory
        if self.create_dir:
            # 创建目录模式 | Create directory mode
            if not self.output_dir:
                self.output_dir = f"{self.accession}.ena.download"
            self.output_path = Path(self.output_dir)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        else:
            # 当前目录模式 | Current directory mode
            if not self.output_dir:
                self.output_dir = "."
            self.output_path = Path(self.output_dir)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 设置基础文件名 | Set base filename
        self.base_name = self.accession
        # 移除文件名前的点 | Remove dot prefix from filename
        self.metadata_file = f"{self.accession}.meta.{self.metadata_format}"
        
        # 验证协议设置 | Validate protocol settings
        if self.protocol not in ["ftp", "aspera"]:
            raise ValueError(f"Unsupported protocol: {self.protocol}. Use 'ftp' or 'aspera'")
        
        if self.method not in ["save", "run"]:
            raise ValueError(f"Unsupported method: {self.method}. Use 'save' or 'run'")
        
        # 验证aspera设置 | Validate aspera settings
        if self.protocol == "aspera":
            if not self.aspera_key:
                raise ValueError("--aspera_key is required when using aspera protocol")
            
            key_file = Path(self.aspera_key)
            if not key_file.exists():
                raise ValueError(f"Aspera key file {key_file} does not exist")
            
            # 检查文件权限 | Check file permissions
            if key_file.stat().st_mode & 0o777 != 0o600:
                raise ValueError(
                    f"Key file permissions are insecure (current: {oct(key_file.stat().st_mode & 0o777)}).\n"
                    f"Run: chmod 600 \"{key_file}\" to fix"
                )
    
    def validate(self):
        """验证配置 | Validate configuration"""
        if not self.accession:
            raise ValueError("Accession number is required")
        
        # 检查元数据格式 | Check metadata format
        if self.metadata_format not in ["tsv", "csv", "xlsx"]:
            raise ValueError(f"Unsupported metadata format: {self.metadata_format}")
