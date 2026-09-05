"""
 ENA下载配置管理模块|ENA Download Configuration Management Module
"""
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

from .utils import expand_input, is_accession, classify_accession

@dataclass
class DownloadConfig:
    """ ENA下载配置类|ENA Download Configuration Class"""

    #  必需参数|Required parameters
    accession: str

    #  输出设置|Output settings
    output_dir: Optional[str] = None
    create_dir: bool = False  #  是否创建目录|Whether to create directory

    #  下载协议设置|Download protocol settings
    protocol: str = "ftp"  #  ftp 或 aspera
    method: str = "save"   #  save 或 run
    aspera_key: Optional[str] = None

    #  元数据设置|Metadata settings
    fields: Optional[List[str]] = None
    metadata_format: str = "tsv"  #  tsv, csv, xlsx

    #  API设置|API settings
    api_url: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    max_retries: int = 3

    #  内部属性|Internal attributes
    base_name: str = ""
    metadata_file: str = ""

    def __post_init__(self):
        """ 初始化后处理|Post-initialization processing"""
        #  展开输入为编号列表(文件或单个编号)|Expand input into accession list (file or single accession)
        self.accessions = self._expand_accessions()

        #  设置输出目录|Set output directory
        if self.create_dir:
            #  创建目录模式|Create directory mode
            if not self.output_dir:
                # 多编号时用通用目录名, 单编号保持原行为|Generic dir name for multiple accessions; keep original behavior for one
                if len(self.accessions) == 1:
                    self.output_dir = f"{self.accessions[0]}.ena.download"
                else:
                    self.output_dir = "ena_batch_download"
            self.output_path = Path(self.output_dir)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        else:
            #  当前目录模式|Current directory mode
            if not self.output_dir:
                self.output_dir = "."
            self.output_path = Path(self.output_dir)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        #  设置基础文件名与元数据文件名(按首个编号, 批量场景由流水线逐编号生成)|Base and metadata names from first accession; batch generates per-accession in pipeline
        self.base_name = self.accessions[0]
        self.metadata_file = f"{self.base_name}_meta.{self.metadata_format}"

        #  验证协议设置|Validate protocol settings
        if self.protocol not in ["ftp", "aspera"]:
            raise ValueError(f" Unsupported protocol: {self.protocol}. Use 'ftp' or 'aspera'")

        if self.method not in ["save", "run"]:
            raise ValueError(f" Unsupported method: {self.method}. Use 'save' or 'run'")

        #  验证aspera设置|Validate aspera settings
        if self.protocol == "aspera":
            if not self.aspera_key:
                raise ValueError(" --aspera_key is required when using aspera protocol")

            key_file = Path(self.aspera_key)
            if not key_file.exists():
                raise ValueError(f" Aspera key file {key_file} does not exist")

            #  检查文件权限|Check file permissions
            if key_file.stat().st_mode & 0o777 != 0o600:
                raise ValueError(
                    f" Key file permissions are insecure (current: {oct(key_file.stat().st_mode & 0o777)}).\n"
                    f" Run: chmod 600 \"{key_file}\" to fix"
                )

    def _expand_accessions(self) -> List[str]:
        """ 展开输入并校验全部编号, 收集所有错误后一次性抛出|Expand input, validate all accessions, and raise once with every error"""
        raw_list = expand_input(str(self.accession).strip())
        if not raw_list:
            raise ValueError(
                f" 输入文件无有效编号(每行一个, 支持#注释)|Input file has no valid accessions (one per line, '#' comments supported): {self.accession}"
            )

        errors = []
        for entry in raw_list:
            if not is_accession(entry):
                errors.append(f" 非法ENA编号格式|Invalid ENA accession format: {entry}")
        if errors:
            raise ValueError("\n".join(errors))

        return raw_list

    def log_input_summary(self, logger):
        """ 记录输入概要(数量/来源/类型), 排错用|Log input summary (count/source/type) for troubleshooting"""
        source = "ID文件|ID file" if len(self.accessions) > 1 else "单编号|single input"
        logger.info(f" 输入编号数|Accession count: {len(self.accessions)} ({source})")
        for entry in self.accessions:
            logger.info(f" 编号|Accession: {entry} ({classify_accession(entry)})")

    def validate(self):
        """ 验证配置|Validate configuration"""
        if not self.accession:
            raise ValueError(" Accession number is required")

        if not getattr(self, 'accessions', None):
            raise ValueError(" Accession list is empty")

        #  检查元数据格式|Check metadata format
        if self.metadata_format not in ["tsv", "csv", "xlsx"]:
            raise ValueError(f" Unsupported metadata format: {self.metadata_format}")

        return True
