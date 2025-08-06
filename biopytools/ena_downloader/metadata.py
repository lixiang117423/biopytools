"""
ENA元数据处理模块 | ENA Metadata Processing Module
"""

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pathlib import Path
from typing import Optional, List, Dict, Any
import warnings

# 默认ENA API配置 | Default ENA API configuration
DEFAULT_ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DEFAULT_FIELDS = {
    "run_accession",
    "study_accession", 
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "submission_accession",
    "tax_id",
    "scientific_name",
    "instrument_model",
    "nominal_length",
    "library_layout",
    "library_source",
    "library_selection",
    "base_count",
    "first_public",
    "last_updated",
    "study_title",
    "experiment_alias",
    "run_alias",
    "fastq_bytes",
    "fastq_md5",
    "fastq_ftp",
    "fastq_aspera",
    "fastq_galaxy",
    "submitted_bytes",
    "submitted_md5",
    "submitted_ftp",
    "submitted_galaxy",
    "submitted_format",
    "sra_bytes",
    "sra_md5",
    "sra_ftp",
    "sample_alias",
    "broker_name",
    "sample_title",
    "nominal_sdev",
    "bam_ftp",
    "bam_bytes",
}

class MetadataDownloader:
    """ENA元数据下载器 | ENA Metadata Downloader"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.session = self._create_session()
    
    def _create_session(self) -> requests.Session:
        """创建HTTP会话 | Create HTTP session"""
        session = requests.Session()
        retries = Retry(
            total=self.config.max_retries,
            backoff_factor=0.5,
            status_forcelist=[500, 502, 503, 504]
        )
        session.mount('https://', HTTPAdapter(max_retries=retries))
        return session
    
    def _prepare_fields(self) -> str:
        """准备API字段参数 | Prepare API fields parameter"""
        if self.config.fields and 'all' not in self.config.fields:
            # 确保包含run_accession | Ensure run_accession is included
            fields_set = {'run_accession'} | set(self.config.fields)
            return ','.join(sorted(fields_set))
        else:
            # 使用所有默认字段 | Use all default fields
            return ','.join(sorted(DEFAULT_FIELDS))
    
    def download_metadata(self) -> Optional[Path]:
        """下载元数据文件 | Download metadata file"""
        self.logger.info(f"开始下载元数据 | Starting metadata download for accession: {self.config.accession}")
        
        # 准备API参数 | Prepare API parameters
        params = {
            'accession': self.config.accession,
            'result': 'read_run',
            'format': 'tsv',
            'fields': self._prepare_fields()
        }
        
        # 确定输出文件路径 | Determine output file path
        output_file = self.config.output_path / self.config.metadata_file
        
        try:
            # 发送API请求 | Send API request
            self.logger.info(f"正在请求ENA API | Requesting ENA API: {self.config.api_url}")
            response = self.session.get(
                self.config.api_url, 
                params=params, 
                timeout=30
            )
            response.raise_for_status()
            
            # 处理响应内容 | Process response content
            file_content = response.text
            
            if not file_content.strip():
                self.logger.error("API返回空内容 | API returned empty content")
                return None
            
            # 根据格式保存文件 | Save file according to format
            if self.config.metadata_format == 'xlsx':
                output_file = self._save_as_excel(file_content, output_file)
            elif self.config.metadata_format == 'csv':
                output_file = self._save_as_csv(file_content, output_file)
            else:  # tsv
                output_file = self._save_as_tsv(file_content, output_file)
            
            self.logger.info(f"元数据文件已保存 | Metadata file saved: {output_file}")
            return output_file
            
        except requests.exceptions.RequestException as e:
            self.logger.error(f"下载元数据失败 | Failed to download metadata: {str(e)}")
            return None
        except IOError as e:
            self.logger.error(f"保存文件失败 | Failed to save file: {str(e)}")
            return None
    
    def _save_as_tsv(self, content: str, output_file: Path) -> Path:
        """保存为TSV格式 | Save as TSV format"""
        output_file.write_text(content, encoding='utf-8')
        return output_file
    
    def _save_as_csv(self, content: str, output_file: Path) -> Path:
        """保存为CSV格式 | Save as CSV format"""
        csv_content = content.replace('\t', ',')
        csv_file = output_file.with_suffix('.csv')
        csv_file.write_text(csv_content, encoding='utf-8')
        return csv_file
    
    def _save_as_excel(self, content: str, output_file: Path) -> Path:
        """保存为Excel格式 | Save as Excel format"""
        try:
            import pandas as pd
            import numpy as np
            warnings.simplefilter('ignore', category=FutureWarning)
            
            # 转换TSV为DataFrame | Convert TSV to DataFrame
            lines = content.strip().split('\n')
            if not lines:
                raise ValueError("Empty content")
            
            # 解析表头和数据 | Parse header and data
            header = lines[0].split('\t')
            data_rows = [line.split('\t') for line in lines[1:] if line.strip()]
            
            df = pd.DataFrame(data_rows, columns=header)
            
            # 处理空值 | Handle empty values
            df.replace('', np.nan, inplace=True)
            
            # 保存为Excel | Save as Excel
            excel_file = output_file.with_suffix('.xlsx')
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                df.to_excel(writer, index=False, sheet_name='ENA_Metadata')
            
            return excel_file
            
        except ImportError:
            self.logger.warning("pandas未安装，无法生成Excel文件，回退到TSV格式 | pandas not installed, cannot generate Excel file, fallback to TSV")
            return self._save_as_tsv(content, output_file.with_suffix('.tsv'))
    
    def get_download_links(self, metadata_file: Path) -> List[str]:
        """从元数据文件提取下载链接 | Extract download links from metadata file"""
        self.logger.info(f"从元数据文件提取下载链接 | Extracting download links from metadata file: {metadata_file}")
        
        try:
            with open(metadata_file, 'r', encoding='utf-8') as f:
                header = f.readline().strip().split('\t')
                
                # 确定链接列 | Determine link column
                link_column = 'fastq_aspera' if self.config.protocol == 'aspera' else 'fastq_ftp'
                
                if link_column not in header:
                    self.logger.error(f"元数据文件中未找到 {link_column} 列 | Column {link_column} not found in metadata file")
                    return []
                
                column_index = header.index(link_column)
                
                # 提取所有链接 | Extract all links
                links = []
                for line in f:
                    if line.strip():
                        columns = line.strip().split('\t')
                        if column_index < len(columns) and columns[column_index]:
                            # 处理多个链接（分号分隔） | Handle multiple links (semicolon separated)
                            for link in columns[column_index].split(';'):
                                if link.strip():
                                    links.append(link.strip())
                
                self.logger.info(f"找到 {len(links)} 个下载链接 | Found {len(links)} download links")
                return links
                
        except FileNotFoundError:
            self.logger.error(f"元数据文件未找到 | Metadata file not found: {metadata_file}")
            return []
        except Exception as e:
            self.logger.error(f"读取元数据文件失败 | Failed to read metadata file: {str(e)}")
            return []
