"""
VCF转换工具函数模块|VCF Converter Utility Functions Module
"""

import gzip
import logging
import sys
from pathlib import Path

# IUPAC核苷酸模糊代码字典|Dictionary of IUPAC ambiguities for nucleotides
AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

# 二倍体SNP二进制编码字典|Dictionary for translating biallelic SNPs into binary format
GEN_BIN = {
    "./.":"?", ".|.":"?",
    "0/0":"0", "0|0":"0",
    "0/1":"1", "0|1":"1", "1/0":"1", "1|0":"1",
    "1/1":"2", "1|1":"2",
}

class ConverterLogger:
    """VCF转换日志管理器|VCF Converter Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "vcf_conversion.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class FileHandler:
    """文件处理工具类|File Handler Utility Class"""
    
    @staticmethod
    def get_opener(filename):
        """根据文件扩展名选择合适的打开方式|Choose appropriate opener based on file extension"""
        return gzip.open if filename.lower().endswith(".gz") else open
    
    @staticmethod
    def safe_open(filename, mode="rt", encoding="utf-8"):
        """安全打开文件|Safely open file"""
        opener = FileHandler.get_opener(filename)
        return opener(filename, mode, encoding=encoding)
