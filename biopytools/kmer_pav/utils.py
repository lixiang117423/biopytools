"""
K-mer分析工具函数模块 | K-mer Analysis Utility Functions Module
"""

import logging
import sys
import subprocess
import shutil
import gzip
from pathlib import Path
from typing import List, Union
from Bio import SeqIO

class KmerLogger:
    """K-mer分析日志管理器 | K-mer Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "kmer_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class KMCRunner:
    """KMC工具运行器 | KMC Tool Runner"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self._check_kmc_installation()
    
    def _check_kmc_installation(self):
        """检查KMC工具是否安装 | Check if KMC tools are installed"""
        for tool, path in [("KMC", self.config.kmc_path), ("kmc_tools", self.config.kmc_tools_path)]:
            if not shutil.which(path):
                raise RuntimeError(
                    f"{tool}工具未找到 | {tool} tool not found: {path}\n"
                    f"请安装KMC工具 | Please install KMC tools:\n"
                    f"  Ubuntu/Debian: sudo apt-get install kmc\n"
                    f"  Conda: conda install -c bioconda kmc\n"
                    f"  源码编译 | Build from source: https://github.com/refresh-bio/KMC"
                )
        self.logger.info("KMC工具检查通过 | KMC tools check passed")
    
    def run_command(self, cmd: List[str], description: str = "") -> bool:
        """运行命令 | Run command"""
        if description:
            self.logger.info(f"执行 | Executing: {description}")
        
        cmd_str = " ".join(cmd)
        self.logger.info(f"命令 | Command: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
                # 移除cwd参数，使用绝对路径 | Remove cwd parameter, use absolute paths
            )
            
            if result.stdout.strip():
                self.logger.info(f"输出 | Output: {result.stdout.strip()}")
            
            self.logger.info(f"✓ {description} 完成 | completed")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ {description} 失败 | failed")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False
        except Exception as e:
            self.logger.error(f"✗ 执行命令时发生异常 | Exception occurred while executing command: {e}")
            return False

class FileDetector:
    """文件格式检测器 | File Format Detector"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def detect_format(self, file_path: str) -> str:
        """检测文件格式 | Detect file format"""
        file_path = Path(file_path)
        
        # 检查文件头 | Check file header
        try:
            # 处理压缩文件 | Handle compressed files
            if file_path.suffix.lower() == '.gz':
                with gzip.open(file_path, 'rt') as f:
                    first_char = f.read(1)
            else:
                with open(file_path, 'r') as f:
                    first_char = f.read(1)
            
            if first_char == '>':
                return 'fasta'
            elif first_char == '@':
                return 'fastq'
            else:
                self.logger.warning(f"无法通过文件头识别格式，使用后缀名判断 | Cannot identify format by header, using file extension: {file_path}")
                return self._detect_by_extension(file_path)
                
        except Exception as e:
            self.logger.warning(f"读取文件头失败，使用后缀名判断 | Failed to read file header, using extension: {e}")
            return self._detect_by_extension(file_path)
    
    def _detect_by_extension(self, file_path: Path) -> str:
        """通过文件扩展名检测格式 | Detect format by file extension"""
        # 处理压缩文件 | Handle compressed files
        if file_path.suffix.lower() in ['.gz', '.bz2']:
            actual_suffix = file_path.stem.split('.')[-1].lower()
        else:
            actual_suffix = file_path.suffix[1:].lower()
        
        if actual_suffix in ['fasta', 'fa', 'fas', 'fna']:
            return 'fasta'
        elif actual_suffix in ['fastq', 'fq']:
            return 'fastq'
        else:
            # 默认为fasta | Default to fasta
            self.logger.warning(f"未知文件格式，默认为FASTA | Unknown file format, defaulting to FASTA: {file_path}")
            return 'fasta'

class FastaSequenceSplitter:
    """FASTA序列分割器 | FASTA Sequence Splitter"""
    
    def __init__(self, logger, output_dir: Path):
        self.logger = logger
        self.output_dir = output_dir
        self.seq_dir = output_dir / "individual_sequences"
        self.seq_dir.mkdir(parents=True, exist_ok=True)
    
    def split_fasta_file(self, fasta_file: str) -> List[dict]:
        """将FASTA文件分割成单个序列文件 | Split FASTA file into individual sequence files"""
        self.logger.info(f"分割FASTA文件 | Splitting FASTA file: {fasta_file}")
        
        sequence_files = []
        file_path = Path(fasta_file)
        
        try:
            # 处理压缩文件 | Handle compressed files
            if file_path.suffix.lower() == '.gz':
                opener = gzip.open
                mode = 'rt'
            else:
                opener = open
                mode = 'r'
            
            with opener(fasta_file, mode) as handle:
                for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                    # 清理序列ID，移除特殊字符 | Clean sequence ID, remove special characters
                    clean_id = "".join(c for c in record.id if c.isalnum() or c in "._-")
                    if not clean_id:
                        clean_id = f"seq_{i}"
                    
                    # 创建单序列文件 | Create individual sequence file
                    seq_filename = f"{file_path.stem}_{clean_id}.fasta"
                    seq_filepath = self.seq_dir / seq_filename
                    
                    with open(seq_filepath, 'w') as out_handle:
                        SeqIO.write(record, out_handle, "fasta")
                    
                    sequence_files.append({
                        'sample_name': f"{file_path.stem}_{clean_id}",
                        'file_path': str(seq_filepath),
                        'original_file': fasta_file,
                        'sequence_id': record.id,
                        'sequence_length': len(record.seq)
                    })
            
            self.logger.info(f"分割完成，共 {len(sequence_files)} 条序列 | Splitting completed, {len(sequence_files)} sequences")
            return sequence_files
            
        except Exception as e:
            self.logger.error(f"分割FASTA文件失败 | Failed to split FASTA file {fasta_file}: {e}")
            return []

def find_sequence_files(input_path: str, pattern: str = None) -> List[str]:
    """查找序列文件 | Find sequence files"""
    input_path = Path(input_path)
    files = []
    
    if input_path.is_file():
        return [str(input_path)]
    
    # 如果是目录 | If it's a directory
    if pattern:
        # 使用用户指定的模式 | Use user-specified pattern
        files = list(input_path.glob(pattern))
    else:
        # 默认查找常见的序列文件 | Default search for common sequence files
        patterns = [
            "*.fasta", "*.fa", "*.fas", "*.fna",
            "*.fastq", "*.fq",
            "*.fasta.gz", "*.fa.gz", "*.fas.gz", "*.fna.gz",
            "*.fastq.gz", "*.fq.gz"
        ]
        
        for pattern in patterns:
            files.extend(input_path.glob(pattern))
    
    return sorted([str(f) for f in files])

def cleanup_files(file_patterns: List[str], base_dir: Path, logger):
    """清理临时文件 | Clean up temporary files"""
    for pattern in file_patterns:
        files_to_remove = list(base_dir.glob(pattern))
        for file_path in files_to_remove:
            try:
                if file_path.is_file():
                    file_path.unlink()
                elif file_path.is_dir():
                    shutil.rmtree(file_path)
                logger.debug(f"删除临时文件 | Removed temporary file: {file_path}")
            except Exception as e:
                logger.warning(f"删除文件失败 | Failed to remove file {file_path}: {e}")
