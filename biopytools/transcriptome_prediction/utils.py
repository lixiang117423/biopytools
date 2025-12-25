"""
🛠️ 转录组预测分析工具函数模块 | Transcriptome Prediction Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import re
from pathlib import Path
from typing import List, Tuple, Dict

class TranscriptomeLogger:
    """📝 转录组分析日志管理器 | Transcriptome Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "transcriptome_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """🔧 设置日志 | Setup logging"""
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
        """📋 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """⚡ 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()  # 使用绝对路径 | Use absolute path
    
    def run(self, cmd: str, description: str = "") -> bool:
        """🚀 执行命令 | Execute command"""
        if description:
            self.logger.info(f"🔄 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"⚠️ 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False

class SequencingTypeDetector:
    """🔍 测序类型检测器 | Sequencing Type Detector"""
    
    @staticmethod
    def detect_sequencing_layout(rna_seq_files: List[str]) -> Tuple[str, List[Tuple[str, ...]]]:
        """
        🧬 检测测序布局和样本组织
        Detect sequencing layout and sample organization
        
        Returns:
            layout: 'paired' or 'single'
            samples: [(sample1_R1, sample1_R2), ...] for paired or [(sample1,), ...] for single
        """
        
        # 按文件名分析样本 | Analyze samples by filename
        samples_dict = {}
        
        for file_path in rna_seq_files:
            file_name = Path(file_path).name
            
            # 检测配对末端模式 | Detect paired-end patterns
            # 常见模式：_R1/_R2, _1/_2, .1/.2, _F/_R
            paired_patterns = [
                (r'(.+)[._]R?1[._]', r'(.+)[._]R?2[._]'),
                (r'(.+)[._]F[._]', r'(.+)[._]R[._]'),
                (r'(.+)\.1\.', r'(.+)\.2\.'),
                (r'(.+)_1\.', r'(.+)_2\.'),
                (r'(.+)_forward', r'(.+)_reverse'),
                (r'(.+)_fwd', r'(.+)_rev')
            ]
            
            sample_base = None
            read_type = None
            
            for r1_pattern, r2_pattern in paired_patterns:
                r1_match = re.search(r1_pattern, file_name, re.IGNORECASE)
                r2_match = re.search(r2_pattern, file_name, re.IGNORECASE)
                
                if r1_match:
                    sample_base = r1_match.group(1)
                    read_type = 'R1'
                    break
                elif r2_match:
                    sample_base = r2_match.group(1)
                    read_type = 'R2'
                    break
            
            if sample_base is None:
                # 单端模式，使用整个文件名作为样本名 | Single-end mode, use whole filename as sample name
                sample_base = re.sub(r'\.(fq|fastq|fa|fasta)(\.gz|\.bz2)?$', '', file_name, flags=re.IGNORECASE)
                read_type = 'single'
            
            if sample_base not in samples_dict:
                samples_dict[sample_base] = {}
            
            samples_dict[sample_base][read_type] = file_path
        
        # 分析样本组织 | Analyze sample organization
        paired_samples = []
        single_samples = []
        
        for sample_base, reads in samples_dict.items():
            if 'R1' in reads and 'R2' in reads:
                paired_samples.append((reads['R1'], reads['R2']))
            elif 'single' in reads:
                single_samples.append((reads['single'],))
            elif 'R1' in reads or 'R2' in reads:
                # 只有一个文件但匹配了配对模式，当作单端处理 | Only one file matched paired pattern, treat as single
                read_file = reads.get('R1') or reads.get('R2')
                single_samples.append((read_file,))
        
        # 决定整体布局 | Determine overall layout
        if paired_samples and not single_samples:
            return 'paired', paired_samples
        elif single_samples and not paired_samples:
            return 'single', single_samples
        elif paired_samples and single_samples:
            # 混合模式，优先配对末端 | Mixed mode, prioritize paired-end
            return 'mixed', paired_samples + single_samples
        else:
            # 降级为单端模式 | Fallback to single-end mode
            all_files = [(f,) for f in rna_seq_files]
            return 'single', all_files

def check_dependencies(config, logger):
    """🔧 检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.hisat2_path, "HISAT2", "🧬"),
        (config.hisat2_build_path, "HISAT2-build", "🏗️"),
        (config.stringtie_path, "StringTie", "🧩"),
        (config.trinity_path, "Trinity", "🔗"),
        (config.samtools_path, "SAMtools", "🔧")
    ]
    
    missing_deps = []
    
    for cmd, name, emoji in dependencies:
        try:
            # 不同工具的版本检查命令不同 | Different tools have different version check commands
            if name == "HISAT2":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "StringTie":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "Trinity":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "SAMtools":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            else:
                result = subprocess.run([cmd, "--help"], 
                                      capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                logger.info(f"✅ {emoji} {name} 可用 | available")
            else:
                missing_deps.append(f"{emoji} {name}")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(f"{emoji} {name}")
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True
