# """
# 🛠️ GTX WGS分析工具函数模块 | GTX WGS Analysis Utility Functions Module 🛠️
# """

# import logging
# import subprocess
# import sys
# import os
# import re
# from pathlib import Path
# from datetime import datetime

# class GTXLogger:
#     """📝 GTX分析日志管理器 | GTX Analysis Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "gtx_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 ✍️ | Setup logging"""
#         if self.log_file.exists():
#             # 备份现有日志文件 💾 | Backup existing log file
#             backup_name = f"gtx_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
#             self.log_file.rename(self.output_dir / backup_name)
        
#         logging.basicConfig(
#             level=logging.INFO,
#             format='[%(asctime)s] %(levelname)s - %(message)s',
#             datefmt='%Y-%m-%d %H:%M:%S',
#             handlers=[
#                 logging.FileHandler(self.log_file),
#                 logging.StreamHandler(sys.stdout)
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def get_logger(self):
#         """获取日志器 📢 | Get logger"""
#         return self.logger

# class CommandRunner:
#     """🚀 命令执行器 | Command Runner"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir.resolve()
    
#     def run(self, cmd: str, description: str = "") -> bool:
#         """执行命令 ▶️ | Execute command"""
#         if description:
#             self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
#         self.logger.info(f"💻 命令 | Command: {cmd}")
#         self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
#         try:
#             result = subprocess.run(
#                 cmd, 
#                 shell=True, 
#                 capture_output=True, 
#                 text=True, 
#                 check=True,
#                 cwd=self.working_dir
#             )
            
#             self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
#             if result.stdout:
#                 self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
#             return True
            
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
#             self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
#             self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
#             self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
#             return False

# class PatternParser:
#     """📋 文件模式解析器 | File Pattern Parser"""
    
#     def __init__(self, logger):
#         self.logger = logger
    
#     def extract_sample_name_from_pattern(self, filename: str, pattern: str) -> str:
#         """根据模式从文件名提取样品名 | Extract sample name from filename based on pattern"""
#         import re
        
#         try:
#             # 先转义所有特殊字符，但保留*作为占位符
#             escaped_pattern = re.escape(pattern)
#             # 现在把转义后的\*替换为捕获组
#             regex_pattern = escaped_pattern.replace(r'\*', '(.+?)')
#             # 添加行首行尾锚点
#             regex_pattern = f"^{regex_pattern}$"
            
#             self.logger.debug(f"模式转换 | Pattern conversion: {pattern} -> {regex_pattern}")
            
#             match = re.match(regex_pattern, filename)
#             if match:
#                 sample_name = match.group(1)
#                 self.logger.debug(f"成功提取样品名 | Successfully extracted sample name: {filename} -> {sample_name}")
#                 return sample_name
#             else:
#                 self.logger.debug(f"文件名 {filename} 不匹配模式 {pattern}")
#                 return None
#         except Exception as e:
#             self.logger.error(f"模式匹配错误 | Pattern matching error: {e}")
#             return None
    
#     def build_paired_filename(self, sample_name: str, pattern: str) -> str:
#         """根据样品名和模式构建配对文件名 | Build paired filename based on sample name and pattern"""
#         return pattern.replace("*", sample_name)

# class FileProcessor:
#     """📁 文件处理器 | File Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
#         self.pattern_parser = PatternParser(logger)
    
#     def find_fastq_files(self):
#         """查找FASTQ文件对 🔍 | Find FASTQ file pairs"""
#         self.logger.info("🔍 搜索输入文件 | Searching input files")
        
#         input_path = Path(self.config.input_dir)
#         r1_files = list(input_path.glob(self.config.read1_pattern))
        
#         if not r1_files:
#             raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到 {self.config.read1_pattern} 文件")
        
#         # 按文件名排序 🔠 | Sort by filename
#         r1_files.sort()
        
#         self.logger.info(f"找到 {len(r1_files)} 个R1文件 | Found {len(r1_files)} R1 files")
#         for r1_file in r1_files[:5]:  # 显示前5个文件作为示例
#             self.logger.info(f"  R1文件示例 | R1 file example: {r1_file.name}")
#         if len(r1_files) > 5:
#             self.logger.info(f"  ... 还有 {len(r1_files) - 5} 个文件 | ... and {len(r1_files) - 5} more files")
        
#         # 验证R2文件存在 ✅ | Validate R2 files exist
#         file_pairs = []
#         for r1_file in r1_files:
#             # 从R1文件名中提取样品名 🏷️ | Extract sample name from R1 filename
#             sample_name = self.pattern_parser.extract_sample_name_from_pattern(
#                 r1_file.name, self.config.read1_pattern
#             )
            
#             if not sample_name:
#                 self.logger.warning(f"⚠️ 无法从文件名提取样品名 | Cannot extract sample name from filename: {r1_file.name}")
#                 continue
            
#             # 构建R2文件路径 🏗️ | Build R2 file path
#             r2_filename = self.pattern_parser.build_paired_filename(sample_name, self.config.read2_pattern)
#             r2_file = input_path / r2_filename
            
#             if not r2_file.exists():
#                 self.logger.warning(f"⚠️ 找不到对应的R2文件 | Cannot find corresponding R2 file: {r2_file}")
#                 continue
            
#             file_pairs.append((sample_name, str(r1_file), str(r2_file)))
#             self.logger.debug(f"✅ 找到配对文件 | Found paired files: {sample_name}")
        
#         total_samples = len(file_pairs)
#         self.logger.info(f"✅ 找到 {total_samples} 个样品需要处理 | Found {total_samples} samples to process")
        
#         # 显示样品名示例 | Show sample name examples
#         if total_samples > 0:
#             self.logger.info("样品列表示例 | Sample list examples:")
#             for i, (sample_name, r1_file, r2_file) in enumerate(file_pairs[:3]):
#                 self.logger.info(f"  {i+1}. {sample_name}")
#                 self.logger.info(f"     R1: {Path(r1_file).name}")
#                 self.logger.info(f"     R2: {Path(r2_file).name}")
#             if total_samples > 3:
#                 self.logger.info(f"  ... 还有 {total_samples - 3} 个样品 | ... and {total_samples - 3} more samples")
        
#         return file_pairs
    
#     def check_output_exists(self, sample_name: str) -> bool:
#         """检查输出文件是否已存在 🧐 | Check if output files already exist"""
#         vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
#         bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
#         return vcf_file.exists() and bam_file.exists()
    
#     def get_file_size(self, file_path: str) -> str:
#         """获取文件大小 📏 | Get file size"""
#         try:
#             size_bytes = os.path.getsize(file_path)
#             # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
#             for unit in ['B', 'KB', 'MB', 'GB']:
#                 if size_bytes < 1024.0:
#                     return f"{size_bytes:.1f} {unit}"
#                 size_bytes /= 1024.0
#             return f"{size_bytes:.1f} TB"
#         except:
#             return "Unknown"

# def check_dependencies(config, logger):
#     """检查依赖软件 🧩 | Check dependencies"""
#     logger.info("🧩 检查依赖软件 | Checking dependencies")
    
#     # 检查GTX程序 💻 | Check GTX program
#     if not os.path.exists(config.gtx_path):
#         error_msg = f"❌ GTX程序不存在 | GTX program does not exist: {config.gtx_path}"
#         logger.error(error_msg)
#         raise RuntimeError(error_msg)
    
#     # 检查GTX程序是否可执行 🏃 | Check if GTX program is executable
#     if not os.access(config.gtx_path, os.X_OK):
#         error_msg = f"❌ GTX程序不可执行 | GTX program is not executable: {config.gtx_path}"
#         logger.error(error_msg)
#         raise RuntimeError(error_msg)
    
#     logger.info("✅ ✓ GTX程序检查通过 | GTX program check passed")
#     return True

"""
🛠️ GTX WGS分析工具函数模块 | GTX WGS Analysis Utility Functions Module 🛠️
"""

import logging
import subprocess
import sys
import os
import re
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional

class GTXLogger:
    """📝 GTX分析日志管理器 | GTX Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "gtx_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 ✍️ | Setup logging"""
        if self.log_file.exists():
            # 备份现有日志文件 💾 | Backup existing log file
            backup_name = f"gtx_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
            self.log_file.rename(self.output_dir / backup_name)
        
        logging.basicConfig(
            level=logging.INFO,
            format='[%(asctime)s] %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 📢 | Get logger"""
        return self.logger

class CommandRunner:
    """🚀 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 ▶️ | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
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
                self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
            return False

class PatternParser:
    """📋 文件模式解析器 | File Pattern Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def extract_sample_name_from_pattern(self, filename: str, pattern: str) -> str:
        """根据模式从文件名提取样品名 | Extract sample name from filename based on pattern"""
        try:
            # 先转义所有特殊字符，但保留*作为占位符
            escaped_pattern = re.escape(pattern)
            # 现在把转义后的\*替换为捕获组
            regex_pattern = escaped_pattern.replace(r'\*', '(.+?)')
            # 添加行首行尾锚点
            regex_pattern = f"^{regex_pattern}$"
            
            self.logger.debug(f"模式转换 | Pattern conversion: {pattern} -> {regex_pattern}")
            
            match = re.match(regex_pattern, filename)
            if match:
                sample_name = match.group(1)
                self.logger.debug(f"成功提取样品名 | Successfully extracted sample name: {filename} -> {sample_name}")
                return sample_name
            else:
                self.logger.debug(f"文件名 {filename} 不匹配模式 {pattern}")
                return None
        except Exception as e:
            self.logger.error(f"模式匹配错误 | Pattern matching error: {e}")
            return None
    
    def build_paired_filename(self, sample_name: str, pattern: str) -> str:
        """根据样品名和模式构建配对文件名 | Build paired filename based on sample name and pattern"""
        return pattern.replace("*", sample_name)

class VCFValidator:
    """📜 VCF文件验证器 | VCF File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def is_valid_vcf(self, vcf_file: Path) -> bool:
        """检查VCF文件是否有效 ✅ | Check if VCF file is valid"""
        if not vcf_file.exists():
            self.logger.debug(f"VCF文件不存在 | VCF file does not exist: {vcf_file}")
            return False
        
        if vcf_file.stat().st_size == 0:
            self.logger.debug(f"VCF文件为空 | VCF file is empty: {vcf_file}")
            return False
        
        try:
            # 检查是否为压缩文件并尝试读取头部
            if str(vcf_file).endswith('.gz'):
                import gzip
                with gzip.open(vcf_file, 'rt') as f:
                    first_line = f.readline()
            else:
                with open(vcf_file, 'r') as f:
                    first_line = f.readline()
            
            # 检查VCF格式头部
            if first_line.startswith('##fileformat=VCF'):
                self.logger.debug(f"VCF文件格式有效 | VCF file format valid: {vcf_file}")
                return True
            else:
                self.logger.debug(f"VCF文件格式无效 | Invalid VCF file format: {vcf_file}")
                return False
                
        except Exception as e:
            self.logger.debug(f"VCF文件验证失败 | VCF file validation failed: {vcf_file}, error: {e}")
            return False
    
    def get_vcf_sample_count(self, vcf_file: Path) -> int:
        """获取VCF文件中的样品数量 📊 | Get sample count in VCF file"""
        try:
            if str(vcf_file).endswith('.gz'):
                import gzip
                open_func = gzip.open
                mode = 'rt'
            else:
                open_func = open
                mode = 'r'
            
            with open_func(vcf_file, mode) as f:
                for line in f:
                    if line.startswith('#CHROM'):
                        # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
                        columns = line.strip().split('\t')
                        sample_count = len(columns) - 9  # 减去前9列固定列
                        self.logger.debug(f"VCF文件包含样品数 | VCF file contains samples: {sample_count}")
                        return sample_count
            
            self.logger.debug(f"未找到样品信息行 | Sample info line not found in: {vcf_file}")
            return 0
            
        except Exception as e:
            self.logger.debug(f"获取VCF样品数失败 | Failed to get VCF sample count: {vcf_file}, error: {e}")
            return 0

class JointValidator:
    """🤝 Joint Calling验证器 | Joint Calling Validator"""
    
    def __init__(self, logger):
        self.logger = logger
        self.vcf_validator = VCFValidator(logger)
    
    def validate_sample_map_format(self, sample_map_file: Path) -> bool:
        """验证样品映射文件格式 📋 | Validate sample mapping file format"""
        if not sample_map_file.exists():
            self.logger.debug(f"样品映射文件不存在 | Sample mapping file does not exist: {sample_map_file}")
            return False
        
        try:
            valid_entries = 0
            with open(sample_map_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:  # 跳过空行
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) != 2:
                        self.logger.warning(f"样品映射文件第{line_num}行格式错误 | Invalid format at line {line_num}: {line}")
                        continue
                    
                    sample_name, vcf_path = parts
                    vcf_file = Path(vcf_path)
                    
                    if not self.vcf_validator.is_valid_vcf(vcf_file):
                        self.logger.warning(f"无效的VCF文件 | Invalid VCF file: {vcf_path}")
                        continue
                    
                    valid_entries += 1
            
            self.logger.debug(f"样品映射文件包含有效条目 | Sample mapping file contains valid entries: {valid_entries}")
            return valid_entries >= 2  # 至少需要2个样品才能进行joint calling
            
        except Exception as e:
            self.logger.error(f"验证样品映射文件失败 | Failed to validate sample mapping file: {e}")
            return False
    
    def estimate_joint_memory_requirement(self, vcf_files: List[Path]) -> str:
        """估算Joint Calling内存需求 💾 | Estimate Joint Calling memory requirement"""
        try:
            total_size = sum(f.stat().st_size for f in vcf_files if f.exists())
            # 经验公式：Joint calling需要约为输入VCF文件总大小的3-5倍内存
            estimated_memory_gb = (total_size * 4) / (1024**3)  # 转换为GB
            
            if estimated_memory_gb < 4:
                return "建议4GB以上内存 | Recommend 4GB+ RAM"
            elif estimated_memory_gb < 16:
                return f"建议{int(estimated_memory_gb*1.5)}GB以上内存 | Recommend {int(estimated_memory_gb*1.5)}GB+ RAM"
            else:
                return f"建议{int(estimated_memory_gb*1.2)}GB以上内存 | Recommend {int(estimated_memory_gb*1.2)}GB+ RAM"
                
        except Exception as e:
            self.logger.debug(f"估算内存需求失败 | Failed to estimate memory requirement: {e}")
            return "建议16GB以上内存 | Recommend 16GB+ RAM"

class ReferenceIndexManager:
    """🧬 参考基因组索引管理器 | Reference Genome Index Manager"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def check_reference_index(self) -> bool:
        """检查参考基因组索引是否存在"""
        reference_path = Path(self.config.reference)
        expected_index_files = [
            f"{reference_path}.fai", 
            f"{reference_path}.bwt", 
            f"{reference_path}.pac", 
            f"{reference_path}.ann", 
            f"{reference_path}.sa",
            f"{reference_path}.amb",
            f"{reference_path}.fai"

        ]
        
        missing_files = []
        for index_file in expected_index_files:
            if not Path(index_file).exists():
                missing_files.append(Path(index_file).name)
        
        if missing_files:
            self.logger.info(f"🔍 缺失索引文件 | Missing index files: {', '.join(missing_files)}")
            return False
        else:
            self.logger.info("✅ 参考基因组索引完整 | Reference index complete")
            return True
    
    def build_reference_index(self) -> bool:
        """构建参考基因组索引"""
        self.logger.info("🔨 构建参考基因组索引 | Building reference index")
        
        reference_path = self.config.reference
        if not Path(reference_path).exists():
            self.logger.error(f"❌ 参考基因组文件不存在: {reference_path}")
            return False
        
        # 预估时间
        file_size_gb = Path(reference_path).stat().st_size / (1024**3)
        if file_size_gb > 1:
            estimated_time = int(file_size_gb * 5)
            self.logger.info(f"⏱️  文件大小: {file_size_gb:.1f} GB，预计需要 ~ {estimated_time} 分钟")
        
        # cmd = f"{self.config.gtx_path} index {reference_path}"
        cmd = f"bwa index {reference_path}"
        success = self.cmd_runner.run(cmd, "BWA构建参考基因组索引")

        cmd = f"samtools faidx {reference_path}"
        success = self.cmd_runner.run(cmd, "SAMtools构建参考基因组索引")
        
        if success:
            self.logger.info("✅ 索引构建完成")
            # 再次验证索引
            return self.check_reference_index()
        else:
            self.logger.error("❌ 索引构建失败")
            return False
    
    def ensure_reference_index(self) -> bool:
        """确保参考基因组索引存在"""
        if self.check_reference_index():
            return True
        else:
            self.logger.info("🚀 开始自动构建缺失的索引文件")
            return self.build_reference_index()

class FileProcessor:
    """📁 文件处理器 | File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.pattern_parser = PatternParser(logger)
        self.vcf_validator = VCFValidator(logger)
        self.joint_validator = JointValidator(logger)
        self.reference_manager = None  # 延迟初始化
    
    def find_fastq_files(self):
        """查找FASTQ文件对 🔍 | Find FASTQ file pairs"""
        self.logger.info("🔍 搜索输入文件 | Searching input files")
        
        input_path = Path(self.config.input_dir)
        r1_files = list(input_path.glob(self.config.read1_pattern))
        
        if not r1_files:
            raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到 {self.config.read1_pattern} 文件")
        
        # 按文件名排序 🔠 | Sort by filename
        r1_files.sort()
        
        self.logger.info(f"找到 {len(r1_files)} 个R1文件 | Found {len(r1_files)} R1 files")
        for r1_file in r1_files[:5]:  # 显示前5个文件作为示例
            self.logger.info(f"  R1文件示例 | R1 file example: {r1_file.name}")
        if len(r1_files) > 5:
            self.logger.info(f"  ... 还有 {len(r1_files) - 5} 个文件 | ... and {len(r1_files) - 5} more files")
        
        # 验证R2文件存在 ✅ | Validate R2 files exist
        file_pairs = []
        for r1_file in r1_files:
            # 从R1文件名中提取样品名 🏷️ | Extract sample name from R1 filename
            sample_name = self.pattern_parser.extract_sample_name_from_pattern(
                r1_file.name, self.config.read1_pattern
            )
            
            if not sample_name:
                self.logger.warning(f"⚠️ 无法从文件名提取样品名 | Cannot extract sample name from filename: {r1_file.name}")
                continue
            
            # 构建R2文件路径 🏗️ | Build R2 file path
            r2_filename = self.pattern_parser.build_paired_filename(sample_name, self.config.read2_pattern)
            r2_file = input_path / r2_filename
            
            if not r2_file.exists():
                self.logger.warning(f"⚠️ 找不到对应的R2文件 | Cannot find corresponding R2 file: {r2_file}")
                continue
            
            file_pairs.append((sample_name, str(r1_file), str(r2_file)))
            self.logger.debug(f"✅ 找到配对文件 | Found paired files: {sample_name}")
        
        total_samples = len(file_pairs)
        self.logger.info(f"✅ 找到 {total_samples} 个样品需要处理 | Found {total_samples} samples to process")
        
        # 显示样品名示例 | Show sample name examples
        if total_samples > 0:
            self.logger.info("样品列表示例 | Sample list examples:")
            for i, (sample_name, r1_file, r2_file) in enumerate(file_pairs[:3]):
                self.logger.info(f"  {i+1}. {sample_name}")
                self.logger.info(f"     R1: {Path(r1_file).name}")
                self.logger.info(f"     R2: {Path(r2_file).name}")
            if total_samples > 3:
                self.logger.info(f"  ... 还有 {total_samples - 3} 个样品 | ... and {total_samples - 3} more samples")
        
        return file_pairs
    
    def check_output_exists(self, sample_name: str) -> bool:
        """检查输出文件是否已存在 🧐 | Check if output files already exist"""
        vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
        bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        return vcf_file.exists() and bam_file.exists()
    
    def get_processed_vcf_files(self) -> List[Path]:
        """获取已处理的VCF文件列表 📜 | Get list of processed VCF files"""
        vcf_files = list(self.config.vcf_output_dir.glob("*.vcf.gz"))
        valid_vcf_files = []
        
        for vcf_file in vcf_files:
            if self.vcf_validator.is_valid_vcf(vcf_file):
                valid_vcf_files.append(vcf_file)
            else:
                self.logger.warning(f"⚠️ 跳过无效VCF文件 | Skipping invalid VCF file: {vcf_file}")
        
        return sorted(valid_vcf_files)
    
    def validate_joint_prerequisites(self) -> Tuple[bool, str]:
        """验证Joint Calling先决条件 🔍 | Validate Joint Calling prerequisites"""
        vcf_files = self.get_processed_vcf_files()
        
        if len(vcf_files) < 2:
            return False, f"需要至少2个有效VCF文件，当前只有{len(vcf_files)}个"
        
        # 检查内存需求
        memory_info = self.joint_validator.estimate_joint_memory_requirement(vcf_files)
        self.logger.info(f"💾 内存需求估算 | Memory requirement estimate: {memory_info}")
        
        # 验证样品映射文件（如果存在）
        if hasattr(self.config, 'sample_map_file') and self.config.sample_map_file.exists():
            if not self.joint_validator.validate_sample_map_format(self.config.sample_map_file):
                return False, "样品映射文件格式无效"
        
        return True, f"找到{len(vcf_files)}个有效VCF文件，满足Joint Calling条件"
    
    def get_file_size(self, file_path: str) -> str:
        """获取文件大小 📏 | Get file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size_bytes < 1024.0:
                    return f"{size_bytes:.1f} {unit}"
                size_bytes /= 1024.0
            return f"{size_bytes:.1f} TB"
        except:
            return "Unknown"

def check_dependencies(config, logger):
    """检查依赖软件 🧩 | Check dependencies"""
    logger.info("🧩 检查依赖软件 | Checking dependencies")
    
    # 检查GTX程序 💻 | Check GTX program
    if not os.path.exists(config.gtx_path):
        error_msg = f"❌ GTX程序不存在 | GTX program does not exist: {config.gtx_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    # 检查GTX程序是否可执行 🏃 | Check if GTX program is executable
    if not os.access(config.gtx_path, os.X_OK):
        error_msg = f"❌ GTX程序不可执行 | GTX program is not executable: {config.gtx_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    # 检查GTX是否支持joint功能 🤝 | Check if GTX supports joint functionality
    if config.enable_joint_calling:
        try:
            result = subprocess.run(
                [config.gtx_path, 'joint', '--help'], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            if result.returncode == 0:
                logger.info("✅ ✓ GTX Joint Calling功能可用 | GTX Joint Calling functionality available")
            else:
                logger.warning("⚠️ GTX Joint Calling功能可能不可用，将在运行时检查 | GTX Joint Calling may not be available, will check at runtime")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.warning(f"⚠️ 无法验证GTX Joint功能：{e} | Cannot verify GTX Joint functionality: {e}")
    
    logger.info("✅ ✓ GTX程序检查通过 | GTX program check passed")
    return True

def cleanup_temp_files(config, logger):
    """清理临时文件 🧹 | Clean up temporary files"""
    logger.info("🧹 清理临时文件 | Cleaning up temporary files")
    
    try:
        temp_patterns = [
            "*.tmp",
            "*.temp", 
            "gtx_*",
            "*.log.bak"
        ]
        
        cleaned_count = 0
        for pattern in temp_patterns:
            temp_files = list(config.tmp_path.glob(pattern))
            for temp_file in temp_files:
                try:
                    if temp_file.is_file():
                        temp_file.unlink()
                        cleaned_count += 1
                    elif temp_file.is_dir():
                        import shutil
                        shutil.rmtree(temp_file)
                        cleaned_count += 1
                except Exception as e:
                    logger.debug(f"清理临时文件失败 | Failed to clean temp file {temp_file}: {e}")
        
        if cleaned_count > 0:
            logger.info(f"✅ 清理了 {cleaned_count} 个临时文件 | Cleaned {cleaned_count} temporary files")
        else:
            logger.info("ℹ️ 没有需要清理的临时文件 | No temporary files to clean")
            
    except Exception as e:
        logger.warning(f"⚠️ 清理临时文件过程中出现错误 | Error during temp file cleanup: {e}")