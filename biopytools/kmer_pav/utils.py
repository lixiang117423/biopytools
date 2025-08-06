"""
🛠️ K-mer PAV分析工具函数模块 | K-mer PAV Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
import glob
import re
from pathlib import Path

class PAVLogger:
    """📝 PAV分析日志管理器 | PAV Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "kmer_pav_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """📋 设置日志 | Setup logging"""
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
        """📜 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """🏃 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()  # 使用绝对路径 | Use absolute path
    
    def run(self, cmd: str, description: str = "", log_file: str = None) -> bool:
        """▶️ 执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        try:
            # 如果指定了日志文件，则重定向输出 | Redirect output if log file specified
            if log_file:
                full_cmd = f"{cmd} 2>&1 | tee {log_file}"
            else:
                full_cmd = cmd
            
            result = subprocess.run(
                full_cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📄 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📄 标准输出 | Stdout: {e.stdout}")
            return False

class SampleDiscovery:
    """🔍 样本发现器 | Sample Discovery"""
    
    def __init__(self, fastq_dir: str, fastq_pattern: str, logger):
        self.fastq_dir = fastq_dir
        self.fastq_pattern = fastq_pattern
        self.logger = logger
    
    def discover_samples(self) -> list:
        """🔎 发现fastq样本 | Discover fastq samples"""
        self.logger.info(f"🔍 从目录中发现样本 | Discovering samples from directory: {self.fastq_dir}")
        self.logger.info(f"🎯 使用模式 | Using pattern: {self.fastq_pattern}")
        
        # 首先检查目录是否存在并可读 | First check if directory exists and is readable
        if not os.path.exists(self.fastq_dir):
            raise ValueError(f"❌ FASTQ目录不存在 | FASTQ directory does not exist: {self.fastq_dir}")
        
        if not os.access(self.fastq_dir, os.R_OK):
            raise ValueError(f"❌ FASTQ目录不可读 | FASTQ directory is not readable: {self.fastq_dir}")
        
        # 使用用户指定的模式查找文件 | Find files using user-specified pattern
        pattern = os.path.join(self.fastq_dir, self.fastq_pattern)
        self.logger.info(f"🔍 完整搜索模式 | Full search pattern: {pattern}")
        
        r1_files = glob.glob(pattern)
        self.logger.info(f"📊 使用模式 '{self.fastq_pattern}' 找到 {len(r1_files)} 个文件 | Found {len(r1_files)} files with pattern '{self.fastq_pattern}'")
        
        # 显示找到的前几个文件用于调试 | Show first few files found for debugging
        if r1_files:
            self.logger.info("📁 找到的文件 | Found files:")
            for i, f in enumerate(r1_files[:5]):
                self.logger.info(f"  {i+1}. {os.path.basename(f)}")
            if len(r1_files) > 5:
                self.logger.info(f"  ... 共{len(r1_files)}个文件 | ... total {len(r1_files)} files")
        
        if not r1_files:
            # 如果用户指定的模式没有找到文件，尝试常见模式 | If user pattern finds no files, try common patterns
            self.logger.warning(f"⚠️ 使用模式 {self.fastq_pattern} 未找到文件，尝试常见模式 | No files found with pattern {self.fastq_pattern}, trying common patterns")
            
            # 先列出目录中所有文件 | First list all files in directory
            try:
                all_files = os.listdir(self.fastq_dir)
                fq_files = [f for f in all_files if f.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fastq'))]
                self.logger.info(f"📂 目录中发现 {len(fq_files)} 个FASTQ相关文件 | Found {len(fq_files)} FASTQ-related files in directory")
                if fq_files:
                    self.logger.info("📋 前几个FASTQ文件 | First few FASTQ files:")
                    for i, f in enumerate(fq_files[:5]):
                        self.logger.info(f"  {i+1}. {f}")
            except Exception as e:
                self.logger.error(f"❌ 无法列出目录内容: {e} | Cannot list directory contents: {e}")
            
            patterns = ["*_1.fq.gz", "*_R1.fastq.gz", "*_R1.fq.gz", "*_1.fastq.gz"]
            for fallback_pattern in patterns:
                full_pattern = os.path.join(self.fastq_dir, fallback_pattern)
                self.logger.info(f"🔄 尝试模式 | Trying pattern: {fallback_pattern}")
                r1_files = glob.glob(full_pattern)
                if r1_files:
                    self.logger.info(f"✅ 使用备用模式找到 {len(r1_files)} 个文件 | Found {len(r1_files)} files with fallback pattern: {fallback_pattern}")
                    self.fastq_pattern = fallback_pattern
                    break
        
        if not r1_files:
            raise ValueError(f"❌ 在目录中未找到fastq文件 | No fastq files found in directory: {self.fastq_dir}")
        
        # 从模式中提取样本名称 | Extract sample names from pattern
        samples = []
        failed_extractions = []
        
        self.logger.info("🔧 开始提取样本名称 | Starting sample name extraction")
        
        # 显示模式转换信息 | Show pattern conversion info
        pattern_base = os.path.basename(self.fastq_pattern)
        # regex_pattern = pattern_base.replace('*', '(.+)').replace('.', r'\.')
        # 先转义所有的点号，然后再替换星号
        regex_pattern = pattern_base.replace('.', r'\.').replace('*', '(.+)')
        self.logger.info(f"🔄 模式转换 | Pattern conversion: '{pattern_base}' -> '{regex_pattern}'")
        
        for i, r1_file in enumerate(r1_files):
            basename = os.path.basename(r1_file)
            # 从模式中获取样本名称 | Get sample name from pattern
            sample_name = self._extract_sample_name(basename, self.fastq_pattern)
            if sample_name:
                samples.append(sample_name)
                if i < 3:  # 只显示前3个样本的详细信息 | Only show details for first 3 samples
                    self.logger.info(f"✅ 从文件 {basename} 提取样本名: {sample_name} | Extracted sample name from {basename}: {sample_name}")
            else:
                failed_extractions.append(basename)
                self.logger.warning(f"⚠️ 无法从文件 {basename} 提取样本名 | Cannot extract sample name from {basename}")
        
        if failed_extractions:
            self.logger.warning(f"⚠️ 以下 {len(failed_extractions)} 个文件无法提取样本名 | Cannot extract sample names from {len(failed_extractions)} files:")
            for f in failed_extractions[:5]:
                self.logger.warning(f"  - {f}")
        
        samples = sorted(list(set(samples)))  # 去重并排序 | Remove duplicates and sort
        
        self.logger.info(f"🎯 发现 {len(samples)} 个样本 | Discovered {len(samples)} samples")
        if len(samples) <= 10:
            for sample in samples:
                self.logger.info(f"  📋 {sample}")
        else:
            for sample in samples[:10]:
                self.logger.info(f"  📋 {sample}")
            self.logger.info(f"... (显示前10个样本) | ... (showing first 10 samples)")
        
        return samples
    
    def _extract_sample_name(self, filename: str, pattern: str) -> str:
        """🔧 从文件名中根据模式提取样本名称 | Extract sample name from filename based on pattern"""
        # 将模式转换为正则表达式 | Convert pattern to regex
        # 例如 "*_1.fq.gz" -> "(.+)_1\.fq\.gz"
        pattern_base = os.path.basename(pattern)
        # 先转义所有的点号，然后再替换星号
        regex_pattern = pattern_base.replace('.', r'\.').replace('*', '(.+)')
        
        # 添加详细的调试信息（只对前3个文件）
        self.logger.info(f"🔍 调试匹配 | Debug matching:")
        self.logger.info(f"  📄 文件名 | Filename: '{filename}'")
        self.logger.info(f"  🎯 正则模式 | Regex pattern: '{regex_pattern}'")
        
        try:
            match = re.match(regex_pattern, filename)
            if match:
                sample_name = match.group(1)
                self.logger.info(f"  ✅ 匹配成功 | Match successful: '{sample_name}'")
                return sample_name
            else:
                self.logger.warning(f"  ❌ 正则匹配失败 | Regex match failed")
                
                # 让我们尝试手动测试一个简单的例子
                test_match = re.match(r'(.+)_1\.fq\.gz', filename)
                if test_match:
                    self.logger.info(f"  🔧 手动测试成功 | Manual test successful: '{test_match.group(1)}'")
                    return test_match.group(1)  # 如果手动测试成功，就返回结果
                else:
                    self.logger.warning(f"  ❌ 手动测试也失败 | Manual test also failed")
                
                return None
        except Exception as e:
            self.logger.warning(f"❌ 正则匹配出错: {e} | Regex matching error: {e}")
            return None
    
    def get_sample_files(self, sample_name: str) -> tuple:
        """📂 获取样本的双端文件路径 | Get paired-end file paths for sample"""
        # 根据第一个文件的模式推断第二个文件的模式 | Infer second file pattern from first file pattern
        pattern_base = os.path.basename(self.fastq_pattern)
        
        # 构建可能的双端文件模式 | Build possible paired-end file patterns
        possible_patterns = []
        
        if "_1." in pattern_base:
            # 如果模式包含_1.，那么对应的应该是_2.
            r2_pattern_base = pattern_base.replace("_1.", "_2.")
            r1_pattern = self.fastq_pattern.replace('*', sample_name)
            r2_pattern = self.fastq_pattern.replace(pattern_base, r2_pattern_base).replace('*', sample_name)
            possible_patterns.append((r1_pattern, r2_pattern))
        elif "_R1." in pattern_base:
            # 如果模式包含_R1.，那么对应的应该是_R2.
            r2_pattern_base = pattern_base.replace("_R1.", "_R2.")
            r1_pattern = self.fastq_pattern.replace('*', sample_name)
            r2_pattern = self.fastq_pattern.replace(pattern_base, r2_pattern_base).replace('*', sample_name)
            possible_patterns.append((r1_pattern, r2_pattern))
        
        # 添加通用的备用模式 | Add generic fallback patterns
        fallback_patterns = [
            (f"{sample_name}_1.fq.gz", f"{sample_name}_2.fq.gz"),
            (f"{sample_name}_R1.fastq.gz", f"{sample_name}_R2.fastq.gz"),
            (f"{sample_name}_R1.fq.gz", f"{sample_name}_R2.fq.gz"),
            (f"{sample_name}_1.fastq.gz", f"{sample_name}_2.fastq.gz")
        ]
        
        # 首先尝试从模式推断的配对文件 | First try inferred paired files from pattern
        for r1_pattern, r2_pattern in possible_patterns:
            self.logger.info(f"🔍 尝试配对文件 | Trying paired files: {os.path.basename(r1_pattern)} & {os.path.basename(r2_pattern)}")
            if os.path.exists(r1_pattern) and os.path.exists(r2_pattern):
                self.logger.info(f"✅ 找到配对文件 | Found paired files for {sample_name}")
                return r1_pattern, r2_pattern
        
        # 然后尝试备用模式 | Then try fallback patterns
        for r1_pattern, r2_pattern in fallback_patterns:
            r1_path = os.path.join(self.fastq_dir, r1_pattern)
            r2_path = os.path.join(self.fastq_dir, r2_pattern)
            
            if os.path.exists(r1_path) and os.path.exists(r2_path):
                self.logger.info(f"✅ 使用备用模式找到配对文件 | Found paired files using fallback pattern for {sample_name}")
                return r1_path, r2_path
        
        raise FileNotFoundError(f"❌ 找不到样本 {sample_name} 的双端fastq文件 | Could not find paired fastq files for sample {sample_name}")

def check_dependencies(config, logger):
    """🔍 检查依赖软件 | Check dependencies"""
    logger.info("🔧 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.unikmer_path, "unikmer")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✅ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True
