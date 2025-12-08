# # """
# # 🛠️ parabricks WGS分析工具函数模块 | parabricks WGS Analysis Utility Functions Module 🛠️
# # """

# # import logging
# # import subprocess
# # import sys
# # import os
# # from pathlib import Path
# # from datetime import datetime
# # import shutil

# # class parabricksLogger:
# #     """📝 parabricks分析日志管理器 | parabricks Analysis Logger Manager"""
    
# #     def __init__(self, output_dir: Path, log_name: str = "parabricks_analysis.log"):
# #         self.output_dir = output_dir
# #         self.log_file = output_dir / log_name
# #         self.setup_logging()
    
# #     def setup_logging(self):
# #         """设置日志 ✏️ | Setup logging"""
# #         if self.log_file.exists():
# #             # 备份现有日志文件 💾 | Backup existing log file
# #             backup_name = f"parabricks_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
# #             self.log_file.rename(self.output_dir / backup_name)
        
# #         logging.basicConfig(
# #             level=logging.INFO,
# #             format='[%(asctime)s] %(levelname)s - %(message)s',
# #             datefmt='%Y-%m-%d %H:%M:%S',
# #             handlers=[
# #                 logging.FileHandler(self.log_file),
# #                 logging.StreamHandler(sys.stdout)
# #             ]
# #         )
# #         self.logger = logging.getLogger(__name__)
    
# #     def get_logger(self):
# #         """获取日志器 📢 | Get logger"""
# #         return self.logger

# # class CommandRunner:
# #     """🚀 命令执行器 | Command Runner"""
    
# #     def __init__(self, logger, working_dir: Path):
# #         self.logger = logger
# #         self.working_dir = working_dir.resolve()
# #         self.container_cmd = None
# #         self.container_path = None
    
# #     def setup_container(self, container_path: Path):
# #         """设置容器环境 🐳 | Setup container environment"""
# #         # 检测使用 apptainer 还是 singularity
# #         for cmd in ['apptainer', 'singularity']:
# #             if subprocess.run(['which', cmd], capture_output=True).returncode == 0:
# #                 self.container_cmd = cmd
# #                 break
        
# #         if not self.container_cmd:
# #             raise RuntimeError("❌ 未找到Apptainer或Singularity | Apptainer or Singularity not found")
        
# #         self.logger.info(f"🐳 使用容器引擎 | Using container engine: {self.container_cmd}")
# #         self.container_path = container_path
    
# #     def run_container(self, cmd: list, description: str = "") -> bool:
# #         """执行容器命令 🐳 | Execute container command"""
# #         if not self.container_cmd or not self.container_path:
# #             raise RuntimeError("❌ 容器环境未初始化 | Container not initialized")
        
# #         if description:
# #             self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
# #         # 收集需要绑定的路径
# #         paths_to_bind = set()
# #         for item in cmd:
# #             item_str = str(item)
# #             if os.path.exists(item_str):
# #                 paths_to_bind.add(Path(item_str).parent.resolve())
        
# #         # 添加工作目录
# #         paths_to_bind.add(self.working_dir)
        
# #         # 构建绑定参数
# #         bind_args = [f"--bind={path}:{path}" for path in sorted(paths_to_bind)]
        
# #         # 构建完整的容器执行命令
# #         container_cmd = [
# #             self.container_cmd, "exec", "--nv",
# #             *bind_args,
# #             str(self.container_path)
# #         ] + [str(c) for c in cmd]
        
# #         self.logger.info(f"💻 容器命令 | Container command: {' '.join(container_cmd)}")
# #         self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
# #         try:
# #             result = subprocess.run(
# #                 container_cmd,
# #                 capture_output=True,
# #                 text=True,
# #                 check=True,
# #                 cwd=self.working_dir
# #             )
            
# #             self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
# #             if result.stdout:
# #                 self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
# #             return True
            
# #         except subprocess.CalledProcessError as e:
# #             self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
# #             self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
# #             self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
# #             self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
# #             return False
    
# #     def run(self, cmd: str, description: str = "") -> bool:
# #         """执行命令 ▶️ | Execute command"""
# #         if description:
# #             self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
# #         self.logger.info(f"💻 命令 | Command: {cmd}")
# #         self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
# #         try:
# #             result = subprocess.run(
# #                 cmd, 
# #                 shell=True, 
# #                 capture_output=True, 
# #                 text=True, 
# #                 check=True,
# #                 cwd=self.working_dir
# #             )
            
# #             self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
# #             if result.stdout:
# #                 self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
# #             return True
            
# #         except subprocess.CalledProcessError as e:
# #             self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
# #             self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
# #             self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
# #             self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
# #             return False

# # class FileProcessor:
# #     """📁 文件处理器 | File Processor"""
    
# #     def __init__(self, config, logger):
# #         self.config = config
# #         self.logger = logger
    
# #     def find_fastq_files(self):
# #         """查找FASTQ文件对 🔍 | Find FASTQ file pairs"""
# #         self.logger.info("🔍 搜索输入文件 | Searching input files")
        
# #         input_path = Path(self.config.input_dir)
# #         r1_files = list(input_path.glob(self.config.read1_pattern))
        
# #         if not r1_files:
# #             raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到 {self.config.read1_pattern} 文件")
        
# #         # 按文件名排序 🔠 | Sort by filename
# #         r1_files.sort()
        
# #         # 验证R2文件存在 ✅ | Validate R2 files exist
# #         file_pairs = []
# #         for r1_file in r1_files:
# #             # 提取样品名 🏷️ | Extract sample name
# #             sample_name = r1_file.name.replace("_1.clean.fq.gz", "")
            
# #             # 构建R2文件路径 🗺️ | Build R2 file path
# #             r2_file = input_path / f"{sample_name}_2.clean.fq.gz"
            
# #             if not r2_file.exists():
# #                 self.logger.warning(f"⚠️ 找不到对应的R2文件 | Cannot find R2 file: {r2_file}")
# #                 continue
            
# #             file_pairs.append((sample_name, str(r1_file), str(r2_file)))
        
# #         total_samples = len(file_pairs)
# #         self.logger.info(f"✅ 找到 {total_samples} 个样品 | Found {total_samples} samples")
        
# #         return file_pairs
    
# #     def check_output_exists(self, sample_name: str) -> bool:
# #         """检查输出文件是否已存在 🧐 | Check if output files already exist"""
# #         vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
# #         bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
# #         return vcf_file.exists() and bam_file.exists()
    
# #     def get_file_size(self, file_path: str) -> str:
# #         """获取文件大小 📏 | Get file size"""
# #         try:
# #             size_bytes = os.path.getsize(file_path)
# #             # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
# #             for unit in ['B', 'KB', 'MB', 'GB']:
# #                 if size_bytes < 1024.0:
# #                     return f"{size_bytes:.1f} {unit}"
# #                 size_bytes /= 1024.0
# #             return f"{size_bytes:.1f} TB"
# #         except:
# #             return "Unknown"

# # def check_dependencies(config, logger):
# #     """检查依赖软件 🧩 | Check dependencies"""
# #     logger.info("🧩 检查依赖软件 | Checking dependencies")
    
# #     # 检查parabricks程序 💻 | Check parabricks program
# #     if not os.path.exists(config.parabricks_path):
# #         error_msg = f"❌ parabricks程序不存在 | parabricks program does not exist: {config.parabricks_path}"
# #         logger.error(error_msg)
# #         raise RuntimeError(error_msg)
    
# #     # 检查parabricks程序是否可执行 🏃 | Check if parabricks program is executable
# #     if not os.access(config.parabricks_path, os.X_OK):
# #         error_msg = f"❌ parabricks程序不可执行 | parabricks program is not executable: {config.parabricks_path}"
# #         logger.error(error_msg)
# #         raise RuntimeError(error_msg)
    
# #     logger.info("✅ ✓ parabricks程序检查通过 | parabricks program check passed")

# #     # 检查bcftools（用于genotypegvcf步骤）
# #     # 如果workflow包含genotypegvcf或all，且启用了joint_calling
# #     if config.workflow in ["genotypegvcf", "all"] and config.joint_calling:
# #         logger.info("  |-- 🔗 检查bcftools (用于genotypegvcf) | Checking bcftools (for genotypegvcf)")
# #         if not shutil.which("bcftools"):
# #             error_msg = ("❌ bcftools未找到 | bcftools not found.\n"
# #                          "   它是genotypegvcf步骤创建VCF索引所必需的 | It is required for VCF indexing in genotypegvcf.\n"
# #                          "   请安装bcftools | Please install bcftools: `conda install -c bioconda bcftools`")
# #             logger.error(error_msg)
# #             raise RuntimeError(error_msg)
# #         logger.info("  ✅ ✓ bcftools检查通过 | bcftools check passed")

# #     return True

# """
# 🛠️ parabricks WGS分析工具函数模块 | parabricks WGS Analysis Utility Functions Module 🛠️
# """

# import logging
# import subprocess
# import sys
# import os
# from pathlib import Path
# from datetime import datetime
# import shutil
# import re # <-- 新增导入
# from typing import Optional # <-- 建议性新增导入

# class parabricksLogger:
#     # ... (此部分代码无需修改) ...
#     """📝 parabricks分析日志管理器 | parabricks Analysis Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "parabricks_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 ✏️ | Setup logging"""
#         if self.log_file.exists():
#             # 备份现有日志文件 💾 | Backup existing log file
#             backup_name = f"parabricks_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
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
#     # ... (此部分代码无需修改) ...
#     """🚀 命令执行器 | Command Runner"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir.resolve()
#         self.container_cmd = None
#         self.container_path = None
    
#     def setup_container(self, container_path: Path):
#         """设置容器环境 🐳 | Setup container environment"""
#         # 检测使用 apptainer 还是 singularity
#         for cmd in ['apptainer', 'singularity']:
#             if subprocess.run(['which', cmd], capture_output=True).returncode == 0:
#                 self.container_cmd = cmd
#                 break
        
#         if not self.container_cmd:
#             raise RuntimeError("❌ 未找到Apptainer或Singularity | Apptainer or Singularity not found")
        
#         self.logger.info(f"🐳 使用容器引擎 | Using container engine: {self.container_cmd}")
#         self.container_path = container_path
    
#     def run_container(self, cmd: list, description: str = "") -> bool:
#         """执行容器命令 🐳 | Execute container command"""
#         if not self.container_cmd or not self.container_path:
#             raise RuntimeError("❌ 容器环境未初始化 | Container not initialized")
        
#         if description:
#             self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
#         # 收集需要绑定的路径
#         paths_to_bind = set()
#         for item in cmd:
#             item_str = str(item)
#             if os.path.exists(item_str):
#                 paths_to_bind.add(Path(item_str).parent.resolve())
        
#         # 添加工作目录
#         paths_to_bind.add(self.working_dir)
        
#         # 构建绑定参数
#         bind_args = [f"--bind={path}:{path}" for path in sorted(paths_to_bind)]
        
#         # 构建完整的容器执行命令
#         container_cmd = [
#             self.container_cmd, "exec", "--nv",
#             *bind_args,
#             str(self.container_path)
#         ] + [str(c) for c in cmd]
        
#         self.logger.info(f"💻 容器命令 | Container command: {' '.join(container_cmd)}")
#         self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
#         try:
#             result = subprocess.run(
#                 container_cmd,
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

# # ================================================================= #
# # vvvvvvvvvvvvvv      这里是修改的核心区域      vvvvvvvvvvvvvv #
# # ================================================================= #

# class FileProcessor:
#     """📁 文件处理器 (优化版) | File Processor (Optimized)"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger

#     def _extract_sample_name(self, filename: str) -> str:
#         """从文件名中提取样本名 (通用)"""
#         name = filename
#         if name.endswith('.gz'):
#             name = name[:-3]
#         if name.endswith('.fq'):
#             name = name[:-3]
#         elif name.endswith('.fastq'):
#             name = name[:-6]
        
#         # 移除 R1/R2 或 1/2 标识符
#         name = re.sub(r'[_.]R[12]$', '', name)
#         name = re.sub(r'[_.]R[12](?=[_.]|$)', '', name, flags=re.IGNORECASE)
#         name = re.sub(r'[_.]1$|[_.]2$', '', name)
        
#         # 移除常见的中间词
#         name = name.replace('.clean', '').replace('_clean', '')
        
#         return name

#     def _find_paired_file(self, r1_file: Path) -> Optional[Path]:
#         """根据R1文件查找对应的R2文件 (通用)"""
#         r1_name = r1_file.name
#         possible_r2_names = set()

#         # 尝试各种替换规则
#         replacements = [
#             ('_1', '_2'), ('.1', '.2'), 
#             ('_R1', '_R2'), ('.R1', '.R2'),
#             ('_r1', '_r2'), ('.r1', '.r2')
#         ]
        
#         for r1_tag, r2_tag in replacements:
#             if r1_tag in r1_name:
#                 # 只替换第一次出现，避免修改样本名中的部分
#                 possible_r2_names.add(r1_name.replace(r1_tag, r2_tag, 1))

#         for r2_name in possible_r2_names:
#             if r2_name == r1_name: continue
#             r2_file = r1_file.with_name(r2_name)
#             if r2_file.exists():
#                 return r2_file
        
#         return None
    
#     def find_fastq_files(self):
#         """查找FASTQ文件对 (更健壮的版本) 🔍 | Find FASTQ file pairs (more robust version)"""
#         self.logger.info("🔍 搜索输入文件 | Searching input files")
        
#         input_path = Path(self.config.input_dir)
        
#         # 1. 查找所有可能的FASTQ文件
#         all_fastq_files = list(input_path.glob('*.fq')) + \
#                           list(input_path.glob('*.fastq')) + \
#                           list(input_path.glob('*.fq.gz')) + \
#                           list(input_path.glob('*.fastq.gz'))

#         if not all_fastq_files:
#             raise FileNotFoundError(f"在 {input_path} 中未找到任何FASTQ文件 (.fq, .fastq, .fq.gz, .fastq.gz)")

#         # 2. 识别R1文件
#         r1_pattern = re.compile(r'.*(_R1|_1)[\._].*', re.IGNORECASE)
#         r1_files = sorted([f for f in all_fastq_files if r1_pattern.match(f.name)])

#         if not r1_files:
#             raise FileNotFoundError(f"在 {input_path} 中未找到符合R1命名规则的文件 (例如 sample_1.fq.gz 或 sample_R1.fq)")

#         # 3. 配对并收集结果
#         file_pairs = []
#         processed_r2s = set() # 防止将R2文件也当作一个样本的R1
        
#         for r1_file in r1_files:
#             if r1_file in processed_r2s:
#                 continue
                
#             r2_file = self._find_paired_file(r1_file)
            
#             if not r2_file:
#                 self.logger.warning(f"⚠️ 找不到对应的R2文件: {r1_file.name}，将跳过此样本。")
#                 continue
            
#             sample_name = self._extract_sample_name(r1_file.name)
            
#             file_pairs.append((sample_name, str(r1_file), str(r2_file)))
#             processed_r2s.add(r2_file)
        
#         if not file_pairs:
#              raise FileNotFoundError(f"在 {input_path} 中成功找到R1文件，但未能完成任何R1/R2配对。请检查R2文件是否存在且命名是否规范。")

#         total_samples = len(file_pairs)
#         self.logger.info(f"✅ 找到 {total_samples} 个样品 | Found {total_samples} samples")
        
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

# # ================================================================= #
# # ^^^^^^^^^^^^^^      这里是修改的核心区域      ^^^^^^^^^^^^^^ #
# # ================================================================= #


# def check_dependencies(config, logger):
#     # ... (此部分代码无需修改) ...
#     """检查依赖软件 🧩 | Check dependencies"""
#     logger.info("🧩 检查依赖软件 | Checking dependencies")
    
#     # 检查parabricks程序 💻 | Check parabricks program
#     if not os.path.exists(config.parabricks_path):
#         error_msg = f"❌ parabricks程序不存在 | parabricks program does not exist: {config.parabricks_path}"
#         logger.error(error_msg)
#         raise RuntimeError(error_msg)
    
#     # 检查parabricks程序是否可执行 🏃 | Check if parabricks program is executable
#     if not os.access(config.parabricks_path, os.X_OK):
#         error_msg = f"❌ parabricks程序不可执行 | parabricks program is not executable: {config.parabricks_path}"
#         logger.error(error_msg)
#         raise RuntimeError(error_msg)
    
#     logger.info("✅ ✓ parabricks程序检查通过 | parabricks program check passed")

#     # 检查bcftools（用于genotypegvcf步骤）
#     # 如果workflow包含genotypegvcf或all，且启用了joint_calling
#     if config.workflow in ["genotypegvcf", "all"] and config.joint_calling:
#         logger.info("  |-- 🔗 检查bcftools (用于genotypegvcf) | Checking bcftools (for genotypegvcf)")
#         if not shutil.which("bcftools"):
#             error_msg = ("❌ bcftools未找到 | bcftools not found.\n"
#                          "   它是genotypegvcf步骤创建VCF索引所必需的 | It is required for VCF indexing in genotypegvcf.\n"
#                          "   请安装bcftools | Please install bcftools: `conda install -c bioconda bcftools`")
#             logger.error(error_msg)
#             raise RuntimeError(error_msg)
#         logger.info("  ✅ ✓ bcftools检查通过 | bcftools check passed")

#     return True

"""
🛠️ parabricks WGS分析工具函数模块 | parabricks WGS Analysis Utility Functions Module 🛠️
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime
import shutil

class parabricksLogger:
    """📝 parabricks分析日志管理器 | parabricks Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "parabricks_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 ✍️ | Setup logging"""
        if self.log_file.exists():
            backup_name = f"parabricks_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
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
        self.container_cmd = None
        self.container_path = None
    
    def setup_container(self, container_path: Path):
        """设置容器环境 🐳 | Setup container environment"""
        for cmd in ['apptainer', 'singularity']:
            if shutil.which(cmd):
                self.container_cmd = cmd
                break
        
        if not self.container_cmd:
            raise RuntimeError("❌ 未找到Apptainer或Singularity | Apptainer or Singularity not found")
        
        self.logger.info(f"🐳 使用容器引擎 | Using container engine: {self.container_cmd}")
        self.container_path = container_path
    
    def run_container(self, cmd: list, description: str = "") -> bool:
        """执行容器命令 🐳 | Execute container command"""
        if not self.container_cmd or not self.container_path:
            raise RuntimeError("❌ 容器环境未初始化，请先调用 setup_container() | Container not initialized")
        
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        paths_to_bind = set()
        for item in cmd:
            item_str = str(item)
            if os.path.exists(item_str) and not os.path.isdir(item_str):
                 paths_to_bind.add(Path(item_str).parent.resolve())
            elif os.path.isdir(item_str):
                 paths_to_bind.add(Path(item_str).resolve())

        paths_to_bind.add(self.working_dir)
        
        bind_args = [f"--bind={path}:{path}" for path in sorted(paths_to_bind)]
        
        container_cmd = [
            self.container_cmd, "exec", "--nv",
            *bind_args,
            str(self.container_path)
        ] + [str(c) for c in cmd]
        
        self.logger.info(f"💻 容器命令 | Container command: {' '.join(container_cmd)}")
        self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                container_cmd, capture_output=True, text=True, check=True, cwd=self.working_dir
            )
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            if result.stdout: self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
            return False
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行常规命令 ▶️ | Execute regular command"""
        if description: self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        try:
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, check=True, cwd=self.working_dir
            )
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            if result.stdout: self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
            return False

class FileProcessor:
    """📁 文件处理器 | File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def find_fastq_files(self):
        """查找FASTQ文件对 🔍 (更智能的版本) | Find FASTQ file pairs (smarter version)"""
        self.logger.info(f"🔍 搜索输入文件，模式: {self.config.read1_pattern} | Searching input files with pattern: {self.config.read1_pattern}")

        input_path = Path(self.config.input_dir)
        
        # --- MODIFIED: 动态样本名提取逻辑 ---
        if '*' not in self.config.read1_pattern or '*' not in self.config.read2_pattern:
            raise ValueError("❌ read1_pattern 和 read2_pattern 都必须包含 '*' 通配符 | read1_pattern and read2_pattern must contain a '*' wildcard")

        r1_prefix, r1_suffix = self.config.read1_pattern.split('*', 1)
        r2_prefix, r2_suffix = self.config.read2_pattern.split('*', 1)

        r1_files = sorted(list(input_path.glob(self.config.read1_pattern)))
        
        if not r1_files:
            raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到匹配 {self.config.read1_pattern} 的文件")
        
        file_pairs = []
        for r1_file in r1_files:
            sample_name = r1_file.name
            if r1_prefix: sample_name = sample_name.removeprefix(r1_prefix)
            if r1_suffix: sample_name = sample_name.removesuffix(r1_suffix)
            
            r2_file_name = f"{r2_prefix}{sample_name}{r2_suffix}"
            r2_file = input_path / r2_file_name
            
            if not r2_file.exists():
                self.logger.warning(f"⚠️ 找不到样品 '{sample_name}' 对应的R2文件: {r2_file.name} | Cannot find corresponding R2 file for sample '{sample_name}': {r2_file.name}")
                continue
            
            file_pairs.append((sample_name, str(r1_file), str(r2_file)))
        
        self.logger.info(f"✅ 找到 {len(file_pairs)} 个样品需要处理 | Found {len(file_pairs)} samples to process")
        return file_pairs
    
    def check_output_exists(self, sample_name: str) -> bool:
        """检查输出文件是否已存在 🧐 | Check if output files already exist"""
        # --- MODIFIED: 根据配置决定检查VCF还是GVCF ---
        if self.config.gvcf:
            vcf_file = self.config.vcf_output_dir / f"{sample_name}.g.vcf.gz"
        else:
            vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
            
        bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        return vcf_file.exists() and bam_file.exists()
    
    def get_file_size(self, file_path: str) -> str:
        """获取文件大小 📏 | Get file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
                if size_bytes < 1024.0:
                    return f"{size_bytes:.1f} {unit}"
                size_bytes /= 1024.0
            return f"{size_bytes:.1f} PB"
        except FileNotFoundError:
            return "Not Found"
        except Exception:
            return "Unknown"

def check_dependencies(config, logger):
    """检查依赖软件 🧩 | Check dependencies"""
    logger.info("🧩 检查依赖软件 | Checking dependencies")
    
    if not os.path.exists(config.parabricks_path):
        error_msg = f"❌ parabricks程序不存在 | parabricks program does not exist: {config.parabricks_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    if not os.access(config.parabricks_path, os.X_OK):
        error_msg = f"❌ parabricks程序不可执行 | parabricks program is not executable: {config.parabricks_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    logger.info("✅ ✓ parabricks程序检查通过 | parabricks program check passed")

    if config.joint_calling:
        logger.info("  |-- 🔗 Joint Calling已启用，检查bcftools | Joint Calling enabled, checking for bcftools")
        if not shutil.which("bcftools"):
            error_msg = ("❌ bcftools未找到或不在系统PATH中 | bcftools not found or not in system PATH.\n"
                         "   它是在Joint Calling后为VCF.gz文件建立索引所必需的 | It is required for indexing the VCF.gz file after Joint Calling.\n"
                         "   请安装bcftools (例如 `conda install -c bioconda bcftools`) | Please install bcftools (e.g., `conda install -c bioconda bcftools`)")
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        logger.info("  ✅ ✓ bcftools检查通过 | bcftools check passed")

    return True