# """
# 🛠️ BLAST分析工具函数模块 | BLAST Analysis Utility Functions Module
# """

# import logging
# import subprocess
# import sys
# import re
# import os
# from pathlib import Path
# from datetime import datetime
# from typing import Dict, List, Optional

# class BLASTLogger:
#     """📝 BLAST分析日志管理器"""
    
#     def __init__(self, output_dir: Path, log_name: str = "blast_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志"""
#         if self.log_file.exists():
#             self.log_file.unlink()
        
#         logging.basicConfig(
#             level=logging.INFO,
#             format='%(asctime)s - %(levelname)s - %(message)s',
#             handlers=[
#                 logging.FileHandler(self.log_file, encoding='utf-8'),
#                 logging.StreamHandler(sys.stdout)
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def get_logger(self):
#         """获取日志器"""
#         return self.logger

# class CommandRunner:
#     """⚡ 命令执行器"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir.resolve()

#     def run(self, cmd: str, description: str = "") -> bool:
#         """执行命令"""
#         if description:
#             self.logger.info(f"执行步骤: {description}")
        
#         self.logger.info(f"命令: {cmd}")
        
#         try:
#             result = subprocess.run(
#                 cmd, shell=True, capture_output=True, text=True, 
#                 check=True  # 移除 cwd=self.working_dir
#             )
#             self.logger.info(f"命令执行成功: {description}")
#             return True
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"命令执行失败: {description}")
#             self.logger.error(f"错误代码: {e.returncode}")
#             self.logger.error(f"错误信息: {e.stderr}")
#             return False

# class SampleMapGenerator:
#     """🗂️ 样品映射文件生成器"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_sample_map(self) -> str:
#         """生成样品映射文件"""
#         if self.config.sample_map_file:
#             self.logger.info(f"使用用户提供的样品映射文件: {self.config.sample_map_file}")
#             return self.config.sample_map_file
        
#         if not self.config.input_path:
#             raise ValueError("既没有提供样品映射文件也没有提供输入路径")
        
#         self.logger.info("自动生成样品映射文件")
        
#         # 获取输入文件列表
#         input_files = get_input_files(self.config.input_path, self.config.input_suffix)
        
#         if not input_files:
#             error_msg = f"在路径中未找到匹配的输入文件: {self.config.input_path}"
#             self.logger.error(error_msg)
#             raise FileNotFoundError(error_msg)
        
#         # 生成样品映射文件路径
#         map_file_path = self.config.output_path / "auto_generated_sample_map.txt"
        
#         self.logger.info(f"生成样品映射文件: {map_file_path}")
        
#         # 写入样品映射文件
#         with open(map_file_path, 'w', encoding='utf-8') as f:
#             f.write("# 自动生成的样品映射文件 | Auto-generated sample mapping file\n")
#             f.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
#             f.write(f"# 输入路径: {self.config.input_path}\n")
#             f.write("# 格式: 文件路径<TAB>样品名称\n")
#             f.write("#" + "="*60 + "\n")
            
#             sample_count = 0
#             for input_file in sorted(input_files):
#                 sample_name = self._extract_sample_name_from_path(input_file)
#                 f.write(f"{input_file}\t{sample_name}\n")
#                 sample_count += 1
        
#         self.logger.info(f"样品映射文件生成完成，包含 {sample_count} 个样品")
#         self.logger.info(f"提示：如需修改样品名称，请编辑文件后重新运行")
        
#         # 更新配置
#         self.config.sample_map_file = str(map_file_path)
#         self.config.auto_generated_map = True
        
#         return str(map_file_path)
    
#     def _extract_sample_name_from_path(self, file_path: str) -> str:
#         """从文件路径提取样品名称"""
#         filename = os.path.basename(file_path)
        
#         if self.config.auto_detect_samples:
#             match = re.search(self.config.sample_name_pattern, filename)
#             if match:
#                 return match.group(1)
        
#         return Path(filename).stem

# class SampleManager:
#     """🧪 样品管理器"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
#         self.sample_mapping = {}
#         self.file_paths = []
    
#     def load_sample_mapping(self) -> Dict[str, str]:
#         """加载样品映射和文件路径"""
#         if not self.config.sample_map_file:
#             raise ValueError("样品映射文件路径未设置")
        
#         self.logger.info(f"从映射文件加载样本和文件路径: {self.config.sample_map_file}")
        
#         try:
#             with open(self.config.sample_map_file, 'r', encoding='utf-8') as f:
#                 for line_no, line in enumerate(f, 1):
#                     line = line.strip()
#                     if not line or line.startswith('#'):
#                         continue
                    
#                     parts = line.split('\t')
#                     if len(parts) >= 2:
#                         file_path = parts[0]
#                         sample_name = parts[1]
                        
#                         if not os.path.exists(file_path):
#                             self.logger.warning(f"文件不存在: {file_path}")
#                             continue
                        
#                         filename = os.path.basename(file_path)
#                         self.sample_mapping[filename] = sample_name
#                         self.file_paths.append(file_path)
#                     else:
#                         self.logger.warning(f"第{line_no}行格式不正确")
            
#             self.logger.info(f"成功加载{len(self.sample_mapping)}个样品映射")
            
#         except Exception as e:
#             self.logger.error(f"无法加载样品映射文件: {e}")
#             raise
            
#         return self.sample_mapping
    
#     def get_file_paths(self) -> List[str]:
#         """获取文件路径列表"""
#         return self.file_paths
    
#     def extract_sample_name(self, filename: str) -> str:
#         """提取样品名称"""
#         if self.sample_mapping and filename in self.sample_mapping:
#             return self.sample_mapping[filename]
#         return Path(filename).stem

# def check_dependencies(config, logger) -> bool:
#     """检查依赖软件"""
#     logger.info("🔍 检查依赖软件")
    
#     dependencies = [
#         (config.makeblastdb_path, "makeblastdb"),
#         (config.get_blast_executable(), config.blast_type.upper())
#     ]
    
#     missing_deps = []
    
#     for cmd, name in dependencies:
#         try:
#             result = subprocess.run([cmd, "-version"], 
#                                   capture_output=True, text=True, timeout=10)
#             if result.returncode == 0:
#                 logger.info(f"{name} 可用")
#             else:
#                 missing_deps.append(name)
#         except (subprocess.TimeoutExpired, FileNotFoundError):
#             missing_deps.append(name)
    
#     if missing_deps:
#         error_msg = f"缺少依赖软件: {', '.join(missing_deps)}"
#         logger.error(error_msg)
#         raise RuntimeError(error_msg)
    
#     return True

# def get_input_files(input_path: str, suffix: str) -> List[str]:
#     """获取输入文件列表"""
#     input_path_obj = Path(input_path)
    
#     # 如果是单个文件，直接返回
#     if input_path_obj.is_file():
#         return [str(input_path_obj)]
    
#     # 如果是目录，扫描文件
#     if input_path_obj.is_dir():
#         if suffix.startswith('*'):
#             pattern = suffix[1:]
#             files = list(input_path_obj.glob(f"*{pattern}"))
#         else:
#             files = list(input_path_obj.glob(suffix))
#         return [str(f) for f in files if f.is_file()]
    
#     return []

"""
🛠️ BLAST分析工具函数模块 | BLAST Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import re
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

class BLASTLogger:
    """📝 BLAST分析日志管理器 | BLAST Analysis Log Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "blast_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup Logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get Logger"""
        return self.logger

class CommandRunner:
    """⚡ 命令执行器 | Command Executor"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute Command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, 
                check=True  # 移除 cwd=self.working_dir | Remove cwd=self.working_dir
            )
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False

class SampleMapGenerator:
    """🗂️ 样品映射文件生成器 | Sample Map File Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_sample_map(self) -> str:
        """生成样品映射文件 | Generate Sample Map File"""
        if self.config.sample_map_file:
            self.logger.info(f"使用用户提供的样品映射文件 | Using user-provided sample map file: {self.config.sample_map_file}")
            return self.config.sample_map_file
        
        if not self.config.input_path:
            raise ValueError("既没有提供样品映射文件也没有提供输入路径 | Neither a sample map file nor an input path was provided")
        
        self.logger.info("自动生成样品映射文件 | Auto-generating sample map file")
        
        # 获取输入文件列表 | Get the list of input files
        input_files = get_input_files(self.config.input_path, self.config.input_suffix)
        
        if not input_files:
            error_msg = f"在路径中未找到匹配的输入文件 | No matching input files found in path: {self.config.input_path}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        # 生成样品映射文件路径 | Generate sample map file path
        map_file_path = self.config.output_path / "auto_generated_sample_map.txt"
        
        self.logger.info(f"生成样品映射文件 | Generating sample map file: {map_file_path}")
        
        # 写入样品映射文件 | Write to the sample map file
        with open(map_file_path, 'w', encoding='utf-8') as f:
            f.write("# 自动生成的样品映射文件 | Auto-generated sample mapping file\n")
            f.write(f"# 生成时间 | Generation Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 输入路径 | Input Path: {self.config.input_path}\n")
            f.write("# 格式: 文件路径<TAB>样品名称 | Format: file_path<TAB>sample_name\n")
            f.write("#" + "="*60 + "\n")
            
            sample_count = 0
            for input_file in sorted(input_files):
                sample_name = self._extract_sample_name_from_path(input_file)
                f.write(f"{input_file}\t{sample_name}\n")
                sample_count += 1
        
        self.logger.info(f"样品映射文件生成完成，包含 {sample_count} 个样品 | Sample map file generation complete, containing {sample_count} samples")
        self.logger.info(f"提示：如需修改样品名称，请编辑文件后重新运行 | Tip: To modify sample names, please edit the file and run again")
        
        # 更新配置 | Update configuration
        self.config.sample_map_file = str(map_file_path)
        self.config.auto_generated_map = True
        
        return str(map_file_path)
    
    def _extract_sample_name_from_path(self, file_path: str) -> str:
        """从文件路径提取样品名称 | Extract Sample Name from File Path"""
        filename = os.path.basename(file_path)
        
        if self.config.auto_detect_samples:
            match = re.search(self.config.sample_name_pattern, filename)
            if match:
                return match.group(1)
        
        return Path(filename).stem

class SampleManager:
    """🧪 样品管理器 | Sample Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.sample_mapping = {}
        self.file_paths = []
    
    def load_sample_mapping(self) -> Dict[str, str]:
        """加载样品映射和文件路径 | Load Sample Mapping and File Paths"""
        if not self.config.sample_map_file:
            raise ValueError("样品映射文件路径未设置 | Sample map file path is not set")
        
        self.logger.info(f"从映射文件加载样本和文件路径 | Loading samples and file paths from map file: {self.config.sample_map_file}")
        
        try:
            with open(self.config.sample_map_file, 'r', encoding='utf-8') as f:
                for line_no, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        file_path = parts[0]
                        sample_name = parts[1]
                        
                        if not os.path.exists(file_path):
                            self.logger.warning(f"文件不存在 | File does not exist: {file_path}")
                            continue
                        
                        filename = os.path.basename(file_path)
                        self.sample_mapping[filename] = sample_name
                        self.file_paths.append(file_path)
                    else:
                        self.logger.warning(f"第{line_no}行格式不正确 | Line {line_no} has incorrect format")
            
            self.logger.info(f"成功加载{len(self.sample_mapping)}个样品映射 | Successfully loaded {len(self.sample_mapping)} sample mappings")
            
        except Exception as e:
            self.logger.error(f"无法加载样品映射文件 | Failed to load sample map file: {e}")
            raise
            
        return self.sample_mapping
    
    def get_file_paths(self) -> List[str]:
        """获取文件路径列表 | Get File Path List"""
        return self.file_paths
    
    def extract_sample_name(self, filename: str) -> str:
        """提取样品名称 | Extract Sample Name"""
        if self.sample_mapping and filename in self.sample_mapping:
            return self.sample_mapping[filename]
        return Path(filename).stem

def check_dependencies(config, logger) -> bool:
    """检查依赖软件 | Check for Dependent Software"""
    logger.info("🔍 检查依赖软件 | Checking for dependent software")
    
    dependencies = [
        (config.makeblastdb_path, "makeblastdb"),
        (config.get_blast_executable(), config.blast_type.upper())
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "-version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"{name} 可用 | {name} is available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependent software: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True

def get_input_files(input_path: str, suffix: str) -> List[str]:
    """获取输入文件列表 | Get Input File List"""
    input_path_obj = Path(input_path)
    
    # 如果是单个文件，直接返回 | If it is a single file, return it directly
    if input_path_obj.is_file():
        return [str(input_path_obj)]
    
    # 如果是目录，扫描文件 | If it is a directory, scan for files
    if input_path_obj.is_dir():
        if suffix.startswith('*'):
            pattern = suffix[1:]
            files = list(input_path_obj.glob(f"*{pattern}"))
        else:
            files = list(input_path_obj.glob(suffix))
        return [str(f) for f in files if f.is_file()]
    
    return []