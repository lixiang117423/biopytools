# """
# 群体遗传分析工具函数模块 | Population Genetics Analysis Utility Functions Module
# """

# import logging
# import subprocess
# import sys
# from pathlib import Path

# class PopGenLogger:
#     """群体遗传分析日志管理器 | Population Genetics Analysis Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "popgen_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 | Setup logging"""
#         if self.log_file.exists():
#             self.log_file.unlink()
        
#         logging.basicConfig(
#             level=logging.INFO,
#             format='%(asctime)s - %(levelname)s - %(message)s',
#             handlers=[
#                 logging.FileHandler(self.log_file),
#                 logging.StreamHandler(sys.stdout)
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def get_logger(self):
#         """获取日志器 | Get logger"""
#         return self.logger

# class CommandRunner:
#     """命令执行器 | Command Runner"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir.resolve()
    
#     def run(self, cmd: str, description: str = "") -> bool:
#         """执行命令 | Execute command"""
#         if description:
#             self.logger.info(f"执行步骤 | Executing step: {description}")
        
#         self.logger.info(f"命令 | Command: {cmd}")
        
#         try:
#             result = subprocess.run(
#                 cmd, shell=True, capture_output=True, text=True, 
#                 check=True, cwd=self.working_dir
#             )
#             self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
#             return True
            
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"命令执行失败 | Command execution failed: {description}")
#             self.logger.error(f"错误信息 | Error message: {e.stderr}")
#             return False

# # def check_dependencies(config, logger):
# #     """检查依赖软件 | Check dependencies"""
# #     logger.info("检查依赖软件 | Checking dependencies")
    
# #     dependencies = [
# #         (config.vcftools_path, "VCFtools"),
# #         (config.plink_path, "PLINK"),
# #         (config.bcftools_path, "BCFtools")
# #     ]
    
# #     if config.calculate_ne:
# #         dependencies.append((config.smcpp_path, "SMC++"))
    
# #     missing_deps = []
    
# #     for cmd, name in dependencies:
# #         try:
# #             result = subprocess.run([cmd, "--version"], 
# #                                   capture_output=True, text=True, timeout=10)
# #             if result.returncode == 0:
# #                 logger.info(f"✓ {name} 可用 | available")
# #             else:
# #                 missing_deps.append(name)
# #         except (subprocess.TimeoutExpired, FileNotFoundError):
# #             missing_deps.append(name)
    
# #     if missing_deps:
# #         error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
# #         logger.error(error_msg)
# #         return False
    
# #     return True

# def check_dependencies(config, logger):
#     """检查依赖软件 | Check dependencies"""
#     logger.info("检查依赖软件 | Checking dependencies")
    
#     # 依赖软件配置：(路径, 名称, 版本命令参数)
#     dependencies = [
#         (config.vcftools_path, "VCFtools", "--version"),
#         (config.plink_path, "PLINK", "--version"),
#         (config.bcftools_path, "BCFtools", "--version")
#     ]
    
#     # 如果需要计算有效群体大小，添加SMC++依赖
#     if config.calculate_ne:
#         dependencies.append((config.smcpp_path, "SMC++", "version"))  # SMC++使用"version"而不是"--version"
    
#     missing_deps = []
    
#     for cmd, name, version_arg in dependencies:
#         try:
#             # 根据不同工具使用不同的版本检查命令
#             result = subprocess.run([cmd, version_arg], 
#                                   capture_output=True, text=True, timeout=10)
            
#             if result.returncode == 0:
#                 logger.info(f"✓ {name} 可用 | available")
                
#                 # 可选：显示版本信息（用于调试）
#                 if logger.level <= 10:  # DEBUG级别
#                     version_info = result.stdout.strip().split('\n')[0] if result.stdout else "版本信息不可用"
#                     logger.debug(f"  {name} 版本: {version_info}")
                    
#             else:
#                 logger.warning(f"✗ {name} 命令执行失败，返回码: {result.returncode}")
#                 if result.stderr:
#                     logger.debug(f"  错误信息: {result.stderr.strip()}")
#                 missing_deps.append(name)
                
#         except subprocess.TimeoutExpired:
#             logger.warning(f"✗ {name} 命令执行超时")
#             missing_deps.append(name)
            
#         except FileNotFoundError:
#             logger.warning(f"✗ {name} 未找到，请检查安装路径: {cmd}")
#             missing_deps.append(name)
            
#         except Exception as e:
#             logger.warning(f"✗ {name} 检查时发生意外错误: {e}")
#             missing_deps.append(name)
    
#     # 报告结果
#     if missing_deps:
#         error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
#         logger.error(error_msg)
        
#         # 提供安装建议
#         logger.error("安装建议 | Installation suggestions:")
#         for dep in missing_deps:
#             if dep == "VCFtools":
#                 logger.error("  VCFtools: conda install -c conda-forge vcftools")
#             elif dep == "PLINK":
#                 logger.error("  PLINK: conda install -c conda-forge plink")
#             elif dep == "BCFtools":
#                 logger.error("  BCFtools: conda install -c conda-forge bcftools")
#             elif dep == "SMC++":
#                 logger.error("  SMC++: pip install smcpp")
        
#         return False
    
#     logger.info("✓ 所有依赖软件检查通过 | All dependencies check passed")
#     return True

"""
群体遗传分析工具函数模块 | Population Genetics Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path

class PopGenLogger:
    """群体遗传分析日志管理器 | Population Genetics Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "popgen_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 确保输出目录存在 | Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 如果日志文件存在，删除旧的 | Remove old log file if exists
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
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, output_dir: Path):
        self.logger = logger
        self.output_dir = output_dir
    
    def run(self, cmd: str, description: str = "", timeout: int = 3600) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            # 改变到输出目录执行命令 | Change to output directory to execute command
            original_cwd = Path.cwd()
            self.output_dir.mkdir(parents=True, exist_ok=True)
            
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8',
                timeout=timeout,
                cwd=str(self.output_dir)
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            # 如果有输出，记录到debug级别 | Log output to debug level if exists
            if result.stdout and result.stdout.strip():
                self.logger.debug(f"标准输出 | Stdout: {result.stdout[:500]}...")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息 | Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired as e:
            self.logger.error(f"命令执行超时 | Command execution timeout: {description}")
            self.logger.error(f"超时时间 | Timeout: {timeout} seconds")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常 | Command execution exception: {description}")
            self.logger.error(f"异常信息 | Exception: {e}")
            return False
    
    def run_with_output(self, cmd: str, description: str = "", timeout: int = 3600) -> str:
        """执行命令并返回输出 | Execute command and return output"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.debug(f"命令 | Command: {cmd}")
        
        try:
            # 改变到输出目录执行命令 | Change to output directory to execute command
            self.output_dir.mkdir(parents=True, exist_ok=True)
            
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8',
                timeout=timeout,
                cwd=str(self.output_dir)
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return ""
        except subprocess.TimeoutExpired as e:
            self.logger.error(f"命令执行超时 | Command execution timeout: {description}")
            self.logger.error(f"超时时间 | Timeout: {timeout} seconds")
            return ""
        except Exception as e:
            self.logger.error(f"命令执行异常 | Command execution exception: {description}")
            self.logger.error(f"异常信息 | Exception: {e}")
            return ""
    
    def check_executable(self, executable_path: str) -> bool:
        """检查可执行文件是否可用 | Check if executable is available"""
        try:
            result = subprocess.run(
                [executable_path, "--version"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            # 某些工具可能不支持--version，尝试--help或-h
            try:
                result = subprocess.run(
                    [executable_path, "--help"], 
                    capture_output=True, 
                    text=True, 
                    timeout=10
                )
                return result.returncode == 0
            except:
                # 最后尝试直接运行工具名称，检查是否存在于PATH中
                return shutil.which(executable_path) is not None

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    # 依赖软件配置：(路径, 名称, 版本命令参数)
    dependencies = [
        (config.vcftools_path, "VCFtools", "--version"),
        (config.plink_path, "PLINK", "--version"),
        (config.bcftools_path, "BCFtools", "--version")
    ]
    
    # 如果需要计算有效群体大小，添加SMC++依赖
    if config.calculate_ne:
        dependencies.append((config.smcpp_path, "SMC++", "version"))  # SMC++使用"version"而不是"--version"
    
    missing_deps = []
    
    for cmd, name, version_arg in dependencies:
        try:
            # 根据不同工具使用不同的版本检查命令
            result = subprocess.run([cmd, version_arg], 
                                  capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
                
                # 可选：显示版本信息（用于调试）
                if logger.level <= 10:  # DEBUG级别
                    version_info = result.stdout.strip().split('\n')[0] if result.stdout else "版本信息不可用"
                    logger.debug(f"  {name} 版本: {version_info}")
                    
            else:
                logger.warning(f"✗ {name} 命令执行失败，返回码: {result.returncode}")
                if result.stderr:
                    logger.debug(f"  错误信息: {result.stderr.strip()}")
                missing_deps.append(name)
                
        except subprocess.TimeoutExpired:
            logger.warning(f"✗ {name} 命令执行超时")
            missing_deps.append(name)
            
        except FileNotFoundError:
            logger.warning(f"✗ {name} 未找到，请检查安装路径: {cmd}")
            missing_deps.append(name)
            
        except Exception as e:
            logger.warning(f"✗ {name} 检查时发生意外错误: {e}")
            missing_deps.append(name)
    
    # 报告结果
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        
        # 提供安装建议
        logger.error("安装建议 | Installation suggestions:")
        for dep in missing_deps:
            if dep == "VCFtools":
                logger.error("  VCFtools: conda install -c conda-forge vcftools")
            elif dep == "PLINK":
                logger.error("  PLINK: conda install -c conda-forge plink")
            elif dep == "BCFtools":
                logger.error("  BCFtools: conda install -c conda-forge bcftools")
            elif dep == "SMC++":
                logger.error("  SMC++: pip install smcpp")
        
        return False
    
    logger.info("✓ 所有依赖软件检查通过 | All dependencies check passed")
    return True

def format_time(seconds):
    """格式化时间 | Format time"""
    if seconds < 60:
        return f"{seconds:.1f}秒 | seconds"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}分钟 | minutes"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}小时 | hours"

def parse_window_sizes(window_str: str):
    """解析窗口大小字符串 | Parse window sizes string"""
    try:
        if ',' in window_str:
            return [int(w.strip()) for w in window_str.split(',')]
        else:
            return [int(window_str)]
    except ValueError as e:
        raise ValueError(f"无效的窗口大小格式 | Invalid window size format: {window_str}")

def validate_file_exists(file_path: str, file_type: str = "文件"):
    """验证文件是否存在 | Validate file exists"""
    if not Path(file_path).exists():
        raise FileNotFoundError(f"{file_type}不存在 | {file_type} does not exist: {file_path}")
    return True

def get_file_size_mb(file_path: str) -> float:
    """获取文件大小(MB) | Get file size in MB"""
    return Path(file_path).stat().st_size / (1024 * 1024)

def create_output_directory(output_dir: str) -> Path:
    """创建输出目录 | Create output directory"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path