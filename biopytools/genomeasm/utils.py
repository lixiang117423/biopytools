"""
基因组组装工具函数模块 | Genome Assembly Utility Functions Module 🛠️
"""

import logging
import subprocess
import sys
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

class AssemblyLogger:
    """基因组组装日志管理器 | Genome Assembly Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "assembly.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / "logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
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
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.commands_log = working_dir / "logs" / "commands.log"
        self.commands_log.parent.mkdir(parents=True, exist_ok=True)
    
    def run(self, cmd: str, description: str = "", timeout: Optional[int] = None) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"📝 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        # 记录命令到日志文件
        with open(self.commands_log, 'a', encoding='utf-8') as f:
            f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {description}\n")
            f.write(f"Directory: {self.working_dir}\n")
            f.write(f"Command: {cmd}\n\n")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                timeout=timeout
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout[:500]}...")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"💥 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"📥 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"⏰ 命令执行超时 | Command execution timeout: {description}")
            return False

# def check_dependencies(config, logger) -> Dict[str, bool]:
#     """检查依赖软件 | Check dependencies"""
#     logger.info("🔍 检查依赖软件 | Checking dependencies")
    
#     dependencies = {
#         'hifiasm': config.hifiasm_path,
#         'bwa': config.bwa_path,
#         'samtools': config.samtools_path,
#         'seqkit': 'seqkit',
#         'fastqc': 'fastqc',
#     }
    
#     # Hi-C相关工具检查
#     if 'Hi-C' in config.detected_data_types:
#         if config.hic_strategy == "complete_juicer":
#             dependencies['juicer'] = config.juicer_path
#             dependencies['juicer_tools'] = config.juicer_tools
#         elif config.hic_strategy == "standard_3ddna":
#             dependencies['3d-dna'] = config.pipeline_3ddna
#             dependencies['juicer_tools'] = config.juicer_tools
#         elif config.hic_strategy == "simplified_salsa2":
#             dependencies['salsa2'] = config.salsa2_path
    
#     results = {}
#     missing_deps = []
    
#     for name, cmd in dependencies.items():
#         try:
#             if name == 'juicer_tools':
#                 # Java jar文件特殊处理
#                 if os.path.exists(cmd):
#                     result = subprocess.run(['java', '-jar', cmd], 
#                                           capture_output=True, text=True, timeout=10)
#                     results[name] = True
#                     logger.info(f"✅ {name} 可用 | available")
#                 else:
#                     results[name] = False
#                     missing_deps.append(name)
#             elif name == '3d-dna':
#                 # 检查脚本文件是否存在且可执行
#                 if os.path.exists(cmd) and os.access(cmd, os.X_OK):
#                     results[name] = True
#                     logger.info(f"✅ {name} 可用 | available")
#                 else:
#                     results[name] = False
#                     missing_deps.append(name)
#             else:
#                 # 标准命令检查
#                 result = subprocess.run([cmd, "--version"], 
#                                       capture_output=True, text=True, timeout=10)
#                 if result.returncode == 0:
#                     results[name] = True
#                     logger.info(f"✅ {name} 可用 | available")
#                 else:
#                     results[name] = False
#                     missing_deps.append(name)
#         except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
#             results[name] = False
#             missing_deps.append(name)
    
#     if missing_deps:
#         error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
#         logger.error(error_msg)
#         logger.info("💡 请检查软件安装和PATH设置 | Please check software installation and PATH settings")
#         raise RuntimeError(error_msg)
    
#     logger.info("🎉 所有依赖软件检查通过 | All dependencies check passed")
#     return results

def check_dependencies(config, logger) -> Dict[str, bool]:
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    # 确保detected_data_types不为None
    if not hasattr(config, 'detected_data_types') or config.detected_data_types is None:
        config.detected_data_types = ['HiFi']  # 默认假设有HiFi
        logger.warning("⚠️ 数据类型检测异常，假设存在HiFi数据")
    
    dependencies = {
        'hifiasm': config.hifiasm_path,
        'bwa': config.bwa_path,
        'samtools': config.samtools_path,
        'seqkit': 'seqkit',
        'fastqc': 'fastqc',
    }
    
    # Hi-C相关工具检查 - 添加安全检查
    if 'Hi-C' in config.detected_data_types:
        if config.hic_strategy == "complete_juicer":
            dependencies['juicer'] = config.juicer_path
            dependencies['juicer_tools'] = config.juicer_tools
        elif config.hic_strategy == "standard_3ddna":
            dependencies['3d-dna'] = config.pipeline_3ddna
            dependencies['juicer_tools'] = config.juicer_tools
        elif config.hic_strategy == "simplified_salsa2":
            dependencies['salsa2'] = config.salsa2_path
    
    results = {}
    missing_deps = []
    
    for name, cmd in dependencies.items():
        try:
            if name == 'juicer_tools':
                # Java jar文件特殊处理
                success = _check_juicer_tools(cmd, logger)
                results[name] = success
                if success:
                    logger.info(f"✅ {name} 可用 | available")
                else:
                    missing_deps.append(name)
            elif name == '3d-dna':
                # 检查脚本文件是否存在且可执行
                if os.path.exists(cmd) and os.access(cmd, os.X_OK):
                    results[name] = True
                    logger.info(f"✅ {name} 可用 | available")
                else:
                    results[name] = False
                    missing_deps.append(name)
            elif name == 'juicer':
                # Juicer特殊检查
                success = _check_juicer(cmd, logger)
                results[name] = success
                if success:
                    logger.info(f"✅ {name} 可用 | available")
                else:
                    missing_deps.append(name)
            else:
                # 标准命令检查 - 改进版本
                success = _check_standard_command(name, cmd, logger)
                results[name] = success
                if success:
                    logger.info(f"✅ {name} 可用 | available")
                else:
                    missing_deps.append(name)
        except Exception as e:
            logger.warning(f"⚠️ 检查 {name} 时出错: {e}")
            results[name] = False
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        logger.info("💡 请检查软件安装和PATH设置 | Please check software installation and PATH settings")
        logger.info("🔧 或使用以下参数指定软件路径 | Or use following parameters to specify tool paths:")
        for dep in missing_deps:
            if dep == 'bwa':
                logger.info(f"   --bwa-path /path/to/bwa")
            elif dep == 'seqkit':
                logger.info(f"   --seqkit-path /path/to/seqkit (需要添加此参数)")
            elif dep == 'juicer':
                logger.info(f"   --juicer-path /path/to/juicer.sh")
            elif dep == 'juicer_tools':
                logger.info(f"   --juicer-tools /path/to/juicer_tools.jar")
        
        # 提供跳过选项的建议
        logger.info("🚫 如需跳过依赖检查，请添加 --skip-dependency-check 参数")
        raise RuntimeError(error_msg)
    
    logger.info("🎉 所有依赖软件检查通过 | All dependencies check passed")
    return results

def _check_standard_command(name: str, cmd: str, logger) -> bool:
    """检查标准命令 | Check standard command"""
    # 不同软件的测试参数
    test_params = {
        'bwa': ['-h', '--help', ''],  # bwa直接调用显示帮助
        'seqkit': ['-h', '--help', ''],  # seqkit直接调用显示帮助
        'samtools': ['--version', '--help'],
        'hifiasm': ['--version', '-h'],
        'fastqc': ['--version', '--help']
    }
    
    params_to_try = test_params.get(name, ['--version', '--help', '-h'])
    
    for param in params_to_try:
        try:
            if param == '':
                # 直接调用命令（某些软件直接调用会显示帮助）
                result = subprocess.run([cmd], capture_output=True, text=True, timeout=10)
            else:
                result = subprocess.run([cmd, param], capture_output=True, text=True, timeout=10)
            
            # 对于某些软件，返回码不为0也可能是正常的（如显示帮助信息）
            if result.returncode == 0 or (param in ['', '-h', '--help'] and result.stderr):
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            continue
    
    return False

def _check_juicer(juicer_cmd: str, logger) -> bool:
    """检查Juicer软件 | Check Juicer software"""
    # 可能的juicer命令名
    possible_commands = [
        juicer_cmd,
        'juicer',
        'juicer.sh',
        os.path.expanduser('~/miniconda3/envs/*/bin/juicer.sh'),  # conda路径
        '/opt/juicer/scripts/juicer.sh'  # 常见安装路径
    ]
    
    for cmd in possible_commands:
        try:
            if '*' in cmd:  # 处理通配符路径
                import glob
                for expanded_cmd in glob.glob(cmd):
                    if os.path.exists(expanded_cmd) and os.access(expanded_cmd, os.X_OK):
                        return True
            else:
                if os.path.exists(cmd) and os.access(cmd, os.X_OK):
                    return True
                # 尝试作为命令执行
                result = subprocess.run([cmd], capture_output=True, text=True, timeout=5)
                # Juicer通常直接调用会显示使用说明
                if result.returncode in [0, 1] or 'juicer' in result.stderr.lower():
                    return True
        except:
            continue
    
    return False

def _check_juicer_tools(juicer_tools_path: str, logger) -> bool:
    """检查Juicer Tools | Check Juicer Tools"""
    # 方式1：作为jar文件检查
    if os.path.exists(juicer_tools_path) and juicer_tools_path.endswith('.jar'):
        try:
            result = subprocess.run(['java', '-jar', juicer_tools_path], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 or 'juicer' in result.stderr.lower():
                return True
        except:
            pass
    
    # 方式2：作为直接可执行命令检查
    try:
        result = subprocess.run([juicer_tools_path, '-h'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 or 'juicer' in result.stderr.lower():
            return True
    except:
        pass
    
    # 方式3：寻找可能的juicer_tools位置
    possible_paths = [
        'juicer_tools',
        'juicer_tools.jar',
        os.path.expanduser('~/miniconda3/envs/*/bin/juicer_tools*'),
        '/opt/juicer/CPU/common/juicer_tools.jar'
    ]
    
    for path in possible_paths:
        try:
            if '*' in path:
                import glob
                for expanded_path in glob.glob(path):
                    if _check_juicer_tools(expanded_path, logger):
                        return True
            else:
                if _check_juicer_tools(path, logger):
                    return True
        except:
            continue
    
    return False

def get_file_stats(file_path: str, logger) -> Dict[str, str]:
    """获取文件统计信息 | Get file statistics"""
    try:
        cmd = f"seqkit stats {file_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        
        # 解析seqkit输出
        lines = result.stdout.strip().split('\n')
        if len(lines) >= 2:
            headers = lines[0].split()
            values = lines[1].split()
            stats = dict(zip(headers, values))
            return stats
        else:
            return {}
    except subprocess.CalledProcessError:
        logger.warning(f"⚠️ 无法获取文件统计 | Cannot get file stats: {file_path}")
        return {}

def estimate_genome_coverage(reads_file: str, genome_size: str, logger) -> float:
    """估算基因组覆盖度 | Estimate genome coverage"""
    try:
        # 获取总碱基数
        stats = get_file_stats(reads_file, logger)
        if 'sum_len' not in stats:
            return 0.0
        
        total_bases = int(stats['sum_len'].replace(',', ''))
        
        # 解析基因组大小
        genome_size = genome_size.lower()
        if genome_size.endswith('k'):
            genome_bp = int(float(genome_size[:-1]) * 1000)
        elif genome_size.endswith('m'):
            genome_bp = int(float(genome_size[:-1]) * 1000000)
        elif genome_size.endswith('g'):
            genome_bp = int(float(genome_size[:-1]) * 1000000000)
        else:
            genome_bp = int(genome_size)
        
        coverage = total_bases / genome_bp
        logger.info(f"📊 预估覆盖度 | Estimated coverage: {coverage:.1f}X")
        return coverage
        
    except Exception as e:
        logger.warning(f"⚠️ 无法估算覆盖度 | Cannot estimate coverage: {e}")
        return 0.0

def create_directory_structure(output_dir: Path, logger):
    """创建目录结构 | Create directory structure"""
    directories = [
        'logs',
        'qc', 
        'assembly',
        'hic_processing',
        'results',
        'temp'
    ]
    
    for dir_name in directories:
        dir_path = output_dir / dir_name
        dir_path.mkdir(parents=True, exist_ok=True)
        logger.debug(f"📁 创建目录 | Created directory: {dir_path}")

def format_time(seconds: float) -> str:
    """格式化时间 | Format time"""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def format_size(bytes_size: int) -> str:
    """格式化文件大小 | Format file size"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024:
            return f"{bytes_size:.1f}{unit}"
        bytes_size /= 1024
    return f"{bytes_size:.1f}PB"
