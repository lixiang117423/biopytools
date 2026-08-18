"""
重复序列分析工具函数模块|Repeat Sequence Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

from ..common.conda_runner import build_conda_command, check_tools

from ..common.conda_runner import build_conda_command, check_tools

class RepeatLogger:
    """ 重复序列分析日志管理器|Repeat Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "repeat_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """ 设置日志|Setup logging"""
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
        """ 获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """ 命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd, description: str = "", timeout: Optional[int] = None) -> bool:
        """ 执行命令|Execute command"""
        if description:
            self.logger.info(f" 执行步骤|Executing step: {description}")
        
        # 列表命令用 shell=False, 字符串命令用 shell=True|List -> shell=False, string -> shell=True
        if isinstance(cmd, (list, tuple)):
            use_shell = False
            display = ' '.join(str(c) for c in cmd)
        else:
            use_shell = True
            display = str(cmd)
        self.logger.info(f" 命令|Command: {display}")
        self.logger.info(f" 工作目录|Working directory: {self.working_dir}")
        
        start_time = time.time()
        
        try:
            result = subprocess.run(
                cmd, 
                shell=use_shell,
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                timeout=timeout
            )
            
            elapsed_time = time.time() - start_time
            self.logger.info(f" 命令执行成功|Command executed successfully: {description} ({elapsed_time:.2f}s)")
            
            if result.stdout:
                self.logger.debug(f" 标准输出|Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            elapsed_time = time.time() - start_time
            self.logger.error(f" 命令执行失败|Command execution failed: {description} ({elapsed_time:.2f}s)")
            self.logger.error(f" 错误代码|Error code: {e.returncode}")
            self.logger.error(f" 错误信息|Error message: {e.stderr}")
            self.logger.error(f" 标准输出|Stdout: {e.stdout}")
            return False
        
        except subprocess.TimeoutExpired:
            elapsed_time = time.time() - start_time
            self.logger.error(f" 命令执行超时|Command execution timeout: {description} ({elapsed_time:.2f}s)")
            return False

class ProgressTracker:
    """ 进度跟踪器|Progress Tracker"""
    
    def __init__(self, logger, total_steps: int = 6):
        self.logger = logger
        self.total_steps = total_steps
        self.current_step = 0
        self.step_names = [
            " 依赖检查|Dependency Check",
            " 从头重复库构建|De novo Repeat Library Construction", 
            " LTR序列分析|LTR Sequence Analysis",
            " 重复库合并|Repeat Library Merging",
            " 重复序列识别|Repeat Sequence Identification",
            " 转座元件分类|Transposable Element Classification"
        ]
    
    def start_step(self, step_name: str = None):
        """开始新步骤|Start new step"""
        self.current_step += 1
        if step_name is None and self.current_step <= len(self.step_names):
            step_name = self.step_names[self.current_step - 1]
        
        progress = (self.current_step / self.total_steps) * 100
        self.logger.info(f" 进度 {progress:.1f}% ({self.current_step}/{self.total_steps})|Progress: {step_name}")
    
    def complete(self):
        """完成所有步骤|Complete all steps"""
        self.logger.info(" 分析完成！| Analysis completed!")

def check_dependencies(config, logger):
    """ 检查依赖软件|Check dependencies"""
    dependencies = []
    
    if not config.skip_modeler:
        dependencies.append((config.repeatmodeler_path, "RepeatModeler", ["-h"], True))
        dependencies.append((config.build_database_path, "BuildDatabase", ["-h"], True))
    
    if not config.skip_ltr:
        dependencies.extend([
            (config.ltr_finder_path, "LTR_FINDER", ["-h"], True),
            (config.ltrharvest_path, "LTRharvest (GenomeTools)", ["ltrharvest", "--help"], True),
            (config.ltr_retriever_path, "LTR_retriever", ["-h"], True)
        ])
    
    dependencies.extend([
        (config.repeatmasker_path, "RepeatMasker", ["-v"], True),
        (config.tesorter_path, "TEsorter", ["-h"], True)
    ])
    
    missing_deps = check_tools(dependencies, logger)
    
    if missing_deps:
        error_msg = f" 缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    logger.info(" 所有依赖软件检查通过|All dependencies check passed")
    return True

def get_genome_stats(genome_file: str, logger) -> Dict[str, any]:
    """ 获取基因组统计信息|Get genome statistics"""
    logger.info(" 分析基因组文件|Analyzing genome file")
    
    stats = {
        'sequences': 0,
        'total_length': 0,
        'n50': 0,
        'max_length': 0,
        'min_length': float('inf'),
        'gc_content': 0
    }
    
    try:
        with open(genome_file, 'r') as f:
            current_seq_length = 0
            lengths = []
            gc_count = 0
            total_bases = 0
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq_length > 0:
                        lengths.append(current_seq_length)
                        stats['total_length'] += current_seq_length
                        stats['max_length'] = max(stats['max_length'], current_seq_length)
                        stats['min_length'] = min(stats['min_length'], current_seq_length)
                        current_seq_length = 0
                    stats['sequences'] += 1
                else:
                    current_seq_length += len(line)
                    gc_count += line.upper().count('G') + line.upper().count('C')
                    total_bases += len(line)
            
            # 处理最后一个序列|Process last sequence
            if current_seq_length > 0:
                lengths.append(current_seq_length)
                stats['total_length'] += current_seq_length
                stats['max_length'] = max(stats['max_length'], current_seq_length)
                stats['min_length'] = min(stats['min_length'], current_seq_length)
            
            # 计算N50|Calculate N50
            if lengths:
                lengths.sort(reverse=True)
                total_half = stats['total_length'] // 2
                cumulative = 0
                for length in lengths:
                    cumulative += length
                    if cumulative >= total_half:
                        stats['n50'] = length
                        break
            
            # 计算GC含量|Calculate GC content
            if total_bases > 0:
                stats['gc_content'] = (gc_count / total_bases) * 100
            
            # 处理空基因组情况|Handle empty genome case
            if stats['min_length'] == float('inf'):
                stats['min_length'] = 0
        
        # 格式化输出大小|Format size output
        size_mb = stats['total_length'] / 1_000_000
        logger.info(f" 基因组统计|Genome statistics:")
        logger.info(f"   序列数|Sequences: {stats['sequences']:,}")
        logger.info(f"   总长度|Total length: {stats['total_length']:,} bp ({size_mb:.2f} MB)")
        logger.info(f"   N50: {stats['n50']:,} bp")
        logger.info(f"   最大序列|Max sequence: {stats['max_length']:,} bp")
        logger.info(f"   最小序列|Min sequence: {stats['min_length']:,} bp")
        logger.info(f"   GC含量|GC content: {stats['gc_content']:.2f}%")
        
    except Exception as e:
        logger.error(f" 基因组统计分析失败|Genome statistics analysis failed: {e}")
    
    return stats
