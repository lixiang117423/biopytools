"""
🌾 EDTA工具函数模块 | EDTA Utility Functions Module
"""

import os
import logging
import subprocess
import sys
import shutil
import time
import signal
from pathlib import Path
from typing import Optional, Dict
from typing import List
from datetime import datetime
from typing import List
import json

class EDTALogger:
    """EDTA分析日志管理器 | EDTA Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "edta_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            backup_file = self.log_file.with_suffix(f".backup.{int(time.time())}.log")
            shutil.copy2(self.log_file, backup_file)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # 记录开始信息
        # self.logger.info("🌾 =" * 50)
        self.logger.info("🌾 EDTA植物基因组TE注释工具启动 | EDTA Plant Genome TE Annotation Tool Started")
        self.logger.info(f"🌾 日志文件 | Log file: {self.log_file}")
        # self.logger.info("🌾 =" * 50)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.process = None
        self.start_time = None
    
    def run(self, cmd: List[str], description: str = "", check: bool = True, 
            show_progress: bool = True) -> subprocess.CompletedProcess:
        """执行命令 | Execute command"""
        
        if description:
            self.logger.info(f"📋 执行步骤 | Executing step: {description}")
        
        cmd_str = " ".join(cmd) if isinstance(cmd, list) else cmd
        self.logger.info(f"🔧 命令 | Command: {cmd_str}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        self.start_time = time.time()
        
        try:
            # 实时输出模式
            self.process = subprocess.Popen(
                cmd,
                shell=False if isinstance(cmd, list) else True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=self.working_dir,
                universal_newlines=True
            )
            
            # 实时显示输出
            output_lines = []
            if show_progress:
                self.logger.info("📊 开始监控进度 | Starting progress monitoring...")
            
            # while True:
            #     output = self.process.stdout.readline()
            #     if output == '' and self.process.poll() is not None:
            #         break
            #     if output:
            #         output_lines.append(output.strip())
            #         if show_progress:
            #             self.logger.info(f"📋 EDTA: {output.strip()}")
            while True:
                output = self.process.stdout.readline()
                if output == '' and self.process.poll() is not None:
                    break
                if output:
                    line = output.strip()
                    output_lines.append(line)
                    if show_progress:
                        # 对于mv、cp、ln等命令，显示完整命令行
                        if any(cmd in line for cmd in ['mv:', 'cp:', 'ln:', 'rm:']):
                            self.logger.info(f"📋 EDTA: {line}")
                        # 对于包含完整路径的错误信息，也完整显示
                        elif 'No such file or directory' in line:
                            self.logger.info(f"📋 EDTA: {line}")
                        # 其他信息正常显示
                        else:
                            self.logger.info(f"📋 EDTA: {line}")
                        
            self.process.wait()
            
            if check and self.process.returncode != 0:
                raise subprocess.CalledProcessError(
                    self.process.returncode, cmd, '\n'.join(output_lines)
                )
            
            # 创建结果对象
            result = subprocess.CompletedProcess(
                cmd, self.process.returncode, 
                stdout='\n'.join(output_lines), stderr=''
            )
            
            elapsed_time = time.time() - self.start_time
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            self.logger.info(f"⏱️ 耗时 | Time elapsed: {elapsed_time:.2f} 秒 | seconds")
            
            return result
            
        except subprocess.CalledProcessError as e:
            elapsed_time = time.time() - self.start_time if self.start_time else 0
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"⏱️ 失败前耗时 | Time before failure: {elapsed_time:.2f} 秒 | seconds")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            
            if hasattr(e, 'stdout') and e.stdout:
                self.logger.error(f"📄 标准输出 | Stdout: {e.stdout}")
            
            if check:
                raise
            return e
    
    def terminate_process(self):
        """终止当前进程 | Terminate current process"""
        if self.process and self.process.poll() is None:
            self.logger.warning("🛑 终止EDTA进程 | Terminating EDTA process...")
            self.process.terminate()
            try:
                self.process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self.logger.warning("⚡ 强制杀死EDTA进程 | Force killing EDTA process...")
                self.process.kill()

def check_edta_dependencies(config, logger) -> bool:
    """检查EDTA依赖软件 | Check EDTA dependencies"""
    logger.info("🔍 检查EDTA依赖软件 | Checking EDTA dependencies")
    
    # 首先检查EDTA是否存在
    try:
        result = subprocess.run(["EDTA.pl", "--help"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info("✅ EDTA 可用 | EDTA available")
        else:
            logger.error("❌ EDTA 不可用 | EDTA not available")
            raise RuntimeError("EDTA 软件未找到或无法运行 | EDTA software not found or cannot run")
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"❌ EDTA 检查失败 | EDTA check failed: {e}")
        raise RuntimeError("EDTA 软件未找到或无法运行 | EDTA software not found or cannot run")
    
    logger.info("✅ 依赖检查完成 | Dependencies check completed")
    return True

def setup_output_directories(output_dir: Path, logger) -> Dict[str, Path]:
    """设置输出目录结构 | Setup output directory structure"""
    logger.info("📁 设置输出目录结构 | Setting up output directory structure")
    
    # 定义目录结构
    directories = {
        "input": output_dir / "00_input",
        "edta_raw": output_dir / "01_edta_raw", 
        "processed": output_dir / "02_processed",
        "reports": output_dir / "03_reports",
        "visualization": output_dir / "04_visualization",
        "logs": output_dir / "logs"
    }
    
    # 创建目录
    for dir_name, dir_path in directories.items():
        dir_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"📂 创建目录 | Created directory: {dir_name} -> {dir_path}")
    
    return directories

def check_resume_capability(output_dir: Path, logger) -> Optional[Dict]:
    """检查断点续跑能力 | Check resume capability"""
    logger.info("🔄 检查断点续跑状态 | Checking resume capability")
    
    resume_file = output_dir / "logs" / "resume_checkpoint.json"
    
    if resume_file.exists():
        try:
            with open(resume_file, 'r', encoding='utf-8') as f:
                checkpoint = json.load(f)
            
            logger.info(f"📋 找到检查点文件 | Found checkpoint file: {resume_file}")
            logger.info(f"📅 上次运行时间 | Last run time: {checkpoint.get('timestamp', 'Unknown')}")
            logger.info(f"📍 上次步骤 | Last step: {checkpoint.get('last_step', 'Unknown')}")
            
            return checkpoint
            
        except Exception as e:
            logger.warning(f"⚠️ 检查点文件读取失败 | Failed to read checkpoint file: {e}")
            return None
    else:
        logger.info("📋 未找到检查点文件，将进行全新分析 | No checkpoint file found, will perform fresh analysis")
        return None

def save_checkpoint(output_dir: Path, step: str, data: Dict, logger):
    """保存检查点 | Save checkpoint"""
    logger.debug(f"💾 保存检查点 | Saving checkpoint: {step}")
    
    resume_file = output_dir / "logs" / "resume_checkpoint.json"
    
    checkpoint = {
        "timestamp": datetime.now().isoformat(),
        "last_step": step,
        "data": data
    }
    
    try:
        with open(resume_file, 'w', encoding='utf-8') as f:
            json.dump(checkpoint, f, indent=2, ensure_ascii=False)
        logger.debug(f"✅ 检查点已保存 | Checkpoint saved: {resume_file}")
    except Exception as e:
        logger.warning(f"⚠️ 检查点保存失败 | Failed to save checkpoint: {e}")
