"""
MSA工具函数模块|MSA Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class MSALogger:
    """MSA日志管理器|MSA Logger Manager"""
    
    def __init__(self, output_prefix: str, log_name: str = "msa_analysis.log"):
        self.log_file = Path(f"{output_prefix}.log")
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def run(self, cmd, description: str = "", output_file: str = None) -> bool:
        """执行命令|Execute command"""
        from ..common.conda_runner import build_conda_command

        if isinstance(cmd, list):
            full_cmd = build_conda_command(str(cmd[0]), [str(c) for c in cmd[1:]])
            shell = False
            display = ' '.join(full_cmd)
        else:
            full_cmd = cmd
            shell = True
            display = str(cmd)

        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")
        
        self.logger.info(f"命令|Command: {display}")
        
        try:
            if output_file:
                with open(output_file, 'w') as f:
                    result = subprocess.run(
                        full_cmd, shell=shell, stdout=f, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    self.logger.error(f"命令执行失败|Command execution failed: {description}")
                    if result.stderr:
                        self.logger.error(f"错误信息|Error message: {result.stderr}")
                    return False
                self.logger.info(f"命令执行成功|Command executed successfully")
                return True
            else:
                result = subprocess.run(
                    full_cmd,
                    shell=shell,
                    capture_output=True,
                    text=True,
                    check=True
                )
            
            self.logger.info(f"命令执行成功|Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    from ..common.conda_runner import check_tools

    # 根据选择的方法检查对应工具|Check tool based on selected method
    tool_map = {
        'mafft': (config.mafft_path, "--version"),
        'clustalo': (config.clustalo_path, "--version"),
        'muscle': (config.muscle_path, "-version"),
        't_coffee': (config.tcoffee_path, "-version")
    }

    tool_cmd, version_arg = tool_map.get(config.method, (None, None))
    if tool_cmd is None:
        error_msg = f"缺少依赖软件|Missing dependency: {config.method.upper()}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    missing = check_tools([(tool_cmd, config.method, [version_arg])], logger)
    if missing:
        error_msg = f"缺少依赖软件|Missing dependency: {config.method.upper()}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

def count_sequences(fasta_file: str) -> int:
    """统计FASTA文件中的序列数|Count sequences in FASTA file"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count
