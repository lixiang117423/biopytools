"""
BLAST分析工具函数模块|BLAST Analysis Utility Functions Module
仅保留conda环境包装函数;样品映射/命令执行等逻辑由main.py的BLASTAnalyzer负责
|Only conda wrapping kept here; sample mapping / command execution live in BLASTAnalyzer (main.py)
"""

import os
import shutil
import re
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中,返回环境名称|Detect if command is in a conda environment

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建conda run命令来运行conda环境中的软件|Build conda run command for conda env software

    Args:
        command: 命令名称或完整路径(建议传完整路径以正确检测env)|Command name or full path (full path recommended for correct env detection)
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表(配合subprocess.run(shell=False))|Complete command list (use with subprocess.run(shell=False))

    Note:
        必须使用--no-capture-output避免conda缓冲输出导致大数据内存溢出
        |Must use --no-capture-output to avoid conda buffering output causing OOM on large data
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args
