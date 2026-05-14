"""
BLAST分析工具函数模块|BLAST Analysis Utility Functions Module
"""

import os
import subprocess
import shutil
import re
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in a conda environment

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
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


def check_tool_availability(tool_path: str) -> bool:
    """检查工具是否可用|Check if a tool is available in the system"""
    return shutil.which(tool_path) is not None


def get_input_files(input_path: str, suffix: str = "*.fa") -> List[str]:
    """获取输入文件列表|Get list of input files

    Args:
        input_path: 输入文件或目录路径|Input file or directory path
        suffix: 文件后缀模式|File suffix pattern for directory scanning

    Returns:
        文件路径列表|List of file paths
    """
    input_path_obj = Path(input_path)

    if input_path_obj.is_file():
        return [str(input_path_obj)]

    if input_path_obj.is_dir():
        if suffix.startswith('*'):
            pattern = suffix[1:]
            files = list(input_path_obj.glob(f"*{pattern}"))
        else:
            files = list(input_path_obj.glob(suffix))
        return [str(f) for f in files if f.is_file()]

    return []


def extract_sample_name_from_path(file_path: str, pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$') -> str:
    """从文件路径中提取样品名称|Extract sample name from file path

    Args:
        file_path: 文件路径|File path
        pattern: 样品名称提取正则表达式|Regex pattern for extracting sample name

    Returns:
        样品名称|Sample name
    """
    filename = os.path.basename(file_path)

    match = re.search(pattern, filename)
    if match:
        return match.group(1)

    return Path(filename).stem


def create_sample_map_file(input_files: List[str], output_file: str, pattern: str) -> int:
    """创建样品映射文件|Create sample mapping file

    Args:
        input_files: 输入文件路径列表|List of input file paths
        output_file: 输出映射文件路径|Output sample map file path
        pattern: 样品名称提取正则表达式|Sample name extraction pattern

    Returns:
        样品数量|Number of samples in the map file
    """
    from datetime import datetime

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Auto-generated sample mapping file\n")
        f.write(f"# Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# Format: file_path<TAB>sample_name\n")
        f.write("#" + "="*60 + "\n")

        sample_count = 0
        for input_file in sorted(input_files):
            sample_name = extract_sample_name_from_path(input_file, pattern)
            f.write(f"{input_file}\t{sample_name}\n")
            sample_count += 1

    return sample_count


def load_sample_mapping(map_file: str) -> dict:
    """加载样品映射文件|Load sample mapping from file

    Args:
        map_file: 样品映射文件路径|Sample mapping file path

    Returns:
        文件名到样品名称的映射字典|Dictionary mapping filename to sample name
    """
    sample_mapping = {}
    file_paths = []

    if not os.path.exists(map_file):
        raise FileNotFoundError(f"样品映射文件不存在|Sample mapping file not found: {map_file}")

    with open(map_file, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                file_path = parts[0]
                sample_name = parts[1]

                if not os.path.exists(file_path):
                    continue

                filename = os.path.basename(file_path)
                sample_mapping[filename] = sample_name
                file_paths.append(file_path)

    return {
        'mapping': sample_mapping,
        'file_paths': file_paths
    }


def format_command(cmd: list) -> str:
    """格式化命令列表为字符串|Format command list to string"""
    return ' '.join(str(x) for x in cmd)


def run_command(cmd: list, logger=None, description: str = "") -> bool:
    """执行命令|Execute command

    Args:
        cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
        logger: 日志器对象|Logger object
        description: 命令描述|Command description

    Returns:
        是否成功|True if successful, False otherwise
    """
    if logger and description:
        logger.info(f"执行|Executing: {description}")

    if logger:
        logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        if logger and description:
            logger.info(f"命令执行成功|Command completed successfully: {description}")
        return True

    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"命令执行失败|Command execution failed: {description}")
            logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                logger.error(f"错误信息|Error output: {e.stderr}")
        return False


def check_blast_dependencies(config, logger) -> bool:
    """检查BLAST依赖软件|Check BLAST dependencies

    Args:
        config: BLAST配置对象|BLAST configuration object
        logger: 日志器对象|Logger object

    Returns:
        是否全部可用|True if all dependencies available

    Raises:
        RuntimeError: 如果有依赖缺失|If any dependency is missing
    """
    logger.info("检查BLAST依赖软件|Checking BLAST dependencies")

    blast_program = getattr(config, f'{config.blast_type}_path', config.blast_type)

    dependencies = [
        (config.makeblastdb_path, "makeblastdb"),
        (blast_program, config.blast_type.upper())
    ]

    missing_deps = []

    for cmd, name in dependencies:
        if check_tool_availability(cmd):
            logger.info(f"  {name}: 可用|available")
        else:
            logger.warning(f"  {name}: 未找到|NOT found")
            missing_deps.append(name)

    if missing_deps:
        error_msg = f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    logger.info("所有依赖软件可用|All dependencies available")
    return True
