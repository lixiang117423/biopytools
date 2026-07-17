"""
基因组组装工具函数模块|Genome Assembly Utility Functions Module

命令执行统一走 conda 包装(§13):单命令 run_command(list+shell=False)、
管道命令 run_pipeline_command(剥conda run+LD_LIBRARY_PATH+shell=True)。
|Command execution unified through conda wrapping (§13): single commands via
run_command (list+shell=False), pipelines via run_pipeline_command (strip conda
run + LD_LIBRARY_PATH + shell=True).
"""

import subprocess
import os
import re
import sys
import shutil
import logging
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中,返回环境名称|Detect if command is in conda env, return env name

    策略|Strategy:
    1. 完整路径含/envs/直接提取(最高优先)|Full path with /envs/ extracted directly
    2. which路径检测|which path detection
    3. 遍历所有conda env兜底|Scan all conda envs as fallback
    """
    # 方法0:完整路径含envs直接提取|Method 0: full path with envs
    if '/envs/' in command:
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 方法1:which路径检测|Method 1: which path
    cmd_path = shutil.which(command)
    if cmd_path:
        if os.path.islink(cmd_path):
            cmd_path = os.path.realpath(cmd_path)
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2:遍历所有conda env|Method 2: scan all envs
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, 'bin', command)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令(必须--no-capture-output,避免缓冲OOM)|Build conda run command

    传递完整路径(非命令名)以便get_conda_env正确提取env(§13.6.1)
    |Pass full path (not command name) so get_conda_env extracts env correctly (§13.6.1)
    """
    # 展开路径(用户可能传含~的路径)|Expand path (user may pass ~ path)
    if '/' in command:
        command = expand_path(command)
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def _conda_lib_env(tool_paths: List[str]) -> dict:
    """
    为管道工具构建LD_LIBRARY_PATH环境(§13.2.1方案A防御性补强)
    |Build LD_LIBRARY_PATH env for pipeline tools (§13.2.1 plan A defensive reinforcement)

    从各工具完整路径推导同级lib目录(dirname/../../lib),拼到LD_LIBRARY_PATH前端。
    |Derive sibling lib dir from each tool full path, prepend to LD_LIBRARY_PATH.
    """
    env = os.environ.copy()
    lib_dirs = []
    for tool_path in tool_paths:
        if not tool_path or '/' not in tool_path:
            continue
        tool_path = expand_path(tool_path)
        tool_dir = os.path.dirname(tool_path)
        lib_dir = os.path.normpath(os.path.join(tool_dir, '..', 'lib'))
        if os.path.isdir(lib_dir) and lib_dir not in lib_dirs:
            lib_dirs.append(lib_dir)
    if lib_dirs:
        existing = env.get('LD_LIBRARY_PATH', '')
        env['LD_LIBRARY_PATH'] = ':'.join(lib_dirs) + (f':{existing}' if existing else '')
    return env


def _build_pipeline_command(commands: List[List[str]]) -> str:
    """
    构建管道命令字符串(剥conda run留裸命令,避免conda run|conda run,§13.2.1方案B)
    |Build pipeline command string (strip conda run to bare commands, §13.2.1 plan B)
    """
    wrapped_parts = []
    for cmd in commands:
        if not cmd:
            continue
        full = build_conda_command(cmd[0], cmd[1:])
        # conda run -n env --no-capture-output command args... → 提取 command args(裸命令)
        if full[0] == 'conda' and len(full) > 5:
            wrapped_parts.append(' '.join(full[5:]))
        else:
            wrapped_parts.append(' '.join(full))
    return ' | '.join(wrapped_parts)


def _pipeline_tool_paths(commands: List[List[str]]) -> List[str]:
    """提取管道中各工具路径(用于LD_LIBRARY_PATH)|Extract tool paths from pipeline"""
    return [cmd[0] for cmd in commands if cmd]


def check_dependencies(logger, hifiasm_bin: str = '~/miniforge3/envs/hifiasm_v.0.25.0/bin/hifiasm') -> bool:
    """
    检查依赖软件|Check dependencies

    Args:
        logger: 日志器|Logger
        hifiasm_bin: hifiasm完整路径(支持~)|hifiasm full path (supports ~)
    """
    logger.info("检查依赖软件|Checking dependencies")
    try:
        cmd = build_conda_command(hifiasm_bin, ['-V'])
        logger.info(f"命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, shell=False)
        if result.returncode == 0:
            logger.info(f"hifiasm 可用|available: {result.stdout.strip()}")
            return True
        logger.error(f"hifiasm版本检查失败|hifiasm version check failed (rc={result.returncode})")
        raise RuntimeError("缺少必需软件: hifiasm|Missing required software: hifiasm")
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"hifiasm未安装或不可用|hifiasm not installed or unavailable: {e}")
        raise RuntimeError("缺少必需软件: hifiasm|Missing required software: hifiasm")


def run_command(cmd: List[str], logger, work_dir: str = None, capture_output: bool = False) -> bool:
    """
    执行单条命令(list形式,自动conda包装,shell=False,§13)|Run single command

    Args:
        cmd: 命令列表,cmd[0]为工具完整路径(支持~)|Command list, cmd[0] is tool full path
        work_dir: 工作目录|Working directory
        capture_output: 是否捕获输出|Whether to capture output
    """
    if not cmd:
        logger.error("空命令|Empty command")
        return False
    wrapped = build_conda_command(cmd[0], cmd[1:])
    cmd_str = ' '.join(wrapped)
    logger.info(f"命令|Command: {cmd_str}")
    try:
        result = subprocess.run(
            wrapped, shell=False, cwd=work_dir,
            capture_output=capture_output, text=True, check=True
        )
        if capture_output and result.stdout:
            logger.debug(f"输出|Output: {result.stdout}")
        logger.info("命令执行成功|Command executed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False
    except FileNotFoundError as e:
        logger.error(f"命令不存在|Command not found: {e}")
        return False


def run_pipeline_command(commands: List[List[str]], logger, work_dir: str = None) -> bool:
    """
    执行管道命令(跨/同env管道:剥conda run+LD_LIBRARY_PATH+shell=True,§13.2.1)
    |Run pipeline command (cross/same env: strip conda run + LD_LIBRARY_PATH + shell=True)

    Args:
        commands: 各段命令列表的列表,每段cmd[0]为工具完整路径|List of command lists
        work_dir: 工作目录|Working directory
    """
    if not commands:
        return False
    pipeline_str = _build_pipeline_command(commands)
    logger.info(f"命令|Command: {pipeline_str}")
    env = _conda_lib_env(_pipeline_tool_paths(commands))
    try:
        result = subprocess.run(
            pipeline_str, shell=True, cwd=work_dir, env=env,
            capture_output=True, text=True, check=True
        )
        if result.stderr:
            logger.debug(f"管道stderr|Pipeline stderr: {result.stderr}")
        logger.info("管道命令执行成功|Pipeline command executed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"管道命令执行失败|Pipeline command failed: {e}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False


def run_shell_command(shell_cmd: str, tool_paths: List[str], logger, work_dir: str = None) -> bool:
    """
    执行shell字符串命令(管道/重定向/复杂命令),设LD_LIBRARY_PATH(§13.2.1方案A)

    用于含awk脚本/重定向/tee等shell元字符的复杂命令;conda工具用完整路径,
    本函数据tool_paths设LD_LIBRARY_PATH实现env隔离(等同方案B剥离+方案A env)。
    |Run shell string command (pipe/redirect/complex), set LD_LIBRARY_PATH (§13.2.1 plan A).
    For commands with awk scripts/redirects/tee; conda tools use full path, this function
    sets LD_LIBRARY_PATH from tool_paths for env isolation.

    Args:
        shell_cmd: shell命令字符串(conda工具用完整路径)|shell command string
        tool_paths: 命令涉及的conda工具完整路径(用于LD_LIBRARY_PATH)|tool full paths
        work_dir: 工作目录|Working directory
    """
    logger.info(f"命令|Command: {shell_cmd}")
    env = _conda_lib_env(tool_paths)
    try:
        result = subprocess.run(
            shell_cmd, shell=True, cwd=work_dir, env=env,
            capture_output=True, text=True, check=True
        )
        if result.stderr:
            logger.debug(f"命令stderr|Command stderr: {result.stderr}")
        logger.info("命令执行成功|Command executed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command failed: {e}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False
    except FileNotFoundError as e:
        logger.error(f"命令不存在|Command not found: {e}")
        return False


def get_fasta_stats(fasta_file: str) -> dict:
    """获取FASTA文件统计信息|Get FASTA file statistics"""
    if not os.path.exists(fasta_file):
        return {}

    stats = {}

    # 计算序列数和总长度|Calculate sequence count and total length
    with open(fasta_file) as f:
        seq_lengths = []
        current_seq = []

        for line in f:
            if line.startswith('>'):
                if current_seq:
                    seq_lengths.append(len(''.join(current_seq)))
                    current_seq = []
            else:
                current_seq.append(line.strip())

        if current_seq:
            seq_lengths.append(len(''.join(current_seq)))

    stats['num_seqs'] = len(seq_lengths)
    stats['total_len'] = sum(seq_lengths)

    # 计算最长序列|Calculate longest sequence
    stats['max_len'] = max(seq_lengths) if seq_lengths else 0

    # 计算N50|Calculate N50
    if seq_lengths:
        sorted_lengths = sorted(seq_lengths, reverse=True)
        total_len = sum(sorted_lengths)
        half_len = total_len / 2

        cum_len = 0
        for length in sorted_lengths:
            cum_len += length
            if cum_len >= half_len:
                stats['n50'] = length
                break

    return stats


def format_time(seconds: int) -> str:
    """格式化时间|Format time"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours}小时 {minutes}分钟 {seconds}秒"


def generate_software_versions_yml(config, output_file: str, logger) -> None:
    """
    生成software_versions.yml文件(§12.5)|Generate software_versions.yml (§12.5)

    Args:
        config: 配置对象(需实现get_software_info)|Config object (must implement get_software_info)
        output_file: 输出文件路径|Output file path
        logger: 日志器|Logger
    """
    import yaml
    from datetime import datetime
    try:
        info = config.get_software_info()
        info['execution'] = {'timestamp': datetime.now().isoformat()}
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, sort_keys=False)
        logger.info(f"版本信息已保存|Version info saved: {output_file}")
    except Exception as e:
        logger.warning(f"生成版本信息文件失败|Failed to generate version info: {e}")


def generate_contig_reads_map_from_gfa(gfa_file: str, output_file: str, logger) -> bool:
    """
    从GFA文件生成contig-reads映射文件|Generate contig-reads mapping file from GFA

    解析GFA文件中的A行(read alignments)，生成contig到reads的映射文件
    Parse A lines (read alignments) in GFA file, generate contig to reads mapping

    Args:
        gfa_file: GFA文件路径|GFA file path
        output_file: 输出映射文件路径|Output mapping file path
        logger: 日志对象|Logger object

    Returns:
        bool: 是否成功|Whether successful
    """
    logger.info("=" * 60)
    logger.info("生成Contig-Reads映射文件|Generating Contig-Reads mapping file")
    logger.info("=" * 60)

    if not os.path.exists(gfa_file):
        logger.warning(f"GFA文件不存在|GFA file not found: {gfa_file}")
        return False

    try:
        # 第一步：解析GFA文件中的A行|Step 1: Parse A lines in GFA file
        logger.info(f"解析GFA文件|Parsing GFA file: {gfa_file}")

        contig_to_reads = {}
        read_count = 0

        with open(gfa_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('A\t'):
                    continue

                # A行格式: A\tcontig_id\tposition\tstrand\tread_name\tread_start\tread_end\ttags
                # A line format: A\tcontig_id\tposition\tstrand\tread_name\tread_start\tread_end\ttags
                parts = line.split('\t')
                if len(parts) >= 5:
                    contig_id = parts[1]
                    read_name = parts[4]

                    if contig_id not in contig_to_reads:
                        contig_to_reads[contig_id] = []

                    # 避免重复添加同一个read到同一个contig
                    # Avoid adding duplicate read to same contig
                    if read_name not in contig_to_reads[contig_id]:
                        contig_to_reads[contig_id].append(read_name)
                        read_count += 1

        logger.info(f"  共处理|Total alignments processed: {read_count:,}")
        logger.info(f"  共发现|Total contigs found: {len(contig_to_reads):,}")

        # 第二步：写入输出文件（按contig ID排序）
        # Step 2: Write to output file (sorted by contig ID)
        logger.info(f"写入映射文件|Writing mapping file: {output_file}")

        # 按数字部分排序|Sort by numeric part
        sorted_contigs = sorted(contig_to_reads.keys(),
                                key=lambda x: int(''.join(filter(str.isdigit, x))) if ''.join(filter(str.isdigit, x)) else 0)

        with open(output_file, 'w') as f_out:
            # 写入header|Write header
            f_out.write("#contig_id\tread_name\n")

            for contig_id in sorted_contigs:
                for read_name in contig_to_reads[contig_id]:
                    f_out.write(f"{contig_id}\t{read_name}\n")

        # 统计信息|Statistics
        logger.info("映射统计|Mapping statistics:")
        logger.info(f"  总contig数|Total contigs: {len(contig_to_reads):,}")

        # 计算每个contig的reads数量|Calculate reads per contig
        reads_per_contig = [len(reads) for reads in contig_to_reads.values()]
        if reads_per_contig:
            logger.info(f"  平均reads数/contig|Avg reads/contig: {sum(reads_per_contig) / len(reads_per_contig):.1f}")
            logger.info(f"  最少reads数|Min reads: {min(reads_per_contig):,}")
            logger.info(f"  最多reads数|Max reads: {max(reads_per_contig):,}")

        logger.info(f"  输出文件|Output file: {output_file}")
        logger.info("=" * 60)
        logger.info("Contig-Reads映射文件生成完成|Contig-Reads mapping file generated successfully")
        logger.info("=" * 60)

        return True

    except Exception as e:
        logger.error(f"生成contig-reads映射失败|Failed to generate contig-reads mapping: {str(e)}")
        return False
