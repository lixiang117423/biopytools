"""Conda环境软件调用公共模块|Conda Environment Software Invocation Common Module

统一的外部命令调用层|Unified external command invocation layer:
- get_conda_env / build_conda_command: conda环境自动检测与包装(§13.2)
- build_pipeline_command / run_pipeline: 管道命令(§13.2.2, 方案B:
  提取conda run中的实际命令, 严禁 conda run | conda run)
- CommandRunner: 标准命令执行器(记录完整命令到INFO, §2.2.1)
- check_tools: 统一依赖检查

使用方式|Usage:
    from biopytools.common.conda_runner import (
        build_conda_command, CommandRunner, check_tools
    )

    cmd = build_conda_command(config.samtools_path, ["index", bam_file])
    success, stdout, stderr = runner.run(cmd, "构建BAM索引|Building BAM index")
"""

import logging
import os
import re
import shutil
import subprocess
from typing import List, Optional, Sequence, Tuple, Union


def _conda_envs_dir() -> Optional[str]:
    """conda环境根目录|Conda envs root directory

    优先用 CONDA_EXE 定位, 回退 ~/miniforge3/envs
    |Locate via CONDA_EXE first, fall back to ~/miniforge3/envs
    """
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        return os.path.join(os.path.dirname(os.path.dirname(conda_exe)), 'envs')
    fallback = os.path.expanduser('~/miniforge3/envs')
    if os.path.isdir(fallback):
        return fallback
    return None


def get_conda_env(command: str, preferred: Optional[str] = None) -> Optional[str]:
    """检测命令所在的conda环境名称|Detect conda env name where the command resides

    检测顺序|Detection order:
    1. preferred 指定环境(存在才用)|preferred env (only if it exists)
    2. 命令为完整路径且含 /envs/<name>/ → 直接提取|full path under envs → extract
    3. shutil.which 路径含 /envs/<name>/ → 提取|which path under envs → extract
    4. 功能域映射表(common/env_map.py) 中该工具所属域环境存在 → 使用
       |domain env from env_map (only if the binary exists there)
    5. 扫描所有conda环境的 bin/<basename>|scan all envs for bin/<basename>

    Args:
        command: 命令名称或完整路径|Command name or full path
                 (e.g., 'samtools' or '~/miniforge3/envs/align/bin/samtools')
        preferred: 优先使用的环境名|Preferred env name (optional)

    Returns:
        conda环境名称或None|conda env name or None
    """
    basename = os.path.basename(command)
    envs_dir = _conda_envs_dir()

    # 1. preferred 环境优先|Preferred env first
    if preferred:
        if envs_dir and os.path.exists(os.path.join(envs_dir, preferred, 'bin', basename)):
            return preferred
        if os.path.isabs(command) and f'/envs/{preferred}/' in command:
            return preferred

    # 2. 完整路径直接提取|Extract directly from full path
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)/', command)
        if match:
            return match.group(1)
        # 非 conda envs 下的绝对路径 → 不包装(如 ~/software/xxx)
        # |Absolute path outside envs → no wrapping (e.g. ~/software/xxx)
        return None

    # 3. which 路径提取|Extract from which path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)/', cmd_path)
        if match:
            return match.group(1)

    # 4. 功能域映射表|Domain env from env_map
    if envs_dir:
        try:
            from .env_map import TOOL_DOMAIN_MAP
            domain = TOOL_DOMAIN_MAP.get(basename)
            if domain and os.path.exists(os.path.join(envs_dir, domain, 'bin', basename)):
                return domain
        except Exception:
            pass

        # 5. 全量扫描|Scan all envs
        try:
            for env_name in sorted(os.listdir(envs_dir)):
                if os.path.exists(os.path.join(envs_dir, env_name, 'bin', basename)):
                    return env_name
        except OSError:
            pass

    return None


def build_conda_command(
    command: str,
    args: Sequence[str],
    preferred_env: Optional[str] = None,
) -> List[str]:
    """构建conda run命令(单工具)|Build conda run command (single tool)

    自动检测conda环境并包装, 强制 --no-capture-output(§13.2.1);
    非conda软件直接调用, 向后兼容|Auto-detect conda env and wrap with
    --no-capture-output; non-conda tools are called directly (backward compatible)

    Args:
        command: 命令名称或完整路径|Command name or full path
                 (必须传完整路径, 禁止 os.path.basename 提取命令名, §13.2.3)
        args: 命令参数列表|Command argument list
        preferred_env: 优先使用的环境名|Preferred env name (optional)

    Returns:
        完整命令列表(配合 subprocess.run(shell=False) 使用)
        |Complete command list (use with subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('~/miniforge3/envs/align/bin/samtools', ['--version'])
        ['conda', 'run', '-n', 'align', '--no-capture-output',
         '~/miniforge3/envs/align/bin/samtools', '--version']
    """
    conda_env = get_conda_env(command, preferred=preferred_env)

    if conda_env:
        # conda环境软件用conda run包装|Wrap conda env tools with conda run
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + list(args)

    # 非conda环境, 直接调用|Non-conda environment, direct call
    return [command] + list(args)


def _extract_actual_command(wrapped_cmd: List[str]) -> List[str]:
    """从conda run包装命令中提取实际命令|Extract actual command from conda run wrapper

    方案B(§13.2.2): 管道中禁止 conda run, 提取实际命令后直接调用
    |Solution B: conda run is forbidden in pipelines; extract and call directly

    Args:
        wrapped_cmd: build_conda_command 的返回值|Result of build_conda_command

    Returns:
        实际命令列表|Actual command list
    """
    if wrapped_cmd and wrapped_cmd[0] == 'conda':
        # ['conda', 'run', '-n', env, '--no-capture-output', command, ...]
        idx = 4
        if len(wrapped_cmd) > idx and wrapped_cmd[idx] == '--no-capture-output':
            idx += 1
        return wrapped_cmd[idx:]
    return list(wrapped_cmd)


def build_pipeline_command(segments: Sequence[Sequence[str]]) -> str:
    """构建管道命令字符串(方案B)|Build pipeline command string (solution B)

    每个管道段先经 build_conda_command 包装, 再提取实际命令,
    避免 "conda run | conda run"(§13.2.2)
    |Each segment is wrapped then unwrapped to avoid conda run in pipes

    Args:
        segments: 管道段列表|List of pipeline segments, 如
                  [['bwa', 'mem', 'ref.fa'], ['samtools', 'view', '-b']]

    Returns:
        管道命令字符串(供 shell=True 使用)|Pipeline command string (for shell=True)
    """
    parts = []
    for seg in segments:
        wrapped = build_conda_command(seg[0], list(seg[1:]))
        parts.append(' '.join(_extract_actual_command(wrapped)))
    return ' | '.join(parts)


def run_pipeline(
    segments: Sequence[Sequence[str]],
    logger: logging.Logger,
    description: str = "",
    output_file: Optional[str] = None,
    working_dir: Optional[str] = None,
    binary_output: bool = False,
) -> Tuple[bool, str]:
    """执行管道命令(方案B, Popen链接)|Run pipeline command (solution B, Popen chaining)

    与 build_pipeline_command 相同的提取逻辑, 但用 Popen 链式执行,
    避免 shell=True; 每个管道段的工具按实际命令(域环境二进制)直接调用
    |Same extraction as build_pipeline_command but chains Popen processes
    instead of shell=True; each segment runs its actual command directly

    Args:
        segments: 管道段列表|List of pipeline segments
        logger: 日志器|Logger
        description: 步骤描述|Step description
        output_file: 末段stdout重定向文件(可空, 为空时末段stdout丢弃)
                     |Redirect last segment stdout to file (optional)
        working_dir: 工作目录|Working directory
        binary_output: output_file 是否为二进制(如BAM)|Whether output_file is binary

    Returns:
        (是否成功, 错误信息)|(success, error message)
    """
    actual_segments = []
    for seg in segments:
        wrapped = build_conda_command(seg[0], list(seg[1:]))
        actual_segments.append(_extract_actual_command(wrapped))

    if description:
        logger.info(f"执行步骤|Executing step: {description}")
    display = ' | '.join(' '.join(seg) for seg in actual_segments)
    logger.info(f"命令|Command: {display}")

    procs: List[subprocess.Popen] = []
    out_fh = None
    try:
        for i, seg in enumerate(actual_segments):
            last = (i == len(actual_segments) - 1)
            stdin = procs[-1].stdout if procs else None
            if last:
                if output_file:
                    mode = 'wb' if binary_output else 'w'
                    out_fh = open(output_file, mode)
                    stdout = out_fh
                else:
                    stdout = subprocess.DEVNULL
            else:
                stdout = subprocess.PIPE

            proc = subprocess.Popen(
                seg,
                stdin=stdin,
                stdout=stdout,
                stderr=subprocess.PIPE,
                cwd=working_dir,
            )
            # 关闭上一进程的stdout, 让上一进程在管道破裂时收到SIGPIPE
            # |Close previous stdout so upstream gets SIGPIPE on broken pipe
            if stdin is not None:
                procs[-1].stdout.close()
            procs.append(proc)

        stderr_tail = ""
        failed = False
        for proc in procs:
            err = (proc.stderr.read() or b"")
            returncode = proc.wait()
            if returncode != 0:
                failed = True
                decoded = err.decode('utf-8', errors='replace')
                if decoded:
                    stderr_tail = decoded
        if out_fh is not None:
            out_fh.close()
            out_fh = None

        if failed:
            logger.error(f"管道执行失败|Pipeline failed: {description}")
            if stderr_tail:
                logger.error(f"错误信息|Error message: {stderr_tail[:500]}")
            return False, stderr_tail

        return True, ""

    except OSError as e:
        if out_fh is not None:
            out_fh.close()
        logger.error(f"管道执行异常|Pipeline execution error: {description}")
        logger.error(f"错误信息|Error message: {e}")
        return False, str(e)


class CommandRunner:
    """命令执行器|Command Runner

    统一执行外部命令: 记录完整命令到INFO级别(§2.2.1),
    支持列表(shell=False)与字符串(shell=True, 已由调用方自行构建)
    |Execute external commands uniformly: log the full command at INFO,
    accept both lists (shell=False) and pre-built strings (shell=True)
    """

    def __init__(self, logger: logging.Logger, working_dir: Optional[str] = None):
        """
        Args:
            logger: 日志器|Logger
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = working_dir

    def run(
        self,
        cmd: Union[List[str], str],
        description: str = "",
        capture_output: bool = True,
        output_file: Optional[str] = None,
        timeout: Optional[int] = None,
        env: Optional[dict] = None,
    ) -> Tuple[bool, str, str]:
        """执行命令|Execute command

        Args:
            cmd: 命令列表(由 build_conda_command 构建)或命令字符串
                 |Command list (built by build_conda_command) or shell string
            description: 步骤描述|Step description
            capture_output: 是否捕获stdout/stderr(未指定output_file时生效)
                            |Whether to capture output (when output_file is None)
            output_file: stdout重定向文件(大数据量输出用, 避免内存缓冲)
                         |Redirect stdout to file (for large outputs)
            timeout: 超时秒数|Timeout in seconds
            env: 环境变量字典|Environment variables

        Returns:
            (是否成功, stdout, stderr)|(success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        if isinstance(cmd, (list, tuple)):
            use_shell = False
            display = ' '.join(str(c) for c in cmd)
        else:
            use_shell = True
            display = str(cmd)
        self.logger.info(f"命令|Command: {display}")

        out_fh = None
        try:
            kwargs = dict(shell=use_shell, text=True, cwd=self.working_dir, env=env)
            if output_file is not None:
                out_fh = open(output_file, 'w')
                kwargs['stdout'] = out_fh
                kwargs['stderr'] = subprocess.PIPE
            else:
                kwargs['capture_output'] = capture_output
                if not capture_output:
                    kwargs['stderr'] = subprocess.PIPE
            if timeout is not None:
                kwargs['timeout'] = timeout

            result = subprocess.run(cmd, **kwargs)

            if out_fh is not None:
                out_fh.close()
                out_fh = None

            stdout = result.stdout or ""
            stderr = result.stderr or ""

            if result.returncode != 0:
                self.logger.error(
                    f"命令执行失败|Command execution failed: {description}"
                )
                self.logger.error(f"错误代码|Error code: {result.returncode}")
                if stderr:
                    self.logger.error(f"错误信息|Error message: {stderr[:500]}")
                return False, stdout, stderr

            return True, stdout, stderr

        except subprocess.TimeoutExpired:
            if out_fh is not None:
                out_fh.close()
            self.logger.error(f"命令超时|Command timed out: {description}")
            return False, "", "Timeout"

        except OSError as e:
            if out_fh is not None:
                out_fh.close()
            self.logger.error(f"命令执行异常|Command execution error: {description}")
            self.logger.error(f"错误信息|Error message: {e}")
            return False, "", str(e)


def check_tools(
    dependencies: Sequence[Sequence],
    logger: logging.Logger,
    timeout: int = 30,
) -> List[str]:
    """统一依赖检查|Unified dependency check

    Args:
        dependencies: 依赖列表|Dependency list, 每项为:
                      (tool_path, 显示名|display name, 版本参数|version args)
                      或 (tool_path, 显示名, 版本参数, 非零返回码是否视为可用)
                      |or 4-tuple with accept_nonzero flag
        logger: 日志器|Logger
        timeout: 单工具检查超时秒数|Per-tool timeout in seconds

    Returns:
        缺失工具显示名列表|List of missing tool display names
    """
    logger.info("检查依赖软件|Checking dependencies")
    missing = []

    for dep in dependencies:
        path = dep[0]
        name = dep[1]
        version_args = list(dep[2])
        accept_nonzero = bool(dep[3]) if len(dep) > 3 else False

        cmd = build_conda_command(path, version_args)
        logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing.append(name)
            logger.error(f"{name} 不可用|{name} not available")
            continue

        if result.returncode == 0 or accept_nonzero:
            info = ((result.stdout or '') + '\n' + (result.stderr or '')).strip()
            version = info.split('\n')[0] if info else 'installed'
            logger.info(f"{name} 可用|{name} available: {version}")
        else:
            missing.append(name)
            logger.error(f"{name} 不可用|{name} not available")

    return missing


__all__ = [
    'get_conda_env',
    'build_conda_command',
    'build_pipeline_command',
    'run_pipeline',
    'CommandRunner',
    'check_tools',
]
