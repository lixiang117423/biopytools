# Conda 环境软件调用详解|Conda Invocation Details

> 📖 本文件为 CLAUDE.md 的按需参考文档,**调用 conda 环境软件 / 排查 conda 命令时读取**。核心铁律见 CLAUDE.md §十三。
> On-demand reference for CLAUDE.md; read when invoking software inside conda envs or debugging conda commands.

---

## 1. 问题背景|Problem Background

conda环境中的软件（特别是Python包）无法通过直接调用执行，因为：

1. **依赖隔离|Dependency Isolation**: conda环境的Python依赖包仅在该环境中可用
2. **路径问题|Path Issues**: `which busco`可能返回非conda环境的路径
3. **错误示例|Error Example**:
   ```bash
   # 直接调用失败|Direct call fails
   $ /miniforge3/envs/BUSCO_v.6.0.0/bin/busco --version
   No module named 'busco'

   # 使用conda run成功|Using conda run succeeds
   $ conda run -n BUSCO_v.6.0.0 busco --version
   BUSCO 6.0.0
   ```

---

## 2. 管道中禁用 conda run 的解决方案|Pipeline Solutions

> 核心禁令见 CLAUDE.md §13.2.1：严禁 `conda run | conda run`。以下为三种替代方案。

### 方案A：管道命令直接调用（设置LD_LIBRARY_PATH）

```python
# 管道中的工具直接调用，设置环境变量
import os
from common.paths import expand_path

env = os.environ.copy()
conda_env_lib = expand_path("~/miniforge3/envs/your_env/lib")
env['LD_LIBRARY_PATH'] = f"{conda_env_lib}:{env.get('LD_LIBRARY_PATH', '')}"

# 直接调用，不使用conda run
proc1 = subprocess.Popen([bwa_bin] + bwa_args, stdout=subprocess.PIPE, env=env)
proc2 = subprocess.Popen([samtools_bin] + samtools_args, stdin=proc1.stdout, env=env)
```

### 方案B：使用_build_pipeline_command辅助函数

```python
# 使用辅助函数处理管道命令（自动提取conda run中的实际命令）
def _build_pipeline_command(self, commands: List[List[str]]) -> str:
    """构建conda环境包装的管道命令"""
    wrapped_commands = []
    for cmd in commands:
        wrapped_cmd = build_conda_command(cmd[0], cmd[1:])
        # 提取conda run中的实际命令，避免 "conda run | conda run"
        if wrapped_cmd[0] == 'conda' and len(wrapped_cmd) > 4:
            actual_cmd = ' '.join(wrapped_cmd[4:])  # 跳过 'conda', 'run', '-n', env_name
            wrapped_commands.append(actual_cmd)
        else:
            wrapped_commands.append(' '.join(wrapped_cmd))
    return ' | '.join(wrapped_commands)

# 使用示例
bwa_cmd = [self.config.bwa_bin, "mem", "-t", "16", asm.fa]
samtools_cmd = [self.samtools_path, "view", "-b", "-o", "output.bam"]
pipeline = self._build_pipeline_command([bwa_cmd, samtools_cmd])
```

### 方案C：管道中的工具使用shutil.which自动检测

```python
# 对于管道中必须使用的工具，用shutil.which自动检测
import shutil

def __init__(self, config, logger):
    # samtools使用系统自动检测，不使用conda run
    self.samtools_path = shutil.which('samtools')
    if not self.samtools_path:
        raise RuntimeError("samtools not found in PATH")

# 管道中使用
samtools_cmd = [self.samtools_path, "view", "-b", "-o", "output.bam"]
```

---

## 3. 完整实现|Full Implementation

> 📌 **2026-08 更新|Update:** 以下实现的权威版本已下沉到公共层
> **biopytools/common/conda_runner.py**（get_conda_env / build_conda_command /
> build_pipeline_command / run_pipeline / CommandRunner / check_tools）。
> 新代码直接 "from ..common.conda_runner import ..." 导入，**禁止在各模块
> utils.py 中复制实现**；此处保留代码仅供理解算法。
> The canonical implementations now live in biopytools/common/conda_runner.py;
> import from there instead of copying per module.

### 3.1 环境检测函数|Environment Detection Function

```python
import shutil
import re
import os
from typing import Optional

def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称

    策略|Strategy:
    1. 首先尝试从which命令路径检测（优先级高）
    2. 如果未找到，搜索所有conda环境（兜底方案）

    Args:
        command: 命令名称或路径 (e.g., 'busco' or '/path/to/busco')

    Returns:
        conda环境名称或None (e.g., 'BUSCO_v.6.0.0' or None)
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs
        # 例如: /miniforge3/envs/BUSCO_v.6.0.0/bin/busco
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None
```

### 3.2 命令构建函数|Command Building Function

```python
from typing import List

def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表 (适用于subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('busco', ['--version'])
        ['conda', 'run', '-n', 'BUSCO_v.6.0.0', '--no-capture-output', 'busco', '--version']

        >>> # 绝对路径且不在conda envs目录下时，直接调用
        >>> # Absolute path not under conda envs: called directly
        >>> build_conda_command('/usr/bin/tool', ['--help'])
        ['/usr/bin/tool', '--help']

    注意|Note:
        返回的列表应配合 subprocess.run(shell=False) 使用
        The returned list must be used with subprocess.run(shell=False)

    ⚠️ 重要|IMPORTANT:
        必须使用--no-capture-output避免conda缓冲输出导致内存问题
        Must use --no-capture-output to avoid conda buffering output causing memory issues
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用
        # 如果command是命令名，conda run会自动找到环境中的版本
        # 添加--no-capture-output避免缓冲输出导致内存问题
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用
        full_cmd = [command] + args

    return full_cmd
```

### 3.3 依赖检查函数|Dependency Check Function

```python
import subprocess
import logging

def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查依赖软件"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        # 使用conda wrapper构建命令
        cmd = build_conda_command(config.busco_path, ['--version'])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode == 0:
            version_info = result.stdout.strip()
            logger.info(f"BUSCO 可用|BUSCO available: {version_info}")
        else:
            raise RuntimeError("BUSCO版本检查失败|BUSCO version check failed")

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"BUSCO不可用|BUSCO not available: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True
```

### 3.4 命令执行器|Command Runner

```python
class CommandRunner:
    """命令执行器|Command Runner"""

    def run(self, cmd: list, description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        注意|Note:
            - 始终传入列表，使用shell=False（更安全）|Always pass a list, use shell=False (safer)
            - 禁止使用shell=True传入列表，仅执行第一个元素|Never use shell=True with list, only first element executes
            - 若必须使用shell执行字符串，请先: cmd_str = " ".join(cmd)|
              If shell string is needed: cmd_str = " ".join(cmd)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )
            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            return False, e.stdout, e.stderr
```

### 3.5 使用示例|Usage Example

```python
# 在分析模块中应用|Apply in analysis module
from .utils import build_conda_command

class BUSCORunner:
    def build_busco_command(self, input_file: str, output_name: str) -> str:
        """构建BUSCO命令"""
        # 构建参数列表
        args = [
            "-i", input_file,
            "-l", self.config.lineage,
            "-o", output_name,
            "-m", self.config.mode,
            "-c", str(self.config.threads)
        ]

        if self.config.offline:
            args.append("--offline")

        # 使用conda wrapper
        cmd_list = build_conda_command(self.config.busco_path, args)

        # 转换为字符串用于shell执行
        return " ".join(cmd_list)
```

---

## 4. 常见错误与禁令|Common Mistakes and Prohibitions

### 4.1 ❌ 禁止提取命令名|Prohibit Command Name Extraction

**错误做法|WRONG:**
```python
# ❌ 错误：提取命令名，丢失路径信息
cmd_name = os.path.basename(cmd[0])  # 错误！丢失了 /envs/ 路径
wrapped_cmd = build_conda_command(cmd_name, cmd[1:])

# 问题：
# 1. cmd[0] = "/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools"
# 2. cmd_name = "samtools" (路径信息丢失！)
# 3. get_conda_env("samtools") 无法从命令名中提取环境名
# 4. 最终生成错误命令: "conda run samtools ..." (缺少 -n <env>)
```

**正确做法|CORRECT:**
```python
# ✅ 正确：传递完整路径
wrapped_cmd = build_conda_command(cmd[0], cmd[1:])

# 结果：
# 1. cmd[0] = "/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools"
# 2. get_conda_env() 从路径提取环境名 "GATK_v.4.6.2.0"
# 3. 生成正确命令: "conda run -n GATK_v.4.6.2.0 samtools ..."
```

### 4.2 ❌ 禁止使用命令名调用|Prohibit Using Command Name Only

**错误代码示例|Error Example:**
```python
# ❌ 错误：只传递命令名
haphic_cmd_name = os.path.basename(self.config.haphic_bin)
wrapped_cmd = build_conda_command(haphic_cmd_name, ["--help"])
# get_conda_env("haphic") → None (可能找到系统PATH中的版本)
# 生成错误命令: conda run haphic --help
```

**正确代码示例|Correct Example:**
```python
# ✅ 正确：传递完整路径
wrapped_cmd = build_conda_command(self.config.haphic_bin, ["--help"])
# get_conda_env("/miniforge3/envs/haphic/bin/haphic") → "haphic"
# 生成正确命令: conda run -n haphic haphic --help
```

---

## 5. 测试验证|Testing

```python
# 单元测试示例|Unit test example
def test_conda_detection():
    """测试conda环境检测"""
    from biopytools.busco.utils import get_conda_env, build_conda_command

    # 测试conda环境中的命令（使用完整路径）
    env = get_conda_env('/miniforge3/envs/BUSCO_v.6.0.0/bin/busco')
    assert env == 'BUSCO_v.6.0.0'

    # 测试命令构建（使用完整路径）
    cmd = build_conda_command('/miniforge3/envs/BUSCO_v.6.0.0/bin/busco', ['--version'])
    assert cmd == ['conda', 'run', '-n', 'BUSCO_v.6.0.0', '--no-capture-output', 'busco', '--version']

    # 测试非conda命令
    env = get_conda_env('/usr/bin/bash')
    assert env is None

    # ❌ 测试错误用法：只传命令名可能失败
    env = get_conda_env('busco')  # 可能返回 None（取决于系统PATH）
    # 注意：不应依赖命令名，始终使用完整路径
```

---

## 6. 故障排查|Troubleshooting

**症状|Symptom:** 生成的命令缺少 `-n <env_name>` 参数

```bash
# 错误命令
conda run samtools view -b -@ 64 -o output.bam -

# 正确命令应该是
conda run -n GATK_v.4.6.2.0 samtools view -b -@ 64 -o output.bam -
```

**排查步骤|Debugging Steps:**

1. **检查调用点|Check Call Site:**
   ```python
   # 添加调试日志
   logger.debug(f"build_conda_command input: command={command}")
   wrapped_cmd = build_conda_command(command, args)
   logger.debug(f"build_conda_command output: {wrapped_cmd}")
   ```

2. **检查command参数|Check command Parameter:**
   - 如果是 `"samtools"` → 错误！使用了命令名
   - 如果是 `"/miniforge3/envs/.../bin/samtools"` → 正确！使用了完整路径

3. **检查get_conda_env返回值|Check get_conda_env Return:**
   ```python
   conda_env = get_conda_env(command)
   logger.debug(f"Detected conda_env: {conda_env}")
   # 如果是 None，说明无法检测环境，通常是路径问题
   ```

4. **修复方法|Fix:**
   - 找到传递 `os.path.basename(...)` 的代码
   - 改为直接传递完整路径
   - 移除所有 `cmd_name = os.path.basename(...)` 相关代码
