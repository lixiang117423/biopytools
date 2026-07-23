# BioPyTools Python代码开发规范文档

## 版本: 2.15
## 日期: 2026-07-23
## 用途: 统一所有生信分析模块的代码结构、命名规范、日志格式

## 目录|Table of Contents

- [重要警告|CRITICAL WARNING](#-重要警告critical-warning)
- [禁止在超算上Git提交|No Git Commits on Supercomputer](#-禁止在超算上执行-git-提交forbidden-git-commits-on-the-supercomputer)
- [一、模块化结构规范](#一模块化结构规范)
- [二、日志格式规范](#二日志格式规范)
- [三、命令行参数规范](#三命令行参数规范)
- [四、命名规范](#四命名规范)
- [五、代码风格规范](#五代码风格规范)
- [六、错误处理规范](#六错误处理规范)
- [七、代码改造检查清单](#七代码改造检查清单)
- [八、标准代码模板](#八标准代码模板)
- [九、常见问题](#九常见问题)
- [十、开发约束与工作流规范](#十开发约束与工作流规范)
- [十一、测试规范与路径管理规范](#十一测试规范与路径管理规范)
- [十二、输出目录和文件命名规范](#十二输出目录和文件命名规范)
- [十三、Conda环境软件调用规范](#十三conda环境软件调用规范)
- [版本历史](#版本历史)

---

## ⚠️ 重要警告|CRITICAL WARNING

### 🚫 禁止硬编码绝对路径|FORBIDDEN: Hardcoded Absolute Paths

**所有新代码和代码修改必须遵守路径管理规范|All new code and modifications MUST follow path management standards**

❌ **严格禁止|STRICTLY FORBIDDEN:**
```python
# 任何包含用户名的绝对路径|Any absolute path with username
tool_path = "/share/org/YZWL/yzwl_lixg/miniforge3/bin/tool"
software_dir = "/share/org/YZWL/yzwl_lixg/software/module"
```

✅ **必须使用|MUST USE:**
```python
# 使用~展开或相对路径|Use ~ expansion or relative paths
tool_path = "~/miniforge3/bin/tool"  # ✅ 可接受|Acceptable
software_dir = "~/software/module"    # ✅ 可接受|Acceptable

# 更好的做法|Better: 使用路径管理系统|Use path management system
from common.paths import get_tool_path, expand_path
tool_path = get_tool_path('tool', '~/miniforge3/bin/tool', 'TOOL_PATH')
```

**后果|Consequences:**
- 🔴 代码审查不通过|Code review will be rejected
- 🔴 无法合并到主分支|Will not be merged to main branch
- 🔴 破坏代码可移植性|Breaks code portability

**参考|Reference:** 详见第十一章|See Section 11 for details

---

## 🚫 禁止在超算上执行 Git 提交|FORBIDDEN: Git Commits on the Supercomputer

> **⚠️ AI 在超算（登录节点/计算节点）上开发时，严禁执行任何 Git 写操作**
> **When developing on the supercomputer, AI MUST NOT perform ANY Git write operations**

❌ **超算上严禁执行|STRICTLY FORBIDDEN on the supercomputer:**
```bash
git add ...        # ❌ 禁止|Forbidden
git commit ...     # ❌ 禁止|Forbidden
git push ...       # ❌ 禁止|Forbidden
git tag ...        # ❌ 禁止|Forbidden
```

✅ **超算上允许执行|ALLOWED on the supercomputer:**
```bash
git status         # ✅ 可查看改动|View changes
git diff           # ✅ 可查看差异|View diff
# 编辑代码、运行测试、提交计算任务均可|Edit code, run tests, submit jobs are all fine
```

**原因|Reason:**
1. 超算出网受限，`git push` 几乎总是失败|Supercomputer has restricted egress; push almost always fails
2. 超算上的 `.git` **不会同步到本地**（`copybiopytools` 已排除 `.git/`），在超算 commit 会产生与本地/GitHub 分叉的历史，需复杂的 rebase 整理|.git on the supercomputer is never synced (copybiopytools excludes .git/); committing there creates divergent history requiring complex rebase
3. **代码提交的唯一入口是本地 Mac**|The only entry point for commits is the local Mac

**正确工作流|Correct Workflow:**
1. **超算**：只编辑代码，**不碰 git 写操作**|Supercomputer: edit code only, no git write operations
2. **同步**：本地执行 `copybiopytools` 把超算工作区同步到本地（仅同步代码文件，不含 `.git/`）|Sync: run copybiopytools locally to pull the supercomputer worktree (code files only, no .git/)
3. **本地 Mac**：由 Claude 检查改动 → commit（遵循第 10.1 节 message 规范）→ push 到 GitHub|Mac: Claude reviews changes → commits (per §10.1 message spec) → pushes to GitHub

> 📌 一句话：**超算只写代码，Mac 才提交。|The supercomputer writes code only; the Mac is the only place commits happen.**

---

## 一、模块化结构规范

### 1.1 标准模块目录结构

每个功能模块应采用以下目录结构：

```
biopytools/
├── module_name/
│   ├── __init__.py          # 模块导出声明
│   ├── config.py            # 配置数据类
│   ├── utils.py             # 工具函数(日志、命令执行等)
│   ├── calculator.py        # 核心计算逻辑,命名不一定叫calculator.py
│   └── main.py              # 命令行入口
```

### 1.2 各文件职责说明

#### __init__.py
- 定义模块导出的类和函数
- 声明模块版本号

```python
"""模块功能中文名|Module Function English Name"""

from .main import ModuleRunner
from .config import ModuleConfig
from .calculator import ModuleCalculator

__version__ = "1.0.0"

__all__ = [
    'ModuleRunner',
    'ModuleConfig',
    'ModuleCalculator'
]
```

#### config.py
- 使用 `@dataclass` 定义配置类
- 包含所有配置参数和验证逻辑

```python
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

@dataclass
class ModuleConfig:
    """配置类中文名|Config Class English Name"""

    # 必需参数|Required parameters
    input_file: str
    output_dir: str

    # 可选参数|Optional parameters
    threads: int = 24
    kmer_size: Optional[int] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not Path(self.input_file).exists():
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if errors:
            raise ValueError("\n".join(errors))
```

#### utils.py
- 日志管理器类
- 通用工具函数
- 命令执行函数

```python
import logging
import sys

class ModuleLogger:
    """模块日志管理器|Module Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
```

#### calculator.py
- 核心分析逻辑
- 继承或组合配置类和日志器
- 不一定叫calculator.py，根据模块名称灵活调整

```python
class ModuleCalculator:
    """模块计算器|Module Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run_analysis(self):
        """运行分析|Run analysis"""
        self.logger.info("开始分析|Starting analysis")
        # 核心逻辑
        results = {}
        return results
```

#### main.py
- 命令行参数解析
- 程序入口函数
- 使用 argparse 或 click

```python
import argparse
import sys
from pathlib import Path

def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="工具功能描述中英文|Tool Description EN"
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入文件中文名|Input file EN')
    parser.add_argument('-o', '--output-dir', default='./output',
                        help='输出目录中文名|Output directory EN')

    return parser.parse_args()

def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = ModuleConfig(
            input_file=args.input,
            output_dir=args.output_dir
        )
        config.validate()

        logger_manager = ModuleLogger()
        logger = logger_manager.get_logger()

        # 运行分析
        logger.info("分析完成|Analysis completed")
        sys.exit(0)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
```

---

## 二、日志格式规范

### 2.1 标准日志格式

**格式定义：**
```
YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 消息中文|Message English
```

**具体要求：**
- 时间戳格式：`年-月-日 时:分:秒.毫秒`
- 使用点号(.)分隔秒和毫秒（不是逗号）
- 毫秒固定3位，不足补0
- 使用 ` - ` (空格-空格) 分隔各部分（不是方括号[]）
- 使用 `|` 分隔中英文

**Python配置：**
```python
log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
date_format = '%Y-%m-%d %H:%M:%S'
```

**实际输出示例：**
```
2026-01-04 13:11:00.587 - INFO - 开始分析|Starting analysis
2026-01-04 13:11:01.234 - INFO - 输入文件|Input file: sample.fastq
2026-01-04 13:11:02.001 - WARNING - 质量较低|Low quality detected
2026-01-04 13:11:03.456 - ERROR - 文件不存在|File not found
```

### 2.2 日志级别使用规范

| 级别 | 用途 | 示例 |
|------|------|------|
| **DEBUG** | 详细调试信息 | 变量值、中间结果 |
| **INFO** | 正常流程信息 | 步骤开始/完成、统计数据 |
| **WARNING** | 警告信息 | 质控警告、参数建议 |
| **ERROR** | 错误信息 | 文件读取失败、格式错误 |
| **CRITICAL** | 严重错误 | 程序无法继续的错误 |

#### 2.2.1 命令执行的日志规范|Command Execution Logging Standards

⚠️ **重要|IMPORTANT**: 所有外部命令执行时，**必须记录完整命令到 INFO 级别**

**原因|Reason:**
1. **调试需求** - 用户需要知道实际执行了什么命令
2. **论文写作** - Methods部分需要记录完整的软件版本和参数
3. **可重复性** - 其他研究者需要能够完全复现分析
4. **黑箱问题** - 避免工具成为"黑箱"，提升透明度

**错误示例|WRONG Example:**
```python
# ❌ 只记录描述，不记录命令
logger.info("执行BWA索引")
# 执行命令
subprocess.run([bwa_bin, 'index', genome.fa])
```

**正确示例|CORRECT Example:**
```python
# ✅ 先记录描述，再记录完整命令
logger.info(f"执行|Executing: BWA索引|BWA indexing")
# 记录完整命令（INFO级别，不是DEBUG）
logger.info(f"命令|Command: {' '.join(cmd)}")
# 执行命令
subprocess.run(cmd)
```

**日志输出示例|Log Output Example:**
```
2026-03-17 20:11:31.311 - INFO - 执行|Executing: BWA索引|BWA indexing
2026-03-17 20:11:31.312 - INFO - 命令|Command: conda run -n bwa_v0.7.17 bwa index /path/to/genome.fa -p bwt
```

**管道命令示例|Pipeline Command Example:**
```python
# ✅ 管道命令也要完整记录
logger.info(f"命令|Command: {' '.join(bwa_cmd)} | {' '.join(samtools_cmd)} > {output_file}")
```

**检查清单|Checklist:**
- [ ] 命令执行前记录了完整命令（INFO级别）
- [ ] 包含了conda run包装（如果使用）
- [ ] 包含了所有参数和文件路径
- [ ] 管道命令用 `|` 连接显示
- [ ] 格式清晰，便于复制到论文Methods部分

### 2.3 超算日志系统分离规范|Job Scheduler Log Separation

在超算系统（如使用 `sub` 函数提交任务）运行时，必须正确配置 stdout 和 stderr handler，以实现日志分离：

> **原则|Principle:** INFO → stdout → .out 文件，WARNING/ERROR → stderr → .err 文件

#### 2.3.1 标准配置|Standard Configuration

```python
import logging
import sys

logger = logging.getLogger("ModuleName")
logger.setLevel(logging.DEBUG)
logger.handlers.clear()
logger.propagate = False  # 避免重复|Avoid duplicates

formatter = logging.Formatter(
    '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# stdout handler - INFO级别|stdout handler - INFO level
# → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
stdout_handler.setFormatter(formatter)
logger.addHandler(stdout_handler)

# stderr handler - WARNING及以上|stderr handler - WARNING and above
# → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setLevel(logging.WARNING)
stderr_handler.setFormatter(formatter)
logger.addHandler(stderr_handler)

# 文件handler - 所有级别|File handler - all levels
# → 本地完整日志|→ Local complete log
file_handler = logging.FileHandler('99_logs/pipeline.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
```

#### 2.3.2 日志流向|Log Flow

```
Python Logger (ModuleName)
    │
    ├─ stdout handler (INFO)
    │   └─→ 系统stdout → 超算捕获 → YYYY-MM-DD-HH-MM-jobname.out
    │       内容：正常流程日志|Content: Normal operation logs
    │
    ├─ stderr handler (WARNING+)
    │   └─→ 系统stderr → 超算捕获 → YYYY-MM-DD-HH-MM-jobname.err
    │       内容：错误和警告|Content: Errors and warnings
    │
    └─ file handler (DEBUG+)
        └─→ 本地文件 → 99_logs/pipeline.log
            内容：完整调试日志|Content: Complete debug logs
```

#### 2.3.3 常见错误|Common Mistakes

❌ **错误1：所有日志都输出到stderr**
```python
# 导致 .out 文件为空，.err 文件包含所有信息
handler = logging.StreamHandler(sys.stderr)  # ❌
```

❌ **错误2：stdout 和 stderr 都输出INFO**
```python
# 导致 .out 和 .err 内容重复
stdout_handler.setLevel(logging.INFO)   # ✅
stderr_handler.setLevel(logging.INFO)    # ❌ 应该是 WARNING
```

❌ **错误3：忘记 propagate=False**
```python
# 导致日志重复输出
logger.propagate = True  # ❌ 会传播到root logger
logger.propagate = False  # ✅ 避免重复
```

### 2.4 中英文对照规范

**所有输出信息必须中英文对照，使用 `|` 分隔**

✅ **正确示例：**
```python
logger.info("开始分析流程|Starting analysis pipeline")
logger.info(f"处理样本|Processing sample: {sample_name}")
logger.warning("检测到低质量区域|Detected low quality region")
logger.error("文件读取失败|Failed to read file")
```

❌ **错误示例：**
```python
logger.info("开始分析流程")  # 缺少英文
logger.info("Starting analysis pipeline")  # 缺少中文
logger.info("开始分析流程 / Starting analysis")  # 错误分隔符
```

---

## 三、命令行参数规范

### 3.1 标准参数命名

| 短参数 | 长参数 | 类型 | 说明 |
|--------|--------|------|------|
| `-i` | `--input` | str | 输入文件或目录 |
| `-o` | `--output-dir` | str | 输出目录 |
| `-g` | `--genome` | str | 基因组FASTA文件 |
| `-t` | `--threads` | int | **线程数(默认12)** |
| `-k` | `--kmer-size` | int | K-mer大小 |
| `-h` | `--help` | flag | 帮助信息 |

### 3.2 CLI包装器规范(Click)

使用Click框架创建CLI包装器，采用懒加载模式：

```python
import click
import sys

def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...module.main import main as module_main
        return module_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)

def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path

@click.command(short_help="简短描述|Short description")
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='输入文件或目录|Input file or directory')
@click.option('-o', '--output-dir',
              default='./output',
              help='输出目录|Output directory')
def module_name(input, output_dir):
    """功能描述|Function description"""

    module_main = _lazy_import_main()

    # 构建参数列表
    args = ['module_name.py']
    args.extend(['-i', input])
    if output_dir != './output':
        args.extend(['-o', output_dir])

    original_argv = sys.argv
    sys.argv = args

    try:
        module_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
```

### 3.3 Help文档格式规范

**示例只展示最简单的用法即可，保持简洁清晰**

✅ **正确示例：**
```python
def module_name(input, output_dir):
    """
    功能描述|Function description

    示例|Examples: biopytools module -i input.txt -o output_dir/
    """
```

❌ **错误示例：**
```python
def module_name(input, output_dir):
    """
    功能描述|Function description

    示例|Examples:

    基本用法|Basic usage:
        biopytools module -i input.txt -o output_dir/

    高级用法|Advanced usage:
        biopytools module -i input.txt -o output_dir/ -t 24 --custom-param value

    超长命令示例|Long command example:
        biopytools module -i input.txt -o output_dir/ \\
            --param1 value1 --param2 value2 --param3 value3
    """
```

**格式要求：**
- `示例|Examples:`后直接跟命令，不换行
- 示例命令与`Examples:`之间用空格分隔
- 命令保持在同一行，不要使用反斜杠换行

**原因：**
- Help文档主要用于快速查看基本用法，不是完整教程
- 过多的示例会使help文档冗长，不易阅读
- 复杂用法应在README或独立文档中详细说明
- 简洁的示例在终端中显示更清晰，不会出现混乱的换行
- 避免缩进导致的额外空格显示

**最佳实践：**
- **只展示最常用的1个基本示例**
- 示例命令保持简短，能够在一行内完整显示
- 所有参数的说明已在参数列表中，示例中不需要重复展示
- 如有复杂用法，在文档中引导用户查看完整文档

---

## 四、命名规范

### 4.1 文件和目录命名

- 使用小写字母和下划线：`module_name.py`
- 模块目录名使用下划线：`bam_stats/`, `merqury_qv/`

### 4.2 变量和函数命名

- 使用小写字母和下划线：`fastq_files`, `calculate_qv()`
- 类名使用驼峰命名：`MerquryQVCalculator`, `ModuleLogger`

### 4.3 常量命名

- 使用大写字母和下划线：`DEFAULT_THREADS`, `MAX_DEPTH`

---

## 五、代码风格规范

### 5.1 字符串和注释

**所有字符串、注释、文档必须中英文对照**

```python
def calculate_qv_value():
    """
    计算QV值|Calculate QV value

    Returns:
        dict: QV统计结果|QV statistics results
    """
    pass
```

### 5.2 禁止使用emoji

❌ **禁止：**
```python
logger.info("开始处理🚀")
logger.info("完成✅")
```

✅ **正确：**
```python
logger.info("开始处理|Starting processing")
logger.info("完成|Completed")
```

### 5.3 数字格式化

**大数字使用M(百万)单位显示**

```python
def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)

# 输出示例
logger.info(f"总reads数|Total reads: {format_number(10000000)}")  # 10.00M
```

---

## 六、错误处理规范

### 6.1 异常捕获

```python
try:
    result = subprocess.run(command, check=True, capture_output=True)
except subprocess.CalledProcessError as e:
    logger.error(f"命令执行失败|Command execution failed: {e.stderr}")
    return False
except FileNotFoundError as e:
    logger.error(f"文件不存在|File not found: {e}")
    return False
```

### 6.2 输入验证

```python
def validate(self):
    """验证配置参数|Validate configuration parameters"""
    errors = []

    if not os.path.exists(self.input_file):
        errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

    if self.threads <= 0:
        errors.append(f"线程数必须为正数|Thread count must be positive")

    if errors:
        raise ValueError("\n".join(errors))

    return True
```

---

## 七、代码改造检查清单

### 7.1 日志格式检查
- [ ] 使用点号(.)分隔秒和毫秒
- [ ] 不使用方括号[]
- [ ] 使用 ` - ` 分隔各部分
- [ ] 所有日志信息中英文对照
- [ ] 不使用emoji
- [ ] 时间戳后面不能换行

### 7.2 代码结构检查
- [ ] 采用模块化结构(config/utils/calculator/main)
- [ ] 使用dataclass定义配置
- [ ] 使用标准的日志管理器类
- [ ] 分离核心逻辑和CLI接口

### 7.3 命名规范检查
- [ ] 所有字符串中英文对照
- [ ] 函数和变量使用snake_case
- [ ] 类名使用CamelCase
- [ ] 不使用emoji

### 7.4 功能检查
- [ ] 输入验证
- [ ] 异常处理
- [ ] 路径存在性检查
- [ ] 大数字使用M单位显示

### 7.5 命令执行日志检查 ⚠️ **重要|IMPORTANT**

- [ ] **所有外部命令执行前记录完整命令（INFO级别）** ⚠️ **必须检查**
- [ ] **禁止只记录描述，不记录命令** ❌ **严格禁止**
- [ ] 命令包含conda run包装（如果使用）
- [ ] 管道命令用 `|` 连接显示
- [ ] 格式清晰，便于复制到论文

**检查方法|Verification Method:**
```bash
# 搜索命令执行点
grep -n "subprocess.run\|Popen" biopytools/module_name/*.py

# 检查是否记录了命令
# 应该在命令执行前看到类似这样的代码：
# logger.info(f"命令|Command: {' '.join(cmd)}")
```

**正确代码示例|CORRECT Code Example:**
```python
def run_command(self, cmd, description):
    # ✅ 先记录描述
    if description:
        self.logger.info(f"执行|Executing: {description}")

    # ✅ 记录完整命令（INFO级别，不是DEBUG）
    self.logger.info(f"命令|Command: {' '.join(cmd)}")

    # 执行命令
    result = subprocess.run(cmd, ...)
```

**错误代码示例|WRONG Code Example:**
```python
def run_command(self, cmd, description):
    # ❌ 只记录描述，不记录命令
    if description:
        self.logger.info(f"执行|Executing: {description}")
        self.logger.debug(f"命令|Command: {' '.join(cmd)}")  # ❌ DEBUG级别看不到！

    # 执行命令
    result = subprocess.run(cmd, ...)
```

### 7.6 路径管理检查 ⚠️ **关键|CRITICAL**
- [ ] **无硬编码绝对路径|No hardcoded absolute paths** (必须检查|Must verify)
- [ ] 所有默认路径使用~展开|All default paths use ~ expansion
- [ ] 工具路径使用get_tool_path()|Tool paths use get_tool_path()
- [ ] 用户输入路径使用expand_path()|User input paths use expand_path()
- [ ] **⚠️ __post_init__中展开所有~路径|Expand all ~ paths in __post_init__** (关键|CRITICAL)
  ```python
  def __post_init__(self):
      self.tool_path = expand_path(self.tool_path)  # 必须|Must
      self.genome_path = expand_path(self.genome_path)  # 必须|Must
  ```
- [ ] 运行迁移脚本验证|Verify with migration script:
  ```bash
  grep -r "/share/org/YZWL/yzwl_lixg/" biopytools/module_name/
  # 应该返回0个结果|Should return 0 results
  ```
- [ ] 验证~路径已正确展开|Verify ~ paths are properly expanded:
  ```python
  # 测试代码|Test code
  config = ModuleConfig(tool_path="~/miniforge3/bin/tool")
  print(config.tool_path)  # 应该显示完整路径|Should show full path, not ~/...
  assert config.tool_path.startswith("/")  # 应该以/开头，不是~|Should start with /, not ~
  ```

**❌ 任何硬编码绝对路径都会导致代码审查失败|Any hardcoded absolute path will cause code review to FAIL**

---

## 八、标准代码模板

### 8.1 完整模块模板

```python
# config.py
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import os

@dataclass
class ModuleConfig:
    """配置类|Config Class"""

    input_dir: str
    output_dir: str = "./output"
    threads: int = 24

    def __post_init__(self):
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        errors = []
        if not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found")
        if errors:
            raise ValueError("\n".join(errors))
```

```python
# utils.py
import logging
import sys

class ModuleLogger:
    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file  # 必须在setup_logging之前赋值|Must assign before setup_logging
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(level=level, format=log_format,
                          datefmt=date_format, handlers=handlers)
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        return self.logger
```

---

## 九、常见问题

### Q1: 为什么日志格式要用点号而不是逗号？
**A:** 点号(.)更符合小数点的视觉习惯。Python logging的默认格式使用逗号，但我们需要通过 `%(msecs)03d` 显式指定点号格式。

### Q2: 中英文对照的顺序是什么？
**A:** 中文在前，英文在后，使用 `|` 分隔。例如：`"开始分析|Starting analysis"`

### Q3: 可以使用emoji吗？
**A:** 不可以。emoji可能在某些终端显示异常，也不利于日志解析。

### Q4: 如何格式化大数字？
**A:** 大于1百万的数字使用M单位，保留2位小数。例如：`10.00M` 而不是 `10000000`

---

## 十、开发约束与工作流规范

### 10.1 版本控制规范|Version Control

> ⚠️ **代码提交只在本地 Mac 进行；超算上严禁任何 Git 写操作**（详见文档顶部"禁止在超算上执行 Git 提交"）。
> Code is committed ONLY on the local Mac; Git write operations are forbidden on the supercomputer (see top-of-document warning).

超算上开发完成后，用 `copybiopytools` 同步到本地，由本地 Claude 统一 commit + push。推荐 commit message 格式如下：

```
<类型>(<模块>): <中文描述>

类型|Type:
  feat     新功能|New feature
  fix      Bug修复|Bug fix
  docs     文档变更|Documentation change
  refactor 重构（无功能变更）|Refactor (no functional change)
  test     测试相关|Test related
  chore    构建/工具变更|Build/tool change
```

**示例|Examples:**
```bash
git commit -m "feat(busco): 添加conda环境自动检测"
git commit -m "fix(path): 修复~路径在__post_init__中未展开的问题"
git commit -m "docs(guides): 更新路径管理规范至v2.4"
```

### 10.2 断点续传规范|Checkpoint Resume

**所有多步骤流程必须支持断点续传**，即已完成的步骤在重新运行时自动跳过，避免重复计算。

**标准实现模式|Standard Implementation:**

```python
def _is_step_completed(self, output_file: str) -> bool:
    """检查步骤是否已完成（通过输出文件存在性判断）|Check if step is done"""
    return Path(output_file).exists()

def run_step(self, step_name: str, output_file: str, run_func, *args, **kwargs):
    """
    带断点续传的步骤执行器|Step runner with checkpoint resume

    Args:
        step_name: 步骤名称|Step name
        output_file: 该步骤的关键输出文件|Key output file for this step
        run_func: 实际执行函数|Actual execution function
    """
    if self._is_step_completed(output_file):
        self.logger.info(f"跳过已完成步骤|Skipping completed step: {step_name}")
        return True

    self.logger.info(f"开始步骤|Starting step: {step_name}")
    return run_func(*args, **kwargs)
```

**使用示例|Usage:**
```python
# 运行jellyfish（若输出已存在则跳过）
self.run_step(
    "jellyfish_count",
    f"{output_dir}/01_jellyfish/{sample_id}.jellyfish.jf",
    self.run_jellyfish, sample_id
)
```

### 10.3 计算节点限制|Compute Node Restrictions

> ⚠️ **登录节点禁止运行大量计算任务！**

- 登录节点（Login Node）仅用于：代码编辑、文件管理、任务提交
- 所有需要大量CPU/内存的任务必须通过作业调度系统（如 `sub`）提交到计算节点
- 开发时的自动化测试只能测试：参数解析、路径验证、小型单元测试（mock subprocess）
- **大型计算测试必须手动在计算节点上运行**

---

## 十一（附）、测试规范|Testing Standards

### 测试目录结构|Test Directory Structure

```
biopytools/
├── module_name/
│   ├── ...
└── tests/
    ├── __init__.py
    ├── test_module_name/
    │   ├── __init__.py
    │   ├── test_config.py       # 配置类测试|Config tests
    │   ├── test_utils.py        # 工具函数测试|Utility tests
    │   └── test_calculator.py   # 核心逻辑测试|Core logic tests
    └── conftest.py              # 共享 fixtures
```

### 测试文件命名|Test File Naming

- 测试文件：`test_<被测模块名>.py`
- 测试函数：`test_<功能描述>()`
- 测试类：`Test<类名>()`

### 登录节点可运行的测试|Tests Safe for Login Node

使用 `unittest.mock` 替代真实子进程调用，避免在登录节点触发计算：

```python
# tests/test_module_name/test_utils.py
import unittest
from unittest.mock import patch, MagicMock
from biopytools.module_name.utils import build_conda_command, get_conda_env

class TestCondaDetection(unittest.TestCase):
    """conda路径检测测试|Conda path detection tests"""

    def test_get_conda_env_from_path(self):
        """测试从路径中检测conda环境|Test conda env detection from path"""
        with patch('shutil.which', return_value='/miniforge3/envs/BUSCO_v.6.0.0/bin/busco'):
            env = get_conda_env('busco')
            self.assertEqual(env, 'BUSCO_v.6.0.0')

    def test_get_conda_env_not_found(self):
        """测试非conda命令返回None|Test non-conda command returns None"""
        with patch('shutil.which', return_value='/usr/bin/bash'):
            with patch.dict('os.environ', {}, clear=True):
                env = get_conda_env('bash')
                self.assertIsNone(env)

    def test_build_conda_command(self):
        """测试conda命令构建|Test conda command building"""
        with patch('biopytools.module_name.utils.get_conda_env', return_value='BUSCO_v.6.0.0'):
            cmd = build_conda_command('busco', ['--version'])
            self.assertEqual(cmd, ['conda', 'run', '-n', 'BUSCO_v.6.0.0', 'busco', '--version'])

class TestPathExpansion(unittest.TestCase):
    """路径展开测试|Path expansion tests"""

    def test_expand_tilde(self):
        """测试~展开|Test ~ expansion"""
        from biopytools.common.paths import expand_path
        result = expand_path("~/test/path")
        self.assertTrue(result.startswith("/"))
        self.assertNotIn("~", result)

    def test_config_validates_input(self):
        """测试配置验证（mock文件存在）|Test config validation with mocked file"""
        from biopytools.module_name.config import ModuleConfig
        with patch('pathlib.Path.exists', return_value=True):
            config = ModuleConfig(input_file="fake.fastq", output_dir="/tmp/test")
            self.assertIsNone(config.validate())  # 不抛出异常|Should not raise
```

### 大型计算测试说明|Large Computation Test Notice

> ⚠️ **涉及真实工具调用（BUSCO、jellyfish等）的集成测试不能在登录节点运行，需手动在计算节点执行。**

---



### 11.1 路径配置优先级

所有工具路径必须遵循以下优先级（从高到低）|All tool paths must follow this priority (high to low):

1. **环境变量|Environment Variables** (最高|Highest)
2. **用户配置文件|User Config File** (`~/.config/biopytools/config.yml`)
3. **代码默认值|Code Defaults** (必须支持~展开|Must support ~ expansion)

### 11.2 禁止硬编码绝对路径

❌ **禁止|Forbidden**:
```python
# 硬编码绝对路径|Hardcoded absolute path
tool_path: str = "/share/org/YZWL/yzwl_lixg/miniforge3/envs/fanc_v.0.9.23b/bin/fanc"
```

✅ **正确|Correct**:
```python
# 使用路径管理工具|Use path management utility
from common.paths import get_tool_path

tool_path: str = get_tool_path(
    'fanc',                                    # 工具名|Tool name
    '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc',  # 默认值|Default (支持~|supports ~)
    'FANC_PATH'                                # 环境变量|Env var (可选|optional)
)
```

### 11.3 使用路径管理模块

导入路径管理工具|Import path management utilities:
```python
from common.paths import get_tool_path, expand_path
```

#### 11.3.0 `common/paths.py` 完整实现|Full Implementation

> 以下为 `common/paths.py` 的参考实现，所有模块直接从此导入|Reference implementation for all modules to import from:

```python
# common/paths.py
"""路径管理工具|Path management utilities"""

import os
import shutil
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def expand_path(path: str) -> str:
    """
    展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 原始路径（可包含~或$VAR）|Raw path (may contain ~ or $VAR)

    Returns:
        展开后的绝对路径|Expanded absolute path

    Examples:
        >>> expand_path("~/miniforge3/bin/tool")
        "/home/user/miniforge3/bin/tool"
        >>> expand_path("$SOFTWARE/bin/tool")
        "/opt/software/bin/tool"
    """
    return os.path.expandvars(os.path.expanduser(path))


def get_tool_path(
    tool_name: str,
    default_path: str,
    env_var: Optional[str] = None
) -> str:
    """
    按优先级获取工具路径|Get tool path by priority

    优先级（高→低）|Priority (high→low):
        1. 环境变量 (env_var)|Environment variable
        2. 用户配置文件 (~/.config/biopytools/config.yml)
        3. 代码默认值 (default_path)|Code default

    Args:
        tool_name: 工具名称（用于配置文件查找）|Tool name (for config file lookup)
        default_path: 默认路径（支持~展开）|Default path (supports ~ expansion)
        env_var: 环境变量名（可选）|Environment variable name (optional)

    Returns:
        展开后的工具路径|Expanded tool path
    """
    # 1. 环境变量优先|Env var takes priority
    if env_var and os.environ.get(env_var):
        path = expand_path(os.environ[env_var])
        logger.debug(f"使用环境变量路径|Using env var path [{env_var}]: {path}")
        return path

    # 2. 用户配置文件|User config file
    config_path = os.path.expanduser("~/.config/biopytools/config.yml")
    if os.path.exists(config_path):
        try:
            import yaml
            with open(config_path) as f:
                config = yaml.safe_load(f) or {}
            tool_path = config.get("tools", {}).get(tool_name)
            if tool_path:
                path = expand_path(tool_path)
                logger.debug(f"使用配置文件路径|Using config file path [{tool_name}]: {path}")
                return path
        except Exception as e:
            logger.warning(f"读取配置文件失败|Failed to read config file: {e}")

    # 3. 代码默认值|Code default
    path = expand_path(default_path)
    logger.debug(f"使用默认路径|Using default path [{tool_name}]: {path}")
    return path
```

**config.py 中使用|Use in config.py**:
```python
from dataclasses import dataclass, field
from common.paths import get_tool_path

@dataclass
class ModuleConfig:
    tool_path: str = field(
        default_factory=lambda: get_tool_path(
            'tool_name',
            '~/default/path/to/tool',
            'TOOL_PATH'  # 可选环境变量|Optional env var
        )
    )
```

**main.py 中使用|Use in main.py**:
```python
from common.paths import expand_path

# 展开用户提供的路径|Expand user-provided path
config.tool_path = expand_path(args.tool_path)
```

### 11.3.1 ⚠️ 关键：必须展开~路径|CRITICAL: Must Expand ~ Paths

**问题|Problem:**
Python不会自动展开字符串中的`~`符号|Python does NOT auto-expand `~` in strings:

```python
# ❌ 错误示例|Wrong Example
meme_path: str = "~/miniforge3/envs/meme_v.5.5.9/bin/meme"

# 在验证时会失败|Will fail in validation:
import shutil
shutil.which(meme_path)  # 返回None|Returns None ❌
os.path.exists(meme_path)  # 返回False|Returns False ❌
```

**解决方案|Solution:**
在`config.py`的`__post_init__`方法中必须展开所有包含`~`的路径|Must expand all paths containing `~` in `__post_init__`:

```python
from dataclasses import dataclass
from ..common.paths import expand_path

@dataclass
class ModuleConfig:
    meme_path: str = '~/miniforge3/envs/meme_v.5.5.9/bin/meme'
    tool_path: str = '~/software/bin/tool'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # ⚠️ 关键：展开所有包含~的路径|CRITICAL: Expand all paths with ~
        self.meme_path = expand_path(self.meme_path)
        self.tool_path = expand_path(self.tool_path)

        # 现在可以正常检测工具|Now tool detection works:
        # shutil.which(self.meme_path)  # 返回完整路径|Returns full path ✓
```

**检查清单|Checklist:**
- [ ] 所有包含`~`的默认路径都在`__post_init__`中调用`expand_path()`
- [ ] 所有从命令行/用户输入的路径都使用`expand_path()`展开
- [ ] `shutil.which()`, `os.path.exists()`等检测前路径已展开

**常见错误|Common Mistakes:**
1. ❌ 只在定义时使用`expand_path()`，不在`__post_init__`中展开
2. ❌ 忘记展开用户通过命令行参数传入的路径
3. ❌ 使用`os.path.expanduser()`而不是统一的`expand_path()`

### 11.4 路径格式要求

- ✅ 支持~展开|Support ~ expansion: `~/miniforge3/bin/tool`
- ✅ 支持环境变量|Support env vars: `$SOFTWARE/tool/bin`
- ✅ 使用相对路径|Use relative paths: `./tools/tool`
- ❌ 禁止用户名路径|Forbidden user path: `/share/org/YZWL/...`

### 11.5 环境变量命名规范

格式|Format: `{TOOL_NAME}_PATH` (全大写|UPPERCASE)

示例|Examples:
- `FANC_PATH`
- `BWA_PATH`
- `SAMTOOLS_PATH`
- `AUGUSTUS_PATH`

### 11.6 配置文件格式

用户配置文件位置|User config file location: `~/.config/biopytools/config.yml`

格式示例|Format example:
```yaml
tools:
  fanc: ~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
  bwa: ~/.local/bin/bwa
  samtools: ~/.local/bin/samtools

databases:
  nr: ~/database/ncbi/nr
```

### 11.7 迁移现有代码

使用迁移脚本|Use migration script:
```bash
cd ~/software/biopytools
python scripts/migrate_paths.py --module module_name --dry-run  # 先试运行|Dry run first
python scripts/migrate_paths.py --module module_name --apply   # 确认后应用|Apply after review
```


## 十二、输出目录和文件命名规范|Output Directory and File Naming Conventions

> **参考标准|Reference Standards:** nf-core, GATK Best Practices, Bioconda Community Guidelines  
> **目的|Purpose:** 统一模块输出结构，提高可读性、可维护性和下游分析便利性

---

### 12.1 通用原则|General Principles

#### 12.1.1 机器友好|Machine-Friendly

- ✅ **禁止空格|No Spaces:** 永远不要在文件名或文件夹名中使用空格
  - ❌ `my output/` → ✅ `my_output/` 或 `my-output/`
- ✅ **避免特殊字符|Avoid Special Characters:** 避免 `&`, `*`, `(`, `)`, `,` `@` 等 Shell 敏感字符
- ✅ **大小写敏感|Case-Sensitive:** 统一使用**小写**（Sample ID 除外）
- ✅ **使用分隔符|Use Separators:** 下划线 `_` 用于文件名内部；流程步骤目录使用 `序号_步骤名` 格式（如 `01_qc/`），普通目录名使用下划线或中划线均可，但**同一项目内保持一致**

#### 12.1.2 人类可读|Human-Readable

- ✅ **见名知意|Self-Descriptive:** 拒绝 `output1`, `step2`, `result.txt`
- ✅ **有序性|Ordered:** 流程步骤使用数字前缀（`01_`, `02_`）

#### 12.1.3 可追溯性|Traceability

- ✅ **样本信息|Sample Info:** 文件名包含样本 ID
- ✅ **处理步骤|Processing Step:** 标注使用的工具或流程
- ✅ **关键参数|Key Parameters:** 必要时标注参数（如 `k21`）

---

### 12.2 输出文件夹命名规范|Output Directory Naming

#### 12.2.1 结构原则|Structure Principles

**多样本分析采用 "By Sample" (按样本)** 层级结构，**单样本或共享分析采用 "By Step" (按步骤)**：

```python
# ✅ 多样本分析|Multi-sample Analysis
output/
├── 00_pipeline_info/          # 全局元数据（可选）
├── sample1/
│   ├── 00_pipeline_info/      # 样本级元数据
│   ├── 01_qc/
│   ├── 02_trimming/
│   └── 99_logs/
├── sample2/
│   ├── 00_pipeline_info/
│   ├── 01_qc/
│   └── 99_logs/
└── sample3/

# ✅ 单样本或共享分析|Single Sample or Shared Analysis
output/
├── 00_pipeline_info/
├── 01_qc/
├── 02_trimming/
├── 03_alignment/
└── 99_logs/

# ❌ 不推荐|Not Recommended
output/
├── qc/
├── sample1/
├── sample2/
└── trimming/
```

#### 12.2.2 命名模式|Naming Pattern

> **格式|Format:** `[序号]_[步骤名称]`

| 组件|Component | 说明|Description | 示例|Example |
|-----|-----------|-----|-------------|------|--------|
| 序号|Number | 两位数字|Two-digit number | `01`, `02`, `03` | |
| 步骤名|Step Name | 小写英文|Lowercase English | `qc`, `alignment` | `01_qc/` |

#### 12.2.3 特殊子目录|Special Subdirectories

| 目录名|Directory | 用途|Purpose | 优先级|Priority |
|---------|-----------|-----|---------|---------|----------|
| `00_pipeline_info/` | 流程元数据|软件版本、参数记录 | P0 | 必须 |
| `99_logs/` | 日志文件|运行日志、错误日志 | P0 | 必须 |
| `work/` | 中间文件|临时文件，可删除 | P1 | 建议 |

> **注|Note:** `00_` 和 `99_` 前缀表示流程元数据和日志分别在流程开始和结束位置

#### 12.2.4 无编号目录|Unnumbered Directories

某些支持性目录**不应**使用数字前缀：

| 场景|Scenario | 示例|Example | 原因|Reason |
|---------|-----------|-----|-------------|------|--------|
| 条件性步骤|Conditional step | `fastk/` (仅Smudgeplot需要) | 不是所有运行都会创建此目录 |
| 预生成输入|Pre-generated input | `reference/`, `index/` | 外部提供的输入文件，非流程输出 |
| 工具缓存|Tool cache | `cache/`, `tmp/` | 中间缓存，非结果文件 |

> **示例|Example:** `fastk/` 目录
> - FastK 是 Smudgeplot 的预处理步骤（非独立分析）
> - 支持通过 `--fastk-table` 使用预生成表（路径不应含编号）
> - 使用 `--skip-smudgeplot` 时不存在（条件性创建）

---

### 12.3 输出文件命名规范|Output File Naming

#### 12.3.1 基本命名模式|Basic Naming Pattern

> **格式|Format:** `{Sample_ID}.{Tool/Step}.{State}.{Extension}`

| 组件|Component | 说明|Description | 示例|Example |
|-----|-----------|-----|-------------|------|--------|
| Sample_ID | 样本编号 | **必须在最前面**，保持一致性 | `R0590-6` | |
| Tool/Step | 工具/步骤 | 使用的软件或处理步骤 | `jellyfish`, `genomescope` | |
| State | 状态（可选）| 处理状态或文件类型 | `sorted`, `filtered`, `model` | |
| Extension | 扩展名 | 标准文件扩展名 | `.jf`, `.png`, `.txt` | |

#### 12.3.2 后缀叠加原则|Suffix Stacking

文件经过多步处理时，后缀可叠加（参考 GATK 风格）：

```python
# 示例：BAM文件处理流程
sample.bam                    # 原始比对
sample.sorted.bam             # 排序后
sample.sorted.markdup.bam     # 去重后
sample.sorted.markdup.recal.bam  # 重校正后
```

#### 12.3.3 扩展名规范|Extension Standards

| 文件类型|File Type | 推荐命名|Recommended Naming | 示例|Example |
|---------|-----------|---------|-------------------|------|--------|
| **FASTQ** | 原始数据 | `{ID}_R{1/2}.fastq.gz` | `S1_R1.fastq.gz` | |
| **质控** | FastQC报告 | `{ID}.fastqc.html` | `S1.fastqc.html` | |
| **清洗后** | 过滤数据 | `{ID}.trimmed.fq.gz` | `S1.trimmed.fq.gz` | |
| **Jellyfish** | K-mer数据库 | `{ID}.jellyfish.jf` | `S1.jellyfish.jf` | |
| **Jellyfish** | 直方图 | `{ID}.jellyfish.histo` | `S1.jellyfish.histo` | |
| **GenomeScope** | 模型文件 | `{ID}.genomescope.model.txt` | `S1.genomescope.model.txt` | |
| **GenomeScope** | 图表 | `{ID}.genomescope.linear.png` | `S1.genomescope.linear.png` | |
| **Smudgeplot** | 结果图 | `{ID}.smudgeplot.linear.png` | `S1.smudgeplot.linear.png` | |
| **日志** | 运行日志 | `{ID}.{tool}.log` | `S1.genomescope.log` | |

#### 12.3.4 常见错误|Common Mistakes

```python
# ❌ 错误示例|Bad Examples
genome_analysis.jf              # 无样本名，太通用
smudgeplot_smudgeplot.png      # 重复前缀
plot.png                        # 无上下文信息
result.txt                      # 毫无意义

# ✅ 正确示例|Good Examples
R0590-6.jellyfish.jf           # 样本名.工具名.扩展名
R0590-6.smudgeplot.linear.png  # 样本名.工具名.类型.扩展名
R0590-6.genomescope.model.txt  # 样本名.工具名.结果类型.扩展名
```

---

### 12.4 临时文件管理|Temporary File Management

#### 12.4.1 临时文件放置|Temporary File Placement

```python
# ✅ 推荐|Recommended
import tempfile

# 使用输出目录下的tmp子目录(超算系统/tmp易爆满,统一用output_dir/tmp)|
# Use tmp subdir under output dir (HPC system /tmp fills up easily)
tmp_root = os.path.join(output_dir, "tmp")
os.makedirs(tmp_root, exist_ok=True)
with tempfile.TemporaryDirectory(prefix=f"module_{sample_id}_", dir=tmp_root) as tmpdir:
    # 临时操作
    process_files(tmpdir)
    # 自动清理(退出with即清理)|Auto-cleaned on with-exit

# ❌ 不推荐|Not Recommended  
output_dir = "results/"
os.makedirs(output_dir + "temp/", exist_ok=True)  # 临时文件混在结果中
```

#### 12.4.2 必须清理的场景|Must Cleanup Scenarios

| 场景|Scenario | 临时文件|Temporary Files | 清理时机|Cleanup Timing |
|---------|-----------|---------|----------------|---------|---------------|
| FastK输入|FastK input | 解压的 `.fq` 文件|Decompressed `.fq` | FastK完成后|After FastK |
| 压缩检查|Compression check | 解压用于检查的文件|Files for checking | 检查完成后|After check |
| 中间步骤|Intermediate | 流程中间文件|Pipeline intermediates | 最终结果生成后|After final output |

#### 12.4.3 实现示例|Implementation Example

```python
import tempfile
import shutil
from contextlib import contextmanager

@contextmanager
def temp_directory(prefix, cleanup=True):
    """临时目录上下文管理器|Temporary directory context manager"""
    temp_dir = tempfile.mkdtemp(prefix=prefix)
    try:
        yield temp_dir
    finally:
        if cleanup and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

# 使用|Usage
with temp_directory(f"genomescope_{sample_id}_") as tmpdir:
    # 解压到临时目录
    decompressed_files = decompress(fastq_files, tmpdir)
    # 运行FastK
    run_fastk(decompressed_files)
    # 上下文退出时自动清理
```

---

### 12.5 版本信息记录|Version Information

#### 12.5.1 禁止文件名包含版本号|No Version in Filename

```python
# ❌ 不推荐|Not Recommended
sample_bwa_v0.7.17.bam         # 破坏下游脚本兼容性
sample_hifiasm_0.19.8.fa       # 难以维护

# ✅ 正确|Correct
sample.bwa.bam                 # 版本信息单独记录
sample.hifiasm.fa
```

#### 12.5.2 software_versions.yml 格式

```yaml
# 00_pipeline_info/software_versions.yml
pipeline:
  name: "biopytools module_name"
  version: "1.0.0"

tools:
  tool_name:
    version: "x.x.x"
    path: "~/miniforge3/envs/tool_env/bin/tool"  # 使用~展开|Use ~ expansion
    command: "tool --version"  # 版本检测命令

parameters:
  param1: value1
  param2: value2

execution:
  start_time: "2026-03-05 10:00:00"
  end_time: "2026-03-05 12:00:00"
  runtime_seconds: 7200
```

#### 12.5.3 生成版本信息的代码示例|Code Example

```python
import subprocess
import yaml
from datetime import datetime
from pathlib import Path

def generate_software_versions_yml(output_dir: str, tools: dict, params: dict, start_time: datetime = None):
    """生成software_versions.yml文件|Generate software_versions.yml file"""

    if start_time is None:
        start_time = datetime.now()

    versions = {}
    for tool_name, tool_path in tools.items():
        try:
            result = subprocess.run(
                [tool_path, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            versions[tool_name] = {
                'version': result.stdout.strip(),
                'path': tool_path
            }
        except Exception:
            versions[tool_name] = {'version': 'unknown', 'path': tool_path}

    end_time = datetime.now()
    runtime_seconds = int((end_time - start_time).total_seconds())

    info = {
        'pipeline': {
            'name': 'biopytools module_name',
            'version': '1.0.0'
        },
        'tools': versions,
        'parameters': params,
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': runtime_seconds
        }
    }
    
    output_file = Path(output_dir) / '00_pipeline_info' / 'software_versions.yml'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        yaml.dump(info, f, default_flow_style=False)
```

---

### 12.6 完整示例|Complete Example

#### 12.6.1 推荐的输出结构|Recommended Output Structure

```
02.output/R0590-6/
├── 00_pipeline_info/                      # 流程元数据
│   ├── software_versions.yml              # 软件版本
│   └── pipeline_params.yaml               # 运行参数
│
├── 01_jellyfish/                           # 步骤1：K-mer计数
│   ├── R0590-6.jellyfish.jf               
│   └── R0590-6.jellyfish.histo           
│
├── 02_genomescope/                         # 步骤2：基因组特征分析
│   ├── R0590-6.genomescope.model.txt      
│   ├── R0590-6.genomescope.summary.txt    
│   ├── R0590-6.genomescope.linear.png     
│   └── R0590-6.genomescope.log.png       
│
├── 03_smudgeplot/                          # 步骤3：倍性分析
│   ├── R0590-6.smudgeplot.kmerpairs.smu  
│   ├── R0590-6.smudgeplot.report.tsv     
│   ├── R0590-6.smudgeplot.linear.png     
│   └── R0590-6.smudgeplot.log10.png      
│
└── 99_logs/                                # 日志文件
    └── genomescope_pipeline.log           
```

#### 12.6.2 genomescope模块改进对照表|genomescope Module Improvement

| 组件|Component | 当前|Current | 改进|Improved | 说明|Notes |
|-----|-----------|-----------|-----------|------|--------|
| **目录结构**|Directory | `jellyfish/` | `01_jellyfish/` | 添加序号前缀|Add number prefix |
| | | `genomescope_output/` | `02_genomescope/` | 简化+序号|Simplify + number |
| | | `smudgeplot_output/` | `03_smudgeplot/` | 简化+序号|Simplify + number |
| **文件命名**|File Names | `genome_analysis.jf` | `R0590-6.jellyfish.jf` | 样本名.工具名|Sample.tool |
| | | `genome_analysis.histo` | `R0590-6.jellyfish.histo` | 样本名.工具名|Sample.tool |
| | | `model.txt` | `R0590-6.genomescope.model.txt` | 明确文件类型|Specify file type |
| | | `plot.png` | `R0590-6.genomescope.linear.png` | 明确比例尺|Specify scale |
| | | `smudgeplot_smudgeplot.png` | `R0590-6.smudgeplot.linear.png` | 去除重复|Remove duplicate |
| **临时文件**|Temp Files | `fastk/*.fq` (68GB) | `<output>/tmp/*.fq` (运行结束清理) | output_dir/tmp子目录|Use output/tmp |
| **版本信息**|Version Info | ❌ 无|Missing | `software_versions.yml` | 添加版本文件|Add version file |
| **日志管理**|Log Management | 散落各处|Scattered | `99_logs/pipeline.log` | 集中管理|Centralize |

---

### 12.7 检查清单|Checklist

在开发新模块或修改现有模块时，请确认：

- [ ] ✅ 所有文件夹名无空格，使用小写
- [ ] ✅ 流程步骤文件夹使用数字前缀（`01_`, `02_`）
- [ ] ✅ 元数据目录使用 `00_pipeline_info/`，日志目录使用 `99_logs/`
- [ ] ✅ 条件性步骤或支持性目录不使用编号（如 `fastk/`）
- [ ] ✅ 所有输出文件包含样本ID作为前缀
- [ ] ✅ 文件名遵循 `{Sample}.{Tool}.{State}.{Ext}` 格式
- [ ] ✅ 扩展名标准且包含压缩格式（`.fastq.gz`, `.bam`）
- [ ] ✅ 临时文件使用 `output_dir/tmp` 子目录,运行结束清理(避免超算系统 `/tmp` 爆满)
- [ ] ✅ 生成 `00_pipeline_info/software_versions.yml`
- [ ] ✅ 日志文件统一放置在 `99_logs/` 目录
- [ ] ✅ 禁止在文件名中包含软件版本号
- [ ] ✅ 避免重复前缀（如 `smudgeplot_smudgeplot.png`）

---

---

## 十三、Conda环境软件调用规范|Conda Environment Software Invocation

### 13.1 问题描述|Problem Description

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

### 13.2 解决方案|Solution

使用`conda run -n <env_name>`包装命令来调用conda环境中的软件。

#### 13.2.0 ⚠️ 重要：conda run必须使用--no-capture-output参数|IMPORTANT: conda run must use --no-capture-output

**⚠️ 严格禁止|STRICTLY FORBIDDEN: 不使用--no-capture-output会导致内存问题**

❌ **错误示例 - 导致CondaMemoryError|WRONG - Causes CondaMemoryError:**
```python
# 🔴 绝对禁止|ABSOLUTELY FORBIDDEN
full_cmd = ['conda', 'run', '-n', 'population_genetics', 'bwa', 'mem', '-t', '64', ref_fa, r1, r2]
# 问题：conda默认缓冲所有输出，对于大数据会导致内存溢出
# Error: CondaMemoryError: The conda process ran out of memory
```

✅ **正确做法|CORRECT Approach:**
```python
# ✅ 正确：添加--no-capture-output参数
full_cmd = ['conda', 'run', '-n', 'population_genetics', '--no-capture-output', 'bwa', 'mem', '-t', '64', ref_fa, r1, r2]
# --no-capture-output告诉conda不要缓冲输出，直接透传到stdout/stderr
```

**为什么必须使用--no-capture-output？|Why must use --no-capture-output?**

1. **避免内存溢出**：conda run默认缓冲所有输出到内存，大数据（如BWA比对32GB Hi-C数据）会导致OOM
2. **流式输出**：--no-capture-output让输出直接流式传递，不占用额外内存
3. **实时进度**：用户可以看到实时进度，而不是等命令结束才显示所有输出

**检查清单|Checklist:**
- [ ] `build_conda_command`函数是否添加了`--no-capture-output`？
- [ ] 所有conda run命令都包含此参数？
- [ ] 文档示例代码已更新？

#### 13.2.1 🚨 严格禁止：管道中使用conda run|STRICTLY FORBIDDEN: conda run in Pipelines

**⚠️ 严禁在Shell管道中使用多个conda run命令|STRICTLY FORBIDDEN: Multiple conda run commands in shell pipelines**

❌ **错误示例 - 会导致BrokenPipeError和内存问题|WRONG - Causes BrokenPipeError and memory issues:**
```python
# 🔴 绝对禁止|ABSOLUTELY FORBIDDEN
cmd1 = build_conda_command(bwa_bin, bwa_args)      # ['conda', 'run', '-n', 'env1', 'bwa', ...]
cmd2 = build_conda_command(samtools_bin, args)     # ['conda', 'run', '-n', 'env2', 'samtools', ...]
pipeline = f"{ ' '.join(cmd1) } | { ' '.join(cmd2) }"  # conda run | conda run
```

**问题|Problems:**
1. **BrokenPipeError**: conda run会缓冲输出，导致管道破裂
2. **内存占用高**: conda run先将所有输出读入内存，再传递给下一个命令
3. **性能低下**: 双重conda包装严重影响性能

✅ **正确做法|CORRECT Approach:**

**方案A：管道命令直接调用（设置LD_LIBRARY_PATH）**
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

**方案B：使用_build_pipeline_command辅助函数**
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

**方案C：管道中的工具使用shutil.which自动检测**
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

**检查清单|Checklist:**
- [ ] 管道中是否有`conda run | conda run`模式？
- [ ] 是否正确设置了`LD_LIBRARY_PATH`？
- [ ] 是否使用了辅助函数或自动检测？
- [ ] 是否避免了`text=True`处理二进制数据（如BAM文件）？

### 13.3 实现模式|Implementation Pattern

#### 13.3.1 环境检测函数|Environment Detection Function

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

#### 13.3.2 命令构建函数|Command Building Function

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

#### 13.3.3 依赖检查函数|Dependency Check Function

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

#### 13.3.4 命令执行器|Command Runner

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

### 13.4 使用示例|Usage Example

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

### 13.5 关键注意事项|Key Considerations

1. **自动检测|Auto Detection**: 代码应自动检测conda环境，无需用户手动指定
2. **向后兼容|Backward Compatible**: 对于非conda环境安装的软件，仍应正常工作
3. **性能影响|Performance Impact**: `conda run`会有轻微的启动开销，通常可忽略
4. **路径传递|Path Passing**: 即使传入完整路径，`conda run`也能正确处理
5. **环境变量|Environment Variables**: 使用`CONDA_EXE`环境变量定位conda基础目录

### 13.6 常见错误与禁令|Common Mistakes and Prohibitions

#### 13.6.1 ❌ 禁止提取命令名|Prohibit Command Name Extraction

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

#### 13.6.2 ❌ 禁止使用命令名调用|Prohibit Using Command Name Only

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

#### 13.6.3 代码审查检查点|Code Review Checklist

在代码审查时，检查以下关键点：

- [ ] **禁止使用** `os.path.basename()` 提取命令名后调用 `build_conda_command()`
- [ ] **必须传递完整路径** 给 `build_conda_command(command, args)`
- [ ] **检查管道命令**中的每个子命令是否传递完整路径
- [ ] **验证生成的命令**包含 `-n <env_name>` 参数

**审查命令示例|Review Command Example:**
```python
# ❌ 代码审查应该拒绝的代码
cmd_name = os.path.basename(tool_path)  # 应该标记为错误
build_conda_command(cmd_name, args)     # 应该标记为错误

# ✅ 代码审查应该通过的代码
build_conda_command(tool_path, args)    # 正确：传递完整路径
```

### 13.7 测试验证|Testing

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

### 13.8 故障排查|Troubleshooting

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

---

---

## 版本历史

| 版本 | 日期 | 主要变更|Major Changes |
|------|------|----------|
| 2.15 | 2026-07-23 | 临时目录统一改造：所有模块临时文件/目录从系统 `/tmp` 改为 `output_dir/tmp` 子目录并运行结束清理，消除超算 `/tmp` 爆满风险；同步更新 12.4.1/12.6/12.7 |
| 2.14 | 2026-06-24 | 新增"禁止在超算上执行 Git 提交"规范（顶部醒目警告 + 修订 10.1）：明确超算只写代码不 commit，提交统一在本地 Mac 由 Claude 完成；配套 `copybiopytools` 增加 `--exclude='.git/'` |
| 2.13 | 2026-04-09 | 文档质量改进：修正版本号不一致(2.11→2.12)；修复测试断言缺少--no-capture-output；裸except改为except Exception；修复示例代码中的硬编码绝对路径（software_versions.yml模板、管道方案A）；补全generate_software_versions_yml的end_time/runtime_seconds字段使其与模板一致；添加目录导航 |
| 2.12 | 2026-03-17 | 完善Conda调用规范：新增13.2.0强制要求--no-capture-output参数，避免conda缓冲输出导致CondaMemoryError；更新build_conda_command示例代码 |
| 2.11 | 2026-03-17 | 新增命令执行日志规范（2.2.1）：强制要求所有外部命令执行前记录完整命令到INFO级别；禁止只记录描述不记录命令；添加代码审查检查点（7.5）；提升工具透明度和可重复性 |
| 2.10 | 2026-03-15 | 完善Conda调用规范：新增13.6常见错误与禁令、13.7更新测试示例、13.8故障排查；明确禁止os.path.basename()提取命令名；添加代码审查检查点 |
| 2.9 | 2026-03-05 | 修复CommandRunner shell=True Bug；补充ModuleLogger模板缺失赋值；补全common/paths.py实现；新增第十章开发约束规范（Git提交、断点续传、节点限制）；新增测试规范章节；统一目录分隔符说明；合并重复版本历史 |
| 2.8 | 2026-03-05 | 添加Conda环境软件调用规范（自动检测并使用conda run）|Add conda environment software invocation spec (auto-detect and use conda run) |
| 2.7 | 2026-03-05 | 添加超算日志系统分离规范（INFO→.out, WARNING→.err）|Add job scheduler log separation spec (INFO→.out, WARNING→.err) |
| 2.6 | 2026-03-05 | 澄清多样本vs单样本目录结构原则|Clarify multi-sample vs single-sample directory structure principles |
| 2.5 | 2026-03-05 | 添加输出目录和文件命名规范（参考nf-core、GATK最佳实践）|Add output directory and file naming conventions (ref: nf-core, GATK best practices) |
| 2.4 | 2026-03-02 | 添加~路径展开的强制要求和验证方法|Add mandatory ~ path expansion requirement and verification methods |
| 2.3 | 2026-03-02 | 添加文档开头的路径管理警告声明和代码审查检查项|Add path management warning at document start and code review checklist |
| 2.2 | 2026-03-02 | 添加路径管理规范（避免硬编码绝对路径）|Add path management spec (avoid hardcoded absolute paths) |
| 2.1 | 2026-01-05 | 添加Help文档格式规范|Add Help doc format spec |
| 2.0 | 2026-01-04 | 更新日志格式规范；添加模块化结构规范；强化中英文对照要求；添加大数字格式化规范 |
| 1.0 | 2024-12-19 | 初始版本|Initial version |