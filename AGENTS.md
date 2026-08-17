# BioPyTools Python代码开发规范文档

## 版本: 2.21
## 日期: 2026-08-17
## 用途: 统一所有生信分析模块的代码结构、命名规范、日志格式

> 📐 **文档结构说明|Doc Structure:** 本文件只保留**核心规则与检查清单**(每次会话全量加载)。
> 完整代码模板、实现细节、故障排查等重内容已拆到 `docs/dev-standards/` 下,**按需读取**(见文末「📚 详细参考文档」)。
> This file keeps only core rules & checklists (loaded every session). Heavy reference material lives in `docs/dev-standards/` — read on demand (see bottom).

## 目录|Table of Contents

- [⚠️ 重要警告|CRITICAL WARNING](#-重要警告critical-warning)
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
- [十四、模块文档规范](#十四模块文档规范)
- [📚 详细参考文档](#-详细参考文档)
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

**后果|Consequences:** 🔴 代码审查不通过 / 🔴 无法合并到主分支 / 🔴 破坏代码可移植性

**参考|Reference:** 详见第十一章节;`common/paths.py` 完整实现见 [docs/dev-standards/11_path_management.md](docs/dev-standards/11_path_management.md)

---

### 🚫 禁止在超算上执行 Git 提交|FORBIDDEN: Git Commits on the Supercomputer

> **⚠️ AI 在超算（登录节点/计算节点）上开发时，严禁执行任何 Git 写操作**

❌ **超算上严禁|STRICTLY FORBIDDEN:** `git add / commit / push / tag`
✅ **超算上允许|ALLOWED:** `git status / diff`，编辑代码、运行测试、提交计算任务

**原因|Reason:** 超算出网受限 push 常失败；超算上 `.git` 不会同步到本地（`copybiopytools` 已排除 `.git/`），在超算 commit 会产生分叉历史；**代码提交的唯一入口是本地 Mac**。

> 📌 一句话：**超算只写代码，Mac 才提交。|The supercomputer writes code only; the Mac is the only place commits happen.**

---

## 一、模块化结构规范

每个功能模块采用如下结构（完整代码模板见 [docs/dev-standards/01_module_template.md](docs/dev-standards/01_module_template.md)）：

```
biopytools/
└── module_name/
    ├── __init__.py      # 模块导出声明 + __version__
    ├── config.py        # @dataclass 配置类 + validate()
    ├── utils.py         # 日志管理器 + 工具函数 + 命令执行
    ├── calculator.py    # 核心计算逻辑(命名可灵活)
    └── main.py          # argparse/click 命令行入口
```

**职责要点|Key Points:**
- `config.py`：用 `@dataclass`，所有 `~` 路径在 `__post_init__` 中 `expand_path()`（见 §十一）
- `utils.py`：`ModuleLogger` 类 + 通用工具；`self.log_file` 必须在 `setup_logging` 之前赋值
- `calculator.py`：组合 config + logger，不含 CLI 逻辑
- `main.py`：解析参数 → 构造 config → validate → 运行；CLI 与核心逻辑分离

---

## 二、日志格式规范

### 2.1 标准日志格式

```
YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 消息中文|Message English
```

- 点号 `.` 分隔秒和毫秒（**不是逗号**），毫秒固定3位
- ` - ` (空格-空格) 分隔各部分（**不是方括号 []**）
- `|` 分隔中英文；时间戳后不换行

```python
log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
date_format = '%Y-%m-%d %H:%M:%S'
```

### 2.2 日志级别

| 级别 | 用途 | 示例 |
|------|------|------|
| **DEBUG** | 详细调试信息 | 变量值、中间结果 |
| **INFO** | 正常流程信息 | 步骤开始/完成、统计数据 |
| **WARNING** | 警告信息 | 质控警告、参数建议 |
| **ERROR** | 错误信息 | 文件读取失败、格式错误 |
| **CRITICAL** | 严重错误 | 程序无法继续的错误 |

#### 2.2.1 ⚠️ 命令执行必须记录完整命令|Command Execution Logging

**所有外部命令执行时，必须记录完整命令到 INFO 级别**（调试需求、论文 Methods 可重复性、避免黑箱）。

```python
# ✅ 先记录描述，再记录完整命令（INFO 级别，不是 DEBUG）
logger.info(f"执行|Executing: BWA索引|BWA indexing")
logger.info(f"命令|Command: {' '.join(cmd)}")
subprocess.run(cmd)

# ❌ 严禁：只记录描述，不记录命令；或把命令放 DEBUG 级别
```

管道命令也要完整记录：`logger.info(f"命令|Command: {' '.join(bwa_cmd)} | {' '.join(samtools_cmd)} > {output_file}")`

### 2.3 超算日志分离|HPC Log Separation

> **原则|Principle:** INFO → stdout → `.out`；WARNING/ERROR → stderr → `.err`
> 超算运行时必须配置 stdout/stderr/file 三 handler 实现 `.out`/`.err` 分离。
> 完整配置代码与常见错误见 [docs/dev-standards/02_logging_detail.md](docs/dev-standards/02_logging_detail.md)

### 2.4 中英文对照

所有输出信息必须中英文对照，`|` 分隔，中文在前。❌ 禁止 emoji（终端显示异常、不利解析）。

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

### 3.2 CLI 包装器(Click)

使用 Click 创建 CLI 包装器，**懒加载**主函数（避免导入失败）：
- `_lazy_import_main()` 延迟导入 `module.main`
- `_validate_path_exists()` 配合 `_is_help_request()` 校验输入路径
- 构造 `sys.argv` 列表 → 临时替换 → 调用 `module_main()` → `finally` 还原

> 完整 Click 模板代码见 [docs/dev-standards/01_module_template.md §2](docs/dev-standards/01_module_template.md)

### 3.3 Help 文档格式

示例只展示最简单的用法，保持简洁：

```python
def module_name(input, output_dir):
    """
    功能描述|Function description

    示例|Examples: biopytools module -i input.txt -o output_dir/
    """
```

- `示例|Examples:` 后直接跟命令，不换行；命令保持一行，不用反斜杠换行
- 只展示最常用的1个基本示例；复杂用法引导到 README/独立文档

---

## 四、命名规范

- **文件/目录**：小写+下划线 `module_name.py`、`bam_stats/`、`merqury_qv/`
- **变量/函数**：snake_case `fastq_files`、`calculate_qv()`
- **类名**：CamelCase `MerquryQVCalculator`、`ModuleLogger`
- **常量**：UPPER_SNAKE `DEFAULT_THREADS`、`MAX_DEPTH`

---

## 五、代码风格规范

### 5.1 中英文对照
所有字符串、注释、文档必须中英文对照：
```python
def calculate_qv_value():
    """计算QV值|Calculate QV value"""
```

### 5.2 禁止 emoji
```python
logger.info("开始处理|Starting processing")   # ✅
logger.info("完成|Completed")                 # ✅
# ❌ logger.info("完成✅")
```

### 5.3 大数字格式化
大于1百万的数字用 M 单位，保留2位小数：
```python
def format_number(num: int) -> str:
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
# logger.info(f"总reads数|Total reads: {format_number(10000000)}")  # 10.00M
```

---

## 六、错误处理规范

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

输入验证：收集所有错误后一次性抛出，不要逐条返回：
```python
def validate(self):
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

### 7.1 日志格式
- [ ] 点号 `.` 分隔秒和毫秒；不用方括号 `[]`；` - ` 分隔各部分
- [ ] 所有日志中英文对照；不用 emoji；时间戳后不换行

### 7.2 代码结构
- [ ] 模块化结构(config/utils/calculator/main)；dataclass 配置；标准日志管理器；CLI 与核心逻辑分离

### 7.3 命名
- [ ] 字符串中英文对照；snake_case 变量/函数；CamelCase 类名；不用 emoji

### 7.4 功能
- [ ] 输入验证；异常处理；路径存在性检查；大数字 M 单位

### 7.5 ⚠️ 命令执行日志（必须）
- [ ] **所有外部命令执行前记录完整命令（INFO 级别）**
- [ ] **禁止只记录描述不记录命令**；conda run 包装也要记录；管道用 `|` 连接
```bash
# 检查方法:应在每个 subprocess.run 前看到 logger.info(f"命令|Command: {' '.join(cmd)}")
grep -n "subprocess.run\|Popen" biopytools/module_name/*.py
```

### 7.6 ⚠️ 路径管理（关键）
- [ ] **无硬编码绝对路径**（`grep -r "/share/org/YZWL/yzwl_lixg/" biopytools/module_name/` 应返回0结果）
- [ ] 默认路径用 `~` 展开；工具路径用 `get_tool_path()`；用户输入用 `expand_path()`
- [ ] **`__post_init__` 中展开所有 `~` 路径**（关键！见 §十一）
```python
def __post_init__(self):
    self.tool_path = expand_path(self.tool_path)   # 必须
    self.genome_path = expand_path(self.genome_path)
```

**❌ 任何硬编码绝对路径都会导致代码审查失败**

---

## 八、标准代码模板

`config.py` / `utils.py` / `calculator.py` / `main.py` 完整模板已集中到参考文档：

👉 **[docs/dev-standards/01_module_template.md](docs/dev-standards/01_module_template.md)**

---

## 九、常见问题

- **Q1 为什么日志用点号不用逗号？** 点号符合小数点视觉习惯；通过 `%(msecs)03d` 显式指定。
- **Q2 中英文顺序？** 中文在前，英文在后，`|` 分隔。
- **Q3 可以用 emoji 吗？** 不可以。
- **Q4 大数字怎么格式化？** >1百万用 M 单位保留2位小数，`10.00M` 而非 `10000000`。

---

## 十、开发约束与工作流规范

### 10.1 版本控制规范|Version Control

> ⚠️ **代码提交只在本地 Mac 进行；超算上严禁任何 Git 写操作**（见文首警告）。
> 超算开发完用 `copybiopytools` 同步到本地，由本地 Claude 统一 commit + push。

commit message 格式：`<类型>(<模块>): <中文描述>`

```
类型|Type: feat 新功能 / fix Bug修复 / docs 文档 / refactor 重构 / test 测试 / chore 构建
```
```bash
git commit -m "feat(busco): 添加conda环境自动检测"
git commit -m "fix(path): 修复~路径在__post_init__中未展开的问题"
```
#### 10.1.1 Mac↔超算↔GitHub 三角工作流|Mac↔HPC↔GitHub Triangle Workflow

> 三个位置|Three locations: **GitHub**=唯一代码源；**Mac**=开发+唯一提交点；**超算**=运行+现场改码(不 commit)。

**循环A：超算写代码|Cycle A: HPC writes code**

```bash
# ① 超算改代码(只写不 commit|edit only, never commit)
# ② 同步回 Mac(执行前 Mac 必须干净|Mac must be clean first)
copybiopytools
# ③ Mac 核对 → 提交 → 推送|verify → commit → push
git diff && git commit -m 'feat(...): ...' && git push
```

**循环B：Mac 写代码|Cycle B: Mac writes code**

```bash
# ① Mac 改 → 提交 → 推送|edit → commit → push
git commit -m 'fix(...): ...' && git push
# ② 超算拉取(执行前超算必须干净|HPC must be clean first)
cd /share/org/YZWL/yzwl_lixg/software/biopytools
git pull --ff-only origin main
```

**两条铁律|Two iron rules:**

1. **同一时刻只有一边在改代码**|Only one side edits at a time — 否则同步互相覆盖
2. **跨边操作前对侧必须干净**|The other side must be clean before any cross-boundary step：
   `copybiopytools` 前 Mac 要 `git status` 干净；超算 `git pull` 前要 `git status --short` 干净（不干净先 `git stash` 或先回传 Mac）

**开工检查|Start-of-day:** Mac `git status -sb`（干净且与 GitHub 同步）；超算 `git status --short` 干净后 `git pull --ff-only origin main`。
**收工|End-of-day:** 超算改过代码就 `copybiopytools`，回 Mac 上 `git diff` → commit → push，不留跨夜未提交状态。

### 10.2 断点续传规范|Checkpoint Resume

**所有多步骤流程必须支持断点续传**——已完成的步骤重新运行时自动跳过：

```python
def _is_step_completed(self, output_file: str) -> bool:
    return Path(output_file).exists()

def run_step(self, step_name, output_file, run_func, *args, **kwargs):
    if self._is_step_completed(output_file):
        self.logger.info(f"跳过已完成步骤|Skipping completed step: {step_name}")
        return True
    self.logger.info(f"开始步骤|Starting step: {step_name}")
    return run_func(*args, **kwargs)
```

### 10.3 计算节点限制|Compute Node Restrictions

> ⚠️ **登录节点禁止运行大量计算任务！** 仅用于代码编辑、文件管理、任务提交。
> 大 CPU/内存任务必须通过作业调度系统（如 `sub`）提交到计算节点。
> 自动化测试只能测：参数解析、路径验证、小型 mock 单元测试；**大型计算测试手动在计算节点运行**。

### 10.4 临时测试/调试输出位置|Ad-hoc Test & Debug Output Location

> **原则|Principle:** 探索性/临时测试（跑个函数看输出、临时脚本验证、肉眼检查结果）**严禁在仓库当前目录生成文件**，所有产物写入 `~/tmp/`。

❌ **禁止|FORBIDDEN:** 在项目根、模块目录、或任何仓库目录下直接写测试输出
```bash
python -c "..." > out.txt        # ❌ 污染仓库 cwd
python my_scratch.py             # ❌ 若脚本把结果写到 cwd
```

✅ **必须|MUST:** 所有 ad-hoc 测试产物写到 `~/tmp/<描述性子目录>/`
```bash
mkdir -p ~/tmp/test_<module>_<purpose>
python my_scratch.py -o ~/tmp/test_<module>_<purpose>/   # ✅ 绝对路径指定输出
# 或直接在 ~/tmp 下运行脚本
cd ~/tmp/test_<module>_<purpose> && python /abs/path/to/module/main.py ...
```

**与现有规则区分|Distinguish from:**
- §11.A 正式 pytest 单测 → 仍在 `biopytools/tests/test_<module>/`（已 gitignore，本条不适用）
- §12.4 模块运行时临时文件 → 仍用 `output_dir/tmp`（真实流程产物，本条不适用）
- 本条仅约束**人/AI 手动跑的探索性测试**的产物位置

**Why:** 避免污染仓库目录、避免误提交、避免 `.gitignore` 漏网；`~/tmp` 是临时区，无需保留的可定期清理。

---

## 十一、测试规范与路径管理规范

### 11.A 测试规范|Testing Standards

**目录结构：** `biopytools/tests/test_<module>/{test_config,test_utils,test_calculator}.py` + `conftest.py`
**命名：** 文件 `test_<被测模块>.py`；函数 `test_<功能>()`；类 `Test<类名>()`

**登录节点安全：** 用 `unittest.mock` 替代真实子进程，避免触发计算。
> ⚠️ 涉及真实工具调用（BUSCO、jellyfish 等）的集成测试不能在登录节点运行。
> mock 测试示例见 [docs/dev-standards/01_module_template.md §4](docs/dev-standards/01_module_template.md)

### 11.B 路径管理规范|Path Management

#### 优先级（高→低）
1. 环境变量 `{TOOL_NAME}_PATH`
2. 用户配置文件 `~/.config/biopytools/config.yml`
3. 代码默认值（必须支持 `~` 展开）

#### 核心规则
```python
from common.paths import get_tool_path, expand_path

# ❌ 禁止:硬编码绝对路径
tool_path = "/share/org/YZWL/yzwl_lixg/.../fanc"
# ✅ 正确:路径管理系统
tool_path = get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc', 'FANC_PATH')
```

#### ⚠️ 必须展开 `~` 路径（关键）
Python **不会**自动展开字符串中的 `~`；`shutil.which("~/...")` 返回 `None`。必须在 `__post_init__` 展开：
```python
@dataclass
class ModuleConfig:
    meme_path: str = '~/miniforge3/envs/meme_v.5.5.9/bin/meme'
    def __post_init__(self):
        self.meme_path = expand_path(self.meme_path)   # ⚠️ 必须
```
- [ ] 所有 `~` 默认路径在 `__post_init__` 调 `expand_path()`
- [ ] 命令行传入的路径也用 `expand_path()` 展开
- [ ] 用统一的 `expand_path()`，不用裸 `os.path.expanduser()`

#### 路径格式
✅ `~/miniforge3/bin/tool` / `$SOFTWARE/tool/bin` / `./tools/tool`　❌ `/share/org/YZWL/...`

#### 配置文件格式 `~/.config/biopytools/config.yml`
```yaml
tools:
  fanc: ~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
  bwa: ~/.local/bin/bwa
databases:
  nr: ~/database/ncbi/nr
```

> `common/paths.py` 完整实现、config/main 用法、迁移脚本见 [docs/dev-standards/11_path_management.md](docs/dev-standards/11_path_management.md)

---

## 十二、输出目录和文件命名规范|Output Directory and File Naming

> **参考|Ref:** nf-core, GATK Best Practices, Bioconda Guidelines

### 12.1 通用原则
- **机器友好：** 禁空格；避免 `& * ( ) , @` 等 Shell 敏感字符；统一小写（Sample ID 除外）；`_` 用于文件名内部
- **人类可读：** 见名知意（拒绝 `output1`/`result.txt`）；流程步骤用数字前缀
- **可追溯：** 文件名含样本 ID、处理步骤、必要时关键参数

### 12.2 输出文件夹命名

#### 默认结构：By-Step（按步骤）
所有样本共享带序号的步骤目录，多样本用**文件名前缀** `{sample}_xxx` 区分；仅当单样本产出大量独立文件、需完全隔离时才用 By-Sample：
```
output/                          # ✅ 推荐:by-step
├── 00_pipeline_info/            # 全局元数据(P0 必须)
├── 01_qc/
│   ├── sample1_1.clean.fq.gz    # 文件名前缀区分样本
│   └── sample2_1.clean.fq.gz
├── 02_alignment/
│   └── sample1.sorted.bam
└── 99_logs/                     # 日志(P0 必须)
```
- 步骤目录格式 `[序号]_[步骤名]`，如 `01_qc/`；序号两位数
- `00_pipeline_info/`（元数据）、`99_logs/`（日志）固定在首尾；`work/` 中间文件可删
- **条件性/支持性目录不加编号**：`fastk/`（条件步骤）、`reference/`、`index/`（预生成输入）、`cache/`、`tmp/`

### 12.3 输出文件命名

**基本模式：`{Sample_ID}.{Tool/Step}.{State}.{Extension}`**（Sample_ID 必须在最前）

**后缀叠加（GATK 风格）：** `sample.bam → sample.sorted.bam → sample.sorted.markdup.bam`

| 文件类型 | 推荐命名 | 示例 |
|---------|---------|------|
| FASTQ 原始 | `{ID}_R{1/2}.fastq.gz` | `S1_R1.fastq.gz` |
| 质控报告 | `{ID}.fastqc.html` | `S1.fastqc.html` |
| 过滤数据 | `{ID}.trimmed.fq.gz` | `S1.trimmed.fq.gz` |
| Jellyfish | `{ID}.jellyfish.{jf,histo}` | `S1.jellyfish.jf` |
| GenomeScope | `{ID}.genomescope.{model.txt,linear.png}` | `S1.genomescope.model.txt` |
| Smudgeplot | `{ID}.smudgeplot.linear.png` | `S1.smudgeplot.linear.png` |
| 日志 | `{ID}.{tool}.log` | `S1.genomescope.log` |

❌ `genome_analysis.jf`（无样本名）/ `smudgeplot_smudgeplot.png`（重复前缀）/ `plot.png` / `result.txt`
✅ `R0590-6.jellyfish.jf` / `R0590-6.smudgeplot.linear.png`

### 12.4 临时文件
统一用 `output_dir/tmp` 子目录并运行结束清理（**避免超算 `/tmp` 爆满**）；禁止混入结果目录。
> `tempfile.TemporaryDirectory` 用法、清理场景、上下文管理器实现见 [docs/dev-standards/12_output_naming.md](docs/dev-standards/12_output_naming.md)

### 12.5 版本信息
- **禁止文件名含版本号**：`sample.bwa.bam` ✅　`sample_bwa_v0.7.17.bam` ❌
- 生成 `00_pipeline_info/software_versions.yml`（代码模板见参考文档）

### 12.6 检查清单
- [ ] 文件夹名无空格、小写；步骤目录用 `01_`/`02_` 数字前缀
- [ ] `00_pipeline_info/` + `99_logs/`；条件性/支持性目录不加编号
- [ ] 输出文件含样本 ID 前缀，遵循 `{Sample}.{Tool}.{State}.{Ext}`
- [ ] 扩展名标准含压缩格式（`.fastq.gz`/`.bam`）；临时文件用 `output_dir/tmp` 并清理
- [ ] 生成 `software_versions.yml`；日志统一 `99_logs/`；文件名不含版本号；避免重复前缀

> 完整结构示例（genomescope 改进对照）、版本信息代码见 [docs/dev-standards/12_output_naming.md](docs/dev-standards/12_output_naming.md)

---

## 十三、Conda环境软件调用规范|Conda Environment Software Invocation

### 13.1 问题
conda 环境中的软件（尤其 Python 包）直接调用会失败（依赖隔离、`which` 路径错乱），须用 `conda run -n <env>` 包装。

### 13.2 核心铁律（必须遵守）

#### ⚠️ 13.2.1 `conda run` 必须加 `--no-capture-output`
不加会导致 conda 缓冲所有输出，大数据（如 BWA 比对 32GB Hi-C）触发 `CondaMemoryError` / OOM。
```python
# ❌ 绝对禁止
full_cmd = ['conda', 'run', '-n', 'env', 'bwa', 'mem', ...]
# ✅ 必须加 --no-capture-output(流式输出,避免内存溢出)
full_cmd = ['conda', 'run', '-n', 'env', '--no-capture-output', 'bwa', 'mem', ...]
```

#### 🚨 13.2.2 严禁管道中使用 conda run
`conda run | conda run` 会 BrokenPipeError + 内存高 + 性能差。管道中的工具改用：方案A 直接调用并设 `LD_LIBRARY_PATH`；方案B `_build_pipeline_command` 提取实际命令；方案C `shutil.which` 自动检测。
> 三种方案完整代码见 [docs/dev-standards/13_conda_invocation.md §2](docs/dev-standards/13_conda_invocation.md)

#### 13.2.3 必须传完整路径，禁止提取命令名
```python
# ❌ os.path.basename() 提取命令名 → 丢失 /envs/ 路径 → get_conda_env 失败 → 缺 -n <env>
build_conda_command(os.path.basename(tool_path), args)
# ✅ 传完整路径 → 自动提取环境名
build_conda_command(tool_path, args)
```

### 13.3 使用方式
统一用 `build_conda_command(command, args)` 自动检测环境并包装；返回列表配合 `subprocess.run(shell=False)`。
```python
from .utils import build_conda_command
cmd = build_conda_command(self.config.busco_path, ["-i", in_fa, "-l", lineage, ...])
result = subprocess.run(cmd, shell=False, ...)
```

### 13.4 关键注意事项
1. 自动检测 conda 环境，无需用户指定
2. 向后兼容：非 conda 软件直接调用
3. `conda run` 有轻微启动开销，通常可忽略
4. 即使传完整路径，`conda run` 也能正确处理
5. 用 `CONDA_EXE` 环境变量定位 conda 基础目录

> `get_conda_env`/`build_conda_command`/`check_dependencies`/`CommandRunner` 完整实现、常见错误、测试、故障排查见 [docs/dev-standards/13_conda_invocation.md](docs/dev-standards/13_conda_invocation.md)

### 13.5 软件→环境速查表|Software→Env Quick Reference

> ⚠️ **超算上找现成软件时，先查速查表再选环境**：
> [docs/conda_env_software_map.md](docs/conda_env_software_map.md)

核心规则|Core rules:
1. **优先**使用 14 个功能域环境（align/pop/asm/hic/annot/repeat/rna/protein/phylo/pan/viz/misc/r/busco）
2. 域环境没有的软件 → 查速查表第二部分的保留独立环境（legacy 强依赖）
3. **禁止**使用 scripts/delete_list.txt 中 154 个待退役环境，新模块也不得依赖它们
4. 新模块引入新软件 → 优先并入现有域环境（配方在 envs/*.yml），禁止新建环境
5. 调用方式：`conda run -n <env> <tool> --no-capture-output`（§13.2.1）

---

## 十四、模块文档规范|Module Documentation Standards

> 📌 **每个模块必须配一份用户文档 `docs/<module>.md`，创建/修改模块时同步写好，不得事后补。**
> 受众含**无生信基础的用户**：技术定义与通俗解释双轨并行。完整模板示例见 [docs/cim.md](docs/cim.md)。

### 14.1 文档固定结构|Fixed template

`docs/<module>.md` 按以下顺序组织（`## 分析流程` 可选，其余必写）：

1. `# 标题`：标题行下第一段用**一句话人话**说明「这工具干什么、解决什么问题」
2. `## 功能概述 | Overview`：要点列表
3. `## 快速开始 | Quick Start`：只放 1 个最简命令示例（同 §3.3 规范）
4. `## 零基础概念速览 | Concepts in plain words`：术语→通俗解释对照表（比喻优先）
5. `## 输入 | Input`：文件格式与要求，附格式示例
6. `## 参数说明 | Parameters`：按功能域分组；**参数表自动生成，禁止手写**（§14.4）
7. `## 分析流程 | Pipeline`：流程示意图（复杂流程推荐）
8. `## 输出 | Output`：目录树 + 关键文件逐一说明
9. `## 结果解读 | Interpreting Results`：每个输出怎么读、好坏判据
10. `## 参数选择建议 | Parameter Guidance`：按场景给建议
11. `## 依赖 | Dependencies`
12. `## 常见问题 | FAQ`：真实踩坑点

### 14.2 通俗化写作要求|Plain-language Writing Rules

- 每个参数组前加 `**通俗理解|In plain words:**` 段落，讲清三件事：**管什么、调大/调小会怎样、什么时候才需要动**（能说「一般不用动」就说）
- 术语首次出现用生活化比喻（如 遗传距离 cM =「开车时间」，物理位置 bp =「门牌号」）
- **人话在前、术语在后**；专业定义保留给进阶用户
- ❌ 禁止只写术语定义；❌ 禁止大段英文；✅ 中文为主、关键术语带英文

### 14.3 中英文对照

标题沿用规范：中文在前、`|` 分隔英文在后（如 `## 结果解读 | Interpreting Results`）；正文 prose 中文为主，表头中英对照。
- **表格列头禁止用 `|` 分隔中英文**（会破坏 Markdown 表格解析），用 `<br>`：`| 分组<br>Group |`
- **需要被页面内锚点引用的标题必须显式声明 id**：`## 标题 | Title { #title-id }`（attr_list）。MkDocs slugify 会剥掉中文，自动锚点不可靠，卡片/目录跳转会变死链

### 14.4 参数表自动生成（禁止手写）|Auto-generated parameter tables

- 「参数说明」的参数表**从 click/argparse 定义自动提取生成**（生成器按需构建），**严禁手写**——手写必然与代码漂移
- 人只写：参数分组逻辑、每组「通俗理解」段落、结果解读、参数建议、FAQ
- 参数有变 → 重跑生成器，文档与代码永远同步

### 14.5 检查清单|Checklist

- [ ] 模块目录存在对应 `docs/<module>.md`，结构齐全（§14.1）
- [ ] 每个参数组有「通俗理解|In plain words」段落；术语有比喻
- [ ] 参数表来自生成器（无手写参数表）
- [ ] 输出文件逐一说明 + 好坏判据（结果解读）
- [ ] 快速开始只有 1 个最简示例（同 §3.3）
- [ ] **审查口径：新模块/新参数无对应文档更新 → 审查不通过**（与硬编码绝对路径同级）

---

## 📚 详细参考文档

> 以下文档**按需读取**（任务相关时再 Read，不主动全读）。核心规则已在本文件全部覆盖。

| 触发场景|Trigger | 参考文档|Reference |
|---------|---------|---------|
| 新建模块、写 config/utils/calculator/main | [docs/dev-standards/01_module_template.md](docs/dev-standards/01_module_template.md) |
| 配置超算 stdout/stderr/file 日志分离 | [docs/dev-standards/02_logging_detail.md](docs/dev-standards/02_logging_detail.md) |
| 改路径管理、迁移旧代码、查 paths.py 实现 | [docs/dev-standards/11_path_management.md](docs/dev-standards/11_path_management.md) |
| 设计输出目录结构、写版本信息 yml | [docs/dev-standards/12_output_naming.md](docs/dev-standards/12_output_naming.md) |
| 调用 conda 环境软件、排查 conda 命令 | [docs/dev-standards/13_conda_invocation.md](docs/dev-standards/13_conda_invocation.md) |
| 写模块用户文档、参数通俗化写法（§14） | [docs/cim.md](docs/cim.md)（cim 模板示例） |
| 查某软件装在哪个 conda 环境（超算找现成软件） | [docs/conda_env_software_map.md](docs/conda_env_software_map.md) |

---

## 版本历史

| 版本 | 日期 | 主要变更|Major Changes |
|------|------|----------|
| 2.21 | 2026-08-17 | 新增 §14「模块文档规范」：每模块必须配 `docs/<module>.md`（固定12节模板+通俗化写作要求+参数表自动生成禁止手写，审查不通过条款）；模板示例 docs/cim.md |
| 2.20 | 2026-08-16 | 新增 §10.1.1「Mac↔超算↔GitHub 三角工作流」：两循环+两铁律(单边改码/跨边先清对侧) |
| 2.19 | 2026-08-16 | 新增 §13.5「软件→环境速查表」引用（docs/conda_env_software_map.md）：超算找现成软件优先 14 域环境、禁用手工废弃环境 |
| 2.18 | 2026-07-28 | 新增 §10.4「临时测试/调试输出位置」：探索性/ad-hoc 测试产物严禁写入仓库 cwd，统一放 `~/tmp/<描述性子目录>/`；与 §11.A 正式单测(tests/)、§12.4 流程临时文件(output_dir/tmp)明确区分 |
| 2.17 | 2026-07-25 | **CLAUDE.md 瘦身拆分**：核心规则+检查清单保留(76KB→约23KB)，完整代码模板/paths.py实现/命名示例/conda故障排查下沉到 `docs/dev-standards/` 五个按需参考文档(01_module_template/02_logging_detail/11_path_management/12_output_naming/13_conda_invocation)，文末加「📚详细参考文档」触发式索引；规则零丢失，仅外移重内容 |
| 2.16 | 2026-07-24 | §12.2.1 目录结构默认改为 **by-step**（多样本共享步骤目录 + 文件名前缀 `{sample}_xxx` 区分），by-sample 降为可选；§12.6.1 加注释标注既有 genomescope 结构为新模块前的示例 |
| 2.15 | 2026-07-23 | 临时目录统一改造：所有模块临时文件/目录从系统 `/tmp` 改为 `output_dir/tmp` 子目录并运行结束清理，消除超算 `/tmp` 爆满风险 |
| 2.14 | 2026-06-24 | 新增"禁止在超算上执行 Git 提交"规范（顶部警告 + 修订 10.1）；配套 `copybiopytools` 增加 `--exclude='.git/'` |
| 2.13 | 2026-04-09 | 文档质量改进：修正版本号不一致；修复测试断言缺少--no-capture-output；裸except改except Exception；修复示例硬编码绝对路径；补全版本信息字段；添加目录导航 |
| 2.12 | 2026-03-17 | 完善 Conda 调用规范：新增13.2.0强制 `--no-capture-output`，避免 CondaMemoryError |
| 2.11 | 2026-03-17 | 新增命令执行日志规范(2.2.1)：所有外部命令执行前记录完整命令到 INFO；代码审查检查点(7.5) |
| 2.10 | 2026-03-15 | 完善 Conda 规范：新增13.6常见错误、13.7测试、13.8故障排查；禁止 os.path.basename() 提取命令名 |
| 2.9 | 2026-03-05 | 修复 CommandRunner shell=True Bug；补全 common/paths.py；新增第十章开发约束、测试规范章节 |
| 2.8 | 2026-03-05 | 添加 Conda 环境软件调用规范（自动检测并使用 conda run）|
| 2.7 | 2026-03-05 | 添加超算日志系统分离规范（INFO→.out, WARNING→.err）|
| 2.6 | 2026-03-05 | 澄清多样本 vs 单样本目录结构原则|
| 2.5 | 2026-03-05 | 添加输出目录和文件命名规范（参考 nf-core、GATK 最佳实践）|
| 2.4 | 2026-03-02 | 添加 ~ 路径展开的强制要求和验证方法|
| 2.3 | 2026-03-02 | 添加文档开头的路径管理警告和代码审查检查项|
| 2.2 | 2026-03-02 | 添加路径管理规范（避免硬编码绝对路径）|
| 2.1 | 2026-01-05 | 添加 Help 文档格式规范|
| 2.0 | 2026-01-04 | 更新日志格式规范；添加模块化结构规范；强化中英文对照；大数字格式化|
| 1.0 | 2024-12-19 | 初始版本|Initial version |
