# 输出命名完整示例|Output Naming Complete Examples

> 📖 本文件为 CLAUDE.md 的按需参考文档,**设计输出目录结构时读取**。核心命名原则与模式表见 CLAUDE.md §十二。
> On-demand reference for CLAUDE.md; read when designing a module's output directory layout.

---

## 1. 临时文件管理|Temporary File Management

### 1.1 临时文件放置|Placement

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

### 1.2 必须清理的场景|Must Cleanup Scenarios

| 场景|Scenario | 临时文件|Temporary Files | 清理时机|Cleanup Timing |
|---------|-----------|---------|----------------|---------|---------------|
| FastK输入|FastK input | 解压的 `.fq` 文件|Decompressed `.fq` | FastK完成后|After FastK |
| 压缩检查|Compression check | 解压用于检查的文件|Files for checking | 检查完成后|After check |
| 中间步骤|Intermediate | 流程中间文件|Pipeline intermediates | 最终结果生成后|After final output |

### 1.3 实现示例|Implementation Example

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

## 2. 版本信息记录|Version Information

### 2.1 software_versions.yml 格式|Format

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

### 2.2 生成版本信息的代码示例|Code Example

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

## 3. 完整输出结构示例|Complete Output Structure

> 注|Note: 以下为 genomescope **既有模块**的实际结构（by-sample，样本目录套步骤）；**新模块默认按 CLAUDE.md §12.2.1 用 by-step**（步骤目录共享、文件名带样本前缀）。此处保留以展示完整的步骤目录与 `{Sample}.{Tool}.{State}.{Ext}` 文件命名规范。

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

### genomescope 模块改进对照表|genomescope Module Improvement

| 组件<br>Component | 当前<br>Current | 改进<br>Improved | 说明<br>Notes |
|---|---|---|---|
| **目录结构**<br>Directory | `jellyfish/` | `01_jellyfish/` | 添加序号前缀<br>Add number prefix |
| | | `genomescope_output/` | `02_genomescope/` | 简化+序号<br>Simplify + number |
| | | `smudgeplot_output/` | `03_smudgeplot/` | 简化+序号<br>Simplify + number |
| **文件命名**<br>File Names | `genome_analysis.jf` | `R0590-6.jellyfish.jf` | 样本名.工具名<br>Sample.tool |
| | | `genome_analysis.histo` | `R0590-6.jellyfish.histo` | 样本名.工具名<br>Sample.tool |
| | | `model.txt` | `R0590-6.genomescope.model.txt` | 明确文件类型<br>Specify file type |
| | | `plot.png` | `R0590-6.genomescope.linear.png` | 明确比例尺<br>Specify scale |
| | | `smudgeplot_smudgeplot.png` | `R0590-6.smudgeplot.linear.png` | 去除重复<br>Remove duplicate |
| **临时文件**<br>Temp Files | `fastk/*.fq` (68GB) | `<output>/tmp/*.fq` (运行结束清理) | output_dir/tmp子目录<br>Use output/tmp |
| **版本信息**<br>Version Info | ❌ 无<br>Missing | `software_versions.yml` | 添加版本文件<br>Add version file |
| **日志管理**<br>Log Management | 散落各处<br>Scattered | `99_logs/pipeline.log` | 集中管理<br>Centralize |
