# BioPyTools 路径配置指南|BioPyTools Path Configuration Guide

## 概述|Overview

BioPyTools 支持多种方式配置工具路径，按优先级排序：
|BioPyTools supports multiple ways to configure tool paths, in priority order:

1. **环境变量|Environment Variables** (最高优先级|Highest priority)
2. **用户配置文件|User Config File** (`~/.config/biopytools/config.yml`)
3. **代码默认值|Code Defaults** (支持~展开|Supports ~ expansion)

---

## 方法1: 环境变量|Method 1: Environment Variables

### 临时设置（当前会话）|Temporary (current session)

```bash
export FANC_PATH=~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
export BWA_PATH=~/.local/bin/bwa
export SAMTOOLS_PATH=~/.local/bin/samtools
```

### 永久设置（推荐）|Permanent (recommended)

在 `~/.bashrc` 或 `~/.zshrc` 中添加：
|Add to `~/.bashrc` or `~/.zshrc`:

```bash
# BioPyTools 工具路径|BioPyTools tool paths
export FANC_PATH=~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
export BWA_PATH=~/.local/bin/bwa
export SAMTOOLS_PATH=~/.local/bin/samtools
export AUGUSTUS_PATH=~/miniforge3/envs/Augustus_v.3.5.0/bin
# ... 添加更多工具|Add more tools
```

然后执行|Then run:
```bash
source ~/.bashrc  # 或|or source ~/.zshrc
```

---

## 方法2: 配置文件|Method 2: Config File

### 创建配置文件|Create config file

```bash
# 创建配置目录|Create config directory
mkdir -p ~/.config/biopytools

# 复制模板|Copy template
cp ~/software/biopytools/biopytools/common/config_template.yml \
   ~/.config/biopytools/config.yml

# 编辑配置|Edit config
vim ~/.config/biopytools/config.yml
```

### 配置文件示例|Config file example

```yaml
tools:
  fanc: ~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
  bwa: ~/.local/bin/bwa
  samtools: ~/.local/bin/samtools
  augustus: ~/miniforge3/envs/Augustus_v.3.5.0/bin
  # ... 更多工具|more tools

databases:
  nr: ~/database/ncbi/nr
  swissprot: ~/database/uniprot/swissprot
```

**特性|Features**:
- 支持~自动展开|Auto ~ expansion
- 支持环境变量|Support environment variables: `$SOFTWARE/fanc/bin/fanc`
- 注释行以#开头|Comments start with #

---

## 方法3: 使用默认值|Method 3: Use Defaults

如果未设置环境变量和配置文件，BioPyTools使用代码中的默认路径：
|If no env var or config, BioPyTools uses built-in defaults:

```python
# 默认路径支持~展开|Default paths support ~ expansion
fanc_path: str = "~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc"  # ✅ 正确|Correct
fanc_path: str = "/share/org/YZWL/yzwl_lixg/miniforge3/..."  # ❌ 错误|Incorrect
```

---

## 快速开始|Quick Start

### 1. 为当前用户配置|Configure for current user

```bash
# 方式A: 使用配置文件（推荐）|Option A: Use config file (recommended)
mkdir -p ~/.config/biopytools
cat > ~/.config/biopytools/config.yml << 'EOF'
tools:
  fanc: ~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc
  bwa: ~/.local/bin/bwa
  samtools: ~/.local/bin/samtools
EOF

# 方式B: 使用环境变量|Option B: Use environment variables
echo 'export FANC_PATH=~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc' >> ~/.bashrc
echo 'export BWA_PATH=~/.local/bin/bwa' >> ~/.bashrc
echo 'export SAMTOOLS_PATH=~/.local/bin/samtools' >> ~/.bashrc
source ~/.bashrc
```

### 2. 验证配置|Verify configuration

```python
from biopytools.common.paths import get_tool_path

# 测试路径|Test path
fanc_path = get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc', 'FANC_PATH')
print(f"FAN-C路径|FAN-C path: {fanc_path}")
```

### 3. 检查工具是否可用|Check if tool is available

```python
from biopytools.common.paths import validate_tool_path

# 验证工具路径|Validate tool path
try:
    validate_tool_path('~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc', 'FAN-C')
    print("FAN-C可用|FAN-C is available")
except ValueError as e:
    print(f"错误|Error: {e}")
```

---

## 常见问题|FAQ

### Q1: 我的工具路径和默认值不同怎么办？

**A**: 有3种方式|Three ways:
1. 设置环境变量（优先级最高）|Set environment variable (highest priority)
2. 编辑配置文件|Edit config file
3. 在命令行参数中指定|Specify in command line arguments

### Q2: 如何在团队中共享配置？

**A**:
- 创建团队配置模板|Create team config template
- 每个成员根据自己路径修改|Each member modifies for their paths
- 可以用环境变量指向共享位置|Use env vars to point to shared locations

### Q3: 配置文件支持相对路径吗？

**A**: 支持！以下路径都有效|Yes! All these are valid:
- `~/miniforge3/bin/fanc` → 展开为用户家目录|Expands to user home
- `$SOFTWARE/fanc/bin/fanc` → 展开环境变量|Expands env var
- `./tools/fanc` → 相对于当前目录|Relative to current dir

---

## 开发者指南|Developer Guide

### 在代码中使用|Use in code

```python
from dataclasses import dataclass, field
from biopytools.common.paths import get_tool_path

@dataclass
class MyConfig:
    # 使用路径工具|Use path utility
    fanc_path: str = field(
        default_factory=lambda: get_tool_path(
            'fanc',                                    # 工具名|Tool name
            '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc',  # 默认值|Default
            'FANC_PATH'                                # 环境变量|Env var (可选|optional)
        )
    )

    def __post_init__(self):
        # 路径已经自动展开|Path already expanded
        self.fanc_path = expand_path(self.fanc_path)
```

### 环境变量命名规范|Environment variable naming

- 格式|Format: `{TOOL_NAME}_PATH` (全大写|UPPERCASE)
- 示例|Examples: `FANC_PATH`, `BWA_PATH`, `SAMTOOLS_PATH`

---

## 更新日志|Changelog

| 版本|Version | 日期|Date | 说明|Description |
|------|------|------|------|
| 1.0 | 2026-03-02 | 初始版本|Initial version |
