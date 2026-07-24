# 模块结构与代码模板|Module Structure & Code Templates

> 📖 本文件为 CLAUDE.md 的按需参考文档,**新建模块时读取**。核心规则见 CLAUDE.md §一、§三、§八。
> On-demand reference for CLAUDE.md; read when creating a new module.

---

## 1. 标准模块目录结构|Standard Module Directory Structure

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

### 各文件职责说明|File Responsibilities

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

## 2. CLI包装器规范(Click)|CLI Wrapper (Click)

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

---

## 3. 标准代码模板(精简版)|Standard Templates (Condensed)

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

## 4. 登录节点安全测试示例|Login-Node-Safe Test Example

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
