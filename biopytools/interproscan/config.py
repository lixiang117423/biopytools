"""
InterProScan注释配置管理模块|InterProScan Annotation Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class InterProScanConfig:
    """InterProScan注释配置类|InterProScan Annotation Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_prefix: str

    # 软件路径|Software path
    interproscan_path: str = '~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh'

    # 输出格式|Output format (支持多个格式逗号分隔|Support multiple formats comma-separated)
    output_format: str = 'TSV,XML'  # 默认同时输出TSV和XML用于解析GO/Pathway

    # 处理参数|Processing parameters
    threads: int = 24
    disable_precalc: bool = True  # 禁用预计算查找服务|Disable precalc lookup service
    goterms: bool = False  # 获取GO术语|Get GO terms
    pathways: bool = False  # 获取pathway信息|Get pathway information

    # 应用列表|Applications list
    applications: Optional[str] = None  # 例如: "Pfam,SMART,Gene3D" | e.g., "Pfam,SMART,Gene3D"

    # 其他选项|Other options
    temp_dir: Optional[str] = None  # 临时目录|Temporary directory

    # 结果整理选项|Result formatting options
    generate_report: bool = True  # 是否生成整理后的报告|Whether to generate formatted report
    report_format: str = 'both'  # 报告格式: 'excel', 'csv', 'both' | Report format (默认both，同时生成Excel和CSV|Default both, generate both Excel and CSV)
    include_summary: bool = True  # 是否包含汇总统计表|Whether to include summary statistics

    # GO数据库选项|GO database options
    go_database_path: Optional[str] = None  # GO术语JSON数据库路径（可选，默认使用内置数据库）|GO term JSON database path (optional, built-in database used by default)

    # Python解释器选项|Python interpreter options
    python_path: Optional[str] = None  # Python解释器路径（可选，用于指定兼容的Python版本）|Python interpreter path (optional, for specifying compatible Python version)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        import sys

        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        # 使用expand_path展开~符号|Use expand_path to expand ~ symbol
        self.interproscan_path = os.path.normpath(expand_path(self.interproscan_path))

        # 设置输出目录|Setup output directory
        output_dir = os.path.dirname(os.path.abspath(self.output_prefix))
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 设置临时目录|Setup temp directory
        if self.temp_dir:
            # 使用expand_path展开~符号|Use expand_path to expand ~ symbol
            self.temp_dir = os.path.normpath(expand_path(self.temp_dir))
            # 移除路径末尾可能存在的斜杠|Remove trailing slash if present
            self.temp_dir = self.temp_dir.rstrip(os.sep)
            Path(self.temp_dir).mkdir(parents=True, exist_ok=True)

        # 检测并设置Python路径|Detect and setup Python path
        if self.python_path:
            self.python_path = expand_path(self.python_path)
        else:
            # 检查系统Python版本|Check system Python version
            py_version = sys.version_info
            # InterProScan 5.70-102.0+ 要求 Python 3.8-3.11|Requires Python 3.8-3.11
            if py_version.major == 3 and py_version.minor >= 12:
                # Python 3.12+ 不兼容，尝试查找兼容版本|Python 3.12+ incompatible, try to find compatible version
                self.logger = None  # 延迟初始化日志|Delayed logger init
                self._find_compatible_python()
            else:
                self.python_path = sys.executable

    def _find_compatible_python(self):
        """查找兼容的Python版本|Find compatible Python version"""
        import shutil
        compatible_versions = ['3.11', '3.10', '3.9', '3.8']

        for version in compatible_versions:
            python_cmd = f'python{version}'
            python_path = shutil.which(python_cmd)
            if python_path:
                self.python_path = python_path
                print(f"警告|WARNING: 系统 Python 版本为 3.12+，InterProScan 不兼容。")
                print(f"使用兼容版本|Using compatible version: {python_path}")
                return

        # 未找到兼容版本，使用系统默认|No compatible version found, use system default
        self.python_path = 'python3'
        print("警告|WARNING: 未找到 Python 3.8-3.11，InterProScan 可能运行失败")
        print("建议|RECOMMENDED: 安装兼容的Python版本或使用 --python-path 参数指定")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        # 检查InterProScan路径|Check InterProScan path
        if not os.path.exists(self.interproscan_path):
            errors.append(f"InterProScan路径不存在|InterProScan path does not exist: {self.interproscan_path}")

        # 检查输出格式|Check output format (支持多格式|Support multiple formats)
        valid_formats = ['TSV', 'GFF3', 'XML', 'HTML', 'JSON', 'TXT']
        formats = [f.strip().upper() for f in self.output_format.split(',')]
        for fmt in formats:
            if fmt not in valid_formats:
                errors.append(
                    f"无效的输出格式|Invalid output format: {fmt} "
                    f"(支持的格式|Supported formats: {', '.join(valid_formats)})"
                )

        # 检查报告格式|Check report format
        if self.report_format.lower() not in ['excel', 'csv', 'both']:
            errors.append(
                f"无效的报告格式|Invalid report format: {self.report_format} "
                f"(支持的格式|Supported formats: excel, csv, both)"
            )

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_output_file(self) -> str:
        """获取输出文件路径|Get output file path (主格式|Main format)"""
        ext_map = {
            'TSV': '.tsv',
            'GFF3': '.gff3',
            'XML': '.xml',
            'HTML': '.html',
            'JSON': '.json',
            'TXT': '.txt'
        }
        # 获取第一个格式作为主输出|Get first format as main output
        main_format = self.output_format.split(',')[0].strip().upper()
        ext = ext_map.get(main_format, '.tsv')
        return f"{self.output_prefix}{ext}"

    def get_output_files(self) -> dict:
        """获取所有输出文件路径|Get all output file paths"""
        ext_map = {
            'TSV': '.tsv',
            'GFF3': '.gff3',
            'XML': '.xml',
            'HTML': '.html',
            'JSON': '.json',
            'TXT': '.txt'
        }
        formats = [f.strip().upper() for f in self.output_format.split(',')]
        return {fmt: f"{self.output_prefix}{ext_map.get(fmt, '.txt')}" for fmt in formats}

    def has_xml_output(self) -> bool:
        """检查是否包含XML输出|Check if XML output is included"""
        return 'XML' in [f.strip().upper() for f in self.output_format.split(',')]
