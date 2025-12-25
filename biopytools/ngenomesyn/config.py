"""
配置管理模块 | Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any

@dataclass
class NGenomeSynConfig:
    """NGenomeSyn配置类 | NGenomeSyn Configuration Class"""

    # 必需参数 | Required parameters
    sample_map: Optional[str] = None
    config_file: Optional[str] = None
    output_dir: str = "./ngenomesyn_output"

    # 比对参数 | Alignment parameters
    aligner: str = "minimap2"
    threads: int = 16
    min_length: int = 5000

    # 可视化参数 | Visualization parameters
    output_prefix: str = "ngenomesyn"
    output_formats: List[str] = None
    ngenomesyn_bin: Optional[str] = None
    convert_bin: Optional[str] = None

    # 高级参数 | Advanced parameters
    chromosomes: Optional[List[int]] = None
    _chromosome_str: Optional[str] = None

    # Minimap2参数
    minimap_preset: str = "asm5"

    # MUMmer参数
    mummer_match_type: str = "mumreference"
    mummer_min_match: int = 20
    mummer_min_cluster: int = 65
    mummer_max_gap: int = 90

    # SyRI参数
    use_syri: bool = False
    syri_bin: Optional[str] = None

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 设置默认输出格式 | Set default output formats
        if self.output_formats is None:
            self.output_formats = ["svg", "png"]

        # 解析染色体参数 | Parse chromosome parameters
        if self._chromosome_str:
            self.chromosomes = self._parse_chromosome_string(self._chromosome_str)

        # 自动查找NGenomeSyn二进制文件 | Auto-find NGenomeSyn binary
        if self.ngenomesyn_bin is None:
            self.ngenomesyn_bin = self._find_ngenomesyn_bin()

        # 自动查找convert命令 | Auto-find convert command
        if self.convert_bin is None and "png" in self.output_formats:
            self.convert_bin = self._find_convert_command()

    def _parse_chromosome_string(self, chr_str: str) -> List[int]:
        """解析染色体字符串 | Parse chromosome string

        支持格式:
        - "1,2,3" : 指定第1、2、3号染色体
        - "1-5" : 指定第1到5号染色体
        - "1,3-5,7" : 混合格式
        """
        chromosomes = []
        for part in chr_str.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    chromosomes.extend(range(start, end + 1))
                except ValueError:
                    raise ValueError(f"无效的染色体范围格式: {part}")
            else:
                try:
                    chromosomes.append(int(part))
                except ValueError:
                    raise ValueError(f"无效的染色体编号: {part}")

        # 去重、排序并验证
        chromosomes = sorted(set(chromosomes))
        if any(c <= 0 for c in chromosomes):
            raise ValueError("染色体编号必须为正整数")

        return chromosomes

    def _find_ngenomesyn_bin(self) -> Optional[str]:
        """查找NGenomeSyn二进制文件 | Find NGenomeSyn binary"""
        import shutil

        # 首先尝试从PATH查找 | Try to find from PATH first
        ngenomesyn_path = shutil.which("NGenomeSyn")
        if ngenomesyn_path:
            return ngenomesyn_path

        # 尝试常见安装路径 | Try common installation paths
        possible_paths = [
            "/share/org/YZWL/yzwl_lixg/tmp/NGenomeSyn-main/bin/NGenomeSyn",
            "/usr/local/bin/NGenomeSyn",
            "./bin/NGenomeSyn",
        ]

        for path in possible_paths:
            if os.path.exists(path) and os.access(path, os.X_OK):
                return path

        return None

    def _find_convert_command(self) -> Optional[str]:
        """查找ImageMagick convert命令 | Find ImageMagick convert command"""
        import shutil
        convert_path = shutil.which("convert")
        return convert_path

    def get_chromosome_suffix(self) -> str:
        """获取染色体后缀用于文件命名 | Get chromosome suffix for file naming"""
        if self.chromosomes:
            return "_chr" + "_".join(map(str, self.chromosomes))
        return ""

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查输入文件 | Check input files
        if not self.sample_map and not self.config_file:
            errors.append("必须提供sample_map或config_file")

        if self.sample_map and not os.path.exists(self.sample_map):
            errors.append(f"Sample map文件不存在: {self.sample_map}")

        if self.config_file and not os.path.exists(self.config_file):
            errors.append(f"配置文件不存在: {self.config_file}")

        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append("线程数必须为正整数")

        if self.min_length <= 0:
            errors.append("最小长度必须为正整数")

        # 检查比对器支持 | Check aligner support
        supported_aligners = ["minimap2", "mummer"]
        if self.aligner not in supported_aligners:
            errors.append(f"不支持的比对器: {self.aligner}，支持的比对器: {', '.join(supported_aligners)}")

        # 检查染色体参数 | Check chromosome parameters
        if self.chromosomes and len(self.chromosomes) == 0:
            errors.append("染色体列表不能为空")

        # 检查NGenomeSyn二进制文件 | Check NGenomeSyn binary
        if not self.ngenomesyn_bin:
            errors.append("找不到NGenomeSyn二进制文件，请安装或指定路径")
        elif not os.path.exists(self.ngenomesyn_bin):
            errors.append(f"NGenomeSyn二进制文件不存在: {self.ngenomesyn_bin}")

        # 检查convert命令（如果需要PNG输出）| Check convert command (if PNG output needed)
        if "png" in self.output_formats and not self.convert_bin:
            errors.append("PNG输出需要ImageMagick convert命令")

        if errors:
            raise ValueError("\n".join(errors))

        return True
