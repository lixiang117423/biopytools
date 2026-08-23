"""xpclr 配置|xpclr configuration.

所有 ~ 路径在 __post_init__ 中展开(§十一);工具路径走 get_tool_path 优先级链。
|All ~ paths expanded in __post_init__; tool path via get_tool_path chain.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class XpclrConfig:
    """XP-CLR 跨群体扫描配置|XP-CLR cross-population scan configuration."""

    # 输入|Input
    input_vcf: str = ""
    samples_a: str = ""           # 群体A样本列表文件(每行一个ID)|Pop A sample list file
    samples_b: str = ""           # 群体B样本列表文件|Pop B sample list file
    output_dir: str = "xpclr_out"
    label: str = "popA_vs_popB"   # 结果文件前缀|Result file prefix
    chroms: Optional[List[str]] = None  # None=VCF header 全部 contig|None=all contigs
    # 窗口参数(默认=工具默认)|Window params (tool defaults)
    size: int = 20000             # 窗口大小(bp)|Window size (bp)
    step: int = 20000             # 滑窗步长(bp)|Step size (bp)
    maxsnps: int = 200            # 窗口最大SNP数|Max SNPs per window
    minsnps: int = 10             # 窗口最小SNP数|Min SNPs per window
    ld: float = 0.95              # LD 加权截断|LD cutoff for weighting
    phased: bool = False          # 数据已 phased|Data phased
    rrate: float = 1e-8           # 每碱基重组率|Recombination rate per base
    # 汇总|Summary
    top_n: int = 50               # Top 候选窗口数|Top candidate windows
    # 工具路径|Tool path
    xpclr_path: str = "~/miniforge3/envs/selective_sweep/bin/xpclr"
    # 后端|Backend
    backend: str = "xpclrs"          # xpclrs=Rust高速版(默认)|xpclrs=Rust fast (default); xpclr=python
    xpclrs_path: str = "~/software/xpclrs/bin/xpclrs"   # Rust 二进制|Rust binary
    threads: int = 12                # 线程数(仅 xpclrs 后端)|Threads (xpclrs backend only)

    def __post_init__(self):
        self.input_vcf = expand_path(self.input_vcf)
        self.samples_a = expand_path(self.samples_a)
        self.samples_b = expand_path(self.samples_b)
        # expand_path 对裸名(无斜杠无点号)不转绝对路径,统一 abspath 兜底
        # |expand_path skips bare names; force abspath for a stable output root
        self.output_dir = os.path.abspath(expand_path(self.output_dir))
        self.xpclr_path = get_tool_path("xpclr", self.xpclr_path, "XPCLR_PATH")
        self.xpclrs_path = get_tool_path("xpclrs", self.xpclrs_path, "XPCLRS_PATH")

    def validate(self) -> bool:
        """收集全部错误后一次性抛出|Collect all errors then raise once (§六)."""
        errors = []
        if not self.input_vcf:
            errors.append("必须提供输入VCF|Input VCF required")
        elif not os.path.isfile(self.input_vcf):
            errors.append(f"输入VCF不存在|Input VCF not found: {self.input_vcf}")
        elif not (os.path.isfile(self.input_vcf + ".tbi")
                  or os.path.isfile(self.input_vcf + ".csi")):
            errors.append(
                f"VCF缺少tabix索引(.tbi/.csi)|VCF missing tabix index: {self.input_vcf}.tbi")
        # xpclr 把含逗号的 --samplesA/B 输入当内联列表,路径含逗号会被误解析
        # |xpclr treats comma-containing input as inline list; comma in path misparsed
        for name, path in (("samples_a", self.samples_a), ("samples_b", self.samples_b)):
            if not path:
                errors.append(f"必须提供{name}样本列表|{name} sample list required")
            elif not os.path.isfile(path):
                errors.append(f"样本列表不存在|Sample list not found [{name}]: {path}")
            elif "," in path:
                errors.append(
                    f"样本列表路径不能含逗号|Sample list path must not contain comma "
                    f"[{name}]: {path}")
        if self.size <= 0 or self.step <= 0:
            errors.append(f"size/step 必须为正|size/step must be positive: {self.size}/{self.step}")
        if self.minsnps < 2:
            errors.append(f"minsnps 必须>=2(工具硬限)|minsnps must be >= 2: {self.minsnps}")
        if self.maxsnps < self.minsnps:
            errors.append(f"maxsnps 必须>=minsnps|maxsnps must be >= minsnps: {self.maxsnps}")
        if not 0 < self.ld <= 1:
            errors.append(f"ld 必须在(0,1]|ld must be in (0,1]: {self.ld}")
        if self.rrate <= 0:
            errors.append(f"rrate 必须为正|rrate must be positive: {self.rrate}")
        if self.top_n <= 0:
            errors.append(f"top_n 必须为正|top_n must be positive: {self.top_n}")
        if not self.label or "/" in self.label:
            errors.append(f"label 不能为空且不能含斜杠|label must be non-empty, no slashes: {self.label!r}")
        if self.backend not in ("xpclrs", "xpclr"):
            errors.append(f"backend 必须是 xpclrs 或 xpclr|backend must be xpclrs or xpclr: {self.backend!r}")
        if self.threads <= 0:
            errors.append(f"threads 必须为正|threads must be positive: {self.threads}")
        # 只校验当前后端工具存在|validate only the active backend's tool
        tool = self.xpclrs_path if self.backend == "xpclrs" else self.xpclr_path
        if not os.path.isfile(tool):
            errors.append(
                f"{self.backend} 可执行不存在|{self.backend} executable not found: {tool}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
