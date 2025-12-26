#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GEMMA GWAS Configuration
配置管理模块
Author: BioTools Development Team
Version: 1.0.0
"""

import os
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class PLINKQCConfig:
    """PLINK质控参数配置"""
    maf: float = 0.05          # 最小等位基因频率
    geno: float = 0.1          # SNP最大缺失率
    mind: float = 0.1          # 样本最大缺失率
    hwe: float = 1e-6          # Hardy-Weinberg平衡p值阈值
    enable: bool = True        # 是否启用质控

    def __post_init__(self):
        """验证参数范围"""
        if not 0 <= self.maf <= 1:
            raise ValueError(f"MAF must be between 0 and 1, got {self.maf}")
        if not 0 <= self.geno <= 1:
            raise ValueError(f"GENO must be between 0 and 1, got {self.geno}")
        if not 0 <= self.mind <= 1:
            raise ValueError(f"MIND must be between 0 and 1, got {self.mind}")
        if not 0 <= self.hwe <= 1:
            raise ValueError(f"HWE must be between 0 and 1, got {self.hwe}")


@dataclass
class GEMMAConfig:
    """GEMMA分析参数配置"""
    lmm_method: int = 4       # LMM检验方法: 1=Wald, 2=LRT, 3=Score, 4=all
    gk_method: int = 1        # 亲缘关系矩阵方法: 1=centered, 2=standardized
    miss: float = 0.05        # GEMMA的缺失率阈值
    maf: float = 0.01         # GEMMA的MAF阈值
    notsnp: bool = False      # 是否跳过SNP输出

    def __post_init__(self):
        """验证参数范围"""
        if self.lmm_method not in [1, 2, 3, 4]:
            raise ValueError(f"LMM method must be 1-4, got {self.lmm_method}")
        if self.gk_method not in [1, 2]:
            raise ValueError(f"GK method must be 1 or 2, got {self.gk_method}")
        if not 0 <= self.miss <= 1:
            raise ValueError(f"MISS must be between 0 and 1, got {self.miss}")
        if not 0 <= self.maf <= 1:
            raise ValueError(f"MAF must be between 0 and 1, got {self.maf}")


@dataclass
class AnalysisConfig:
    """整体分析配置"""
    # 输入文件
    vcf: str = ""
    pheno: str = ""

    # 输出目录
    outdir: str = ""

    # PCA参数
    n_pca: int = 10

    # 线程数
    threads: int = 12

    # GEMMA程序路径
    gemma_path: str = "/share/org/YZWL/yzwl_lixg/.local/bin/gemma"

    # 子配置
    plink_qc: PLINKQCConfig = field(default_factory=PLINKQCConfig)
    gemma: GEMMAConfig = field(default_factory=GEMMAConfig)

    def __post_init__(self):
        """验证配置"""
        if self.n_pca <= 0:
            raise ValueError(f"Number of PCA must be positive, got {self.n_pca}")
        if self.threads <= 0:
            raise ValueError(f"Threads must be positive, got {self.threads}")

    def validate_inputs(self) -> tuple[bool, str]:
        """
        验证输入文件

        Returns:
            tuple: (是否有效, 错误消息)
        """
        if not self.vcf:
            return False, "VCF file path is required"

        if not self.pheno:
            return False, "Phenotype file path is required"

        if not os.path.exists(self.vcf):
            return False, f"VCF file not found: {self.vcf}"

        if not os.path.exists(self.pheno):
            return False, f"Phenotype file not found: {self.pheno}"

        if not os.path.exists(self.gemma_path):
            return False, f"GEMMA program not found: {self.gemma_path}"

        return True, ""

    def get_plink_qc_args(self) -> list[str]:
        """获取PLINK质控参数列表"""
        if not self.plink_qc.enable:
            return []

        return [
            f"--maf {self.plink_qc.maf}",
            f"--geno {self.plink_qc.geno}",
            f"--mind {self.plink_qc.mind}",
            f"--hwe {self.plink_qc.hwe}"
        ]

    def get_gemma_qc_args(self) -> list[str]:
        """获取GEMMA质控参数列表"""
        args = []
        args.append(f"-miss {self.gemma.miss}")
        args.append(f"-maf {self.gemma.maf}")
        return args

    def get_gemma_output_dir(self) -> str:
        """获取GEMMA输出目录路径"""
        return os.path.join(self.outdir, "output")
