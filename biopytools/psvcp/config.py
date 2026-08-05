"""PSVCP 配置管理模块|PSVCP Configuration Module

参数暴露照搬原始 1Genome_construct_Pangenome.py:仅 genome_dir / genome_list / -t。
assemblytics、nucmer、ins 阈值为原始硬编码值(内部常量,不上 CLI)。
|Param exposure mirrors the original script: only genome_dir / genome_list / -t.
assemblytics / nucmer / insertion thresholds are the original hardcoded values
(internal constants, not exposed on CLI).
"""

import os
from dataclasses import dataclass, field

from ..common.paths import expand_path, get_tool_path


@dataclass
class PSVCPConfig:
    """PSVCP 泛基因组构建配置|PSVCP pangenome construction configuration"""

    # ===== 必需输入|Required input =====
    genome_dir: str           # 含 {name}.fa + {name}.gff/.gff3 的目录|dir with {name}.fa + {name}.gff/.gff3
    genome_list: str          # 文本,行1=ref,其余=query(顺序即并入顺序)|line1=ref, rest=queries

    # ===== biopytools 新增(替代原始 CWD)|biopytools additions (replace original CWD) =====
    output_dir: str = "~/psvcp_out"
    threads: int = 12         # = 原始 -t (nucmer -t)|= original -t (nucmer -t)
    force: bool = False       # 忽略断点续传重跑|ignore checkpoint, rerun

    # ===== 原始硬编码值(内部常量,不暴露 CLI)|original hardcoded values =====
    # assemblytics(原始 1000 50 10000000)
    ASSEMBLYTICS_UNIQUE_LENGTH: int = 1000   # 非 2.0.1 默认 10000|NOT 2.0.1 default 10000
    ASSEMBLYTICS_MIN_SIZE: int = 50
    ASSEMBLYTICS_MAX_SIZE: int = 10000000
    # nucmer(原始 --maxgap 500 --mincluster 1000 --diagdiff 20)
    NUCMER_MAXGAP: int = 500
    NUCMER_MINCLUSTER: int = 1000
    NUCMER_DIAGDIFF: int = 20
    # 插入阈值(原始 4ins_more_50.py 的 >50)|insertion threshold (original >50)
    MIN_INSERTION_SIZE: int = 50
    # 原始工作目录名,照搬|original working-dir name, kept as-is
    PAN_DIR_NAME: str = "pan_dir_result"

    # ===== 工具路径(单一 env psvcp_v.1.0.1)|tool paths (single env) =====
    nucmer_path: str = field(
        default_factory=lambda: get_tool_path(
            'nucmer', '~/miniforge3/envs/psvcp_v.1.0.1/bin/nucmer', 'NUCMER_PATH')
    )
    assemblytics_path: str = field(
        default_factory=lambda: get_tool_path(
            'assemblytics', '~/miniforge3/envs/psvcp_v.1.0.1/bin/assemblytics', 'ASSEMBLYTICS_PATH')
    )
    bedtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'bedtools', '~/miniforge3/envs/psvcp_v.1.0.1/bin/bedtools', 'BEDTOOLS_PATH')
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'samtools', '~/miniforge3/envs/psvcp_v.1.0.1/bin/samtools', 'SAMTOOLS_PATH')
    )
    rscript_path: str = field(
        default_factory=lambda: get_tool_path(
            'rscript', '~/miniforge3/envs/psvcp_v.1.0.1/bin/Rscript', 'RSCRIPT_PATH')
    )
    python_path: str = field(
        default_factory=lambda: get_tool_path(
            'python3', '~/miniforge3/envs/psvcp_v.1.0.1/bin/python3', 'PSVCP_PYTHON_PATH')
    )

    def __post_init__(self):
        """展开所有 ~ 路径|Expand all ~ paths (§11.3.1 关键)"""
        self.genome_dir = os.path.normpath(os.path.abspath(expand_path(self.genome_dir)))
        self.genome_list = os.path.normpath(os.path.abspath(expand_path(self.genome_list)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        for t in ('nucmer_path', 'assemblytics_path', 'bedtools_path',
                  'samtools_path', 'rscript_path', 'python_path'):
            setattr(self, t, expand_path(getattr(self, t)))

    @property
    def pan_dir(self) -> str:
        """泛基因组工作目录|pangenome working dir (= output_dir/pan_dir_result)"""
        return os.path.join(self.output_dir, self.PAN_DIR_NAME)

    def read_genome_list(self):
        """读取 genome_list,返回非空行列表|Read genome_list, return non-empty stripped lines"""
        with open(self.genome_list, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip()]

    @staticmethod
    def find_genome_gff(genome_dir: str, name: str):
        """在 genome_dir 下按 .gff > .gff3 查找基因组注释,返回首个命中路径,均无则 None
        |Find genome annotation under genome_dir (.gff preferred over .gff3);
        return first existing path, or None if neither present"""
        stem = name[:-3] if name.endswith('.fa') else name
        for ext in ('.gff', '.gff3'):
            candidate = os.path.join(genome_dir, stem + ext)
            if os.path.isfile(candidate):
                return candidate
        return None

    def validate(self):
        """验证配置|Validate configuration (§六: 一次性收集错误)"""
        errors = []

        # genome_dir|genome directory
        if not os.path.isdir(self.genome_dir):
            errors.append(f"genome 目录不存在|genome directory not found: {self.genome_dir}")

        # genome_list|genome list
        if not os.path.isfile(self.genome_list):
            errors.append(f"genome_list 不存在|genome_list not found: {self.genome_list}")
        else:
            names = self.read_genome_list()
            if len(names) < 2:
                errors.append(
                    f"genome_list 至少需要 2 个基因组(1 ref + ≥1 query)|"
                    f"genome_list needs >=2 genomes (1 ref + >=1 query), 当前|got: {len(names)}"
                )
            # 逐行查 .fa + .gff/.gff3|check each name has .fa + .gff/.gff3
            for name in names:
                fa = os.path.join(self.genome_dir, name)
                if not os.path.isfile(fa):
                    errors.append(f"基因组 fasta 不存在|genome fasta not found: {fa}")
                gff = self.find_genome_gff(self.genome_dir, name)
                if gff is None:
                    stem = name[:-3] if name.endswith('.fa') else name
                    errors.append(
                        f"基因组 gff 不存在|genome gff not found: "
                        f"{os.path.join(self.genome_dir, stem)}.gff/.gff3"
                    )

        # threads|threads
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))
        return True
