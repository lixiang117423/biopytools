"""
GenomeSyn2配置管理模块|GenomeSyn2 Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class GenomeSyn2Config:
    """GenomeSyn2配置类|GenomeSyn2 Configuration Class"""

    # ========== 比对模式配置 | Alignment Mode Configuration ==========

    # 比对软件类型|Alignment software type
    align: Optional[str] = None  # mummer, minimap2, blastp, mmseqs, diamond

    # 基因组文件目录|Genome files directory
    genome_dir: Optional[str] = None

    # 基因注释文件目录|Gene annotation files directory
    gene_dir: Optional[str] = None

    # 输出目录|Output directory
    outdir: Optional[str] = None

    # 线程数|Number of threads
    threads: int = 12

    # ========== VCF模式配置 | VCF Mode Configuration ==========

    # VCF文件路径|VCF file path
    vcf: Optional[str] = None

    # Bin大小|Bin size for SNP analysis
    bin: int = 50000

    # SNP一致性文件|SNP identity file
    identity: Optional[str] = None

    # SNP密度文件|SNP density file
    density: Optional[str] = None

    # ========== 绘图模式配置 | Plotting Mode Configuration ==========

    # 配置文件路径|Configuration file path
    conf: Optional[str] = None

    # 注释配置|Annotation configuration flag
    anno: bool = False

    # ========== 文件生成模式配置 | File Generation Mode ==========

    # 文件类型|File type for generation
    type: Optional[str] = None  # fa, prot, anno

    # 路径|Path for file generation
    path: Optional[str] = None

    # 输出文件名|Output file name
    out: Optional[str] = None

    # ========== 软件路径配置 | Software Path Configuration ==========

    # GenomeSyn2 Perl脚本路径|GenomeSyn2 Perl script path
    genomesyn2_pl: str = '~/miniforge3/envs/genomesyn2/bin/GenomeSyn2.pl'

    # Perl路径|Perl path
    perl_path: str = '~/miniforge3/envs/genomesyn2/bin/perl'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.genomesyn2_pl = expand_path(self.genomesyn2_pl)
        self.perl_path = expand_path(self.perl_path)

        # 如果指定了输出目录，确保其存在
        if self.outdir:
            self.output_path = Path(self.outdir)
            self.output_path.mkdir(parents=True, exist_ok=True)
        else:
            self.output_path = None

        # 标准化文件路径|Normalize file paths
        if self.genome_dir:
            self.genome_dir = os.path.normpath(os.path.abspath(self.genome_dir))
        if self.gene_dir:
            self.gene_dir = os.path.normpath(os.path.abspath(self.gene_dir))
        if self.vcf:
            self.vcf = os.path.normpath(os.path.abspath(self.vcf))
        if self.identity:
            self.identity = os.path.normpath(os.path.abspath(self.identity))
        if self.density:
            self.density = os.path.normpath(os.path.abspath(self.density))
        if self.conf:
            self.conf = os.path.normpath(os.path.abspath(self.conf))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查GenomeSyn2 Perl脚本是否存在
        if not os.path.exists(self.genomesyn2_pl):
            errors.append(
                f"GenomeSyn2 Perl脚本不存在|GenomeSyn2 Perl script not found: "
                f"{self.genomesyn2_pl}"
            )

        # 检查Perl是否存在
        if not os.path.exists(self.perl_path):
            errors.append(
                f"Perl不存在|Perl not found: {self.perl_path}"
            )

        # 验证比对模式参数
        if self.align:
            if not self.genome_dir:
                errors.append(
                    "比对模式需要指定--genome参数|"
                    "Alignment mode requires --genome parameter"
                )
            else:
                if not os.path.exists(self.genome_dir):
                    errors.append(
                        f"基因组目录不存在|Genome directory not found: {self.genome_dir}"
                    )

            # 蛋白质比对需要基因注释目录
            if self.align in ['blastp', 'mmseqs', 'diamond']:
                if not self.gene_dir:
                    errors.append(
                        f"{self.align}模式需要指定--gene参数|"
                        f"{self.align} mode requires --gene parameter"
                    )
                elif not os.path.exists(self.gene_dir):
                    errors.append(
                        f"基因注释目录不存在|Gene directory not found: {self.gene_dir}"
                    )

            # 验证比对软件类型
            valid_aligners = ['mummer', 'minimap2', 'blastp', 'mmseqs', 'diamond']
            if self.align not in valid_aligners:
                errors.append(
                    f"无效的比对软件|Invalid aligner: {self.align} "
                    f"(必须为|must be one of {', '.join(valid_aligners)})"
                )

        # 验证VCF模式参数
        if self.vcf:
            if not os.path.exists(self.vcf):
                errors.append(
                    f"VCF文件不存在|VCF file not found: {self.vcf}"
                )

        # 验证SNP文件模式参数
        if self.identity and self.density:
            if not os.path.exists(self.identity):
                errors.append(
                    f"SNP一致性文件不存在|SNP identity file not found: {self.identity}"
                )
            if not os.path.exists(self.density):
                errors.append(
                    f"SNP密度文件不存在|SNP density file not found: {self.density}"
                )

        # 验证文件生成模式参数
        if self.type:
            valid_types = ['fa', 'prot', 'anno']
            if self.type not in valid_types:
                errors.append(
                    f"无效的文件类型|Invalid file type: {self.type} "
                    f"(必须为|must be one of {', '.join(valid_types)})"
                )
            if not self.path:
                errors.append(
                    "文件生成模式需要指定--path参数|"
                    "File generation mode requires --path parameter"
                )
            elif not os.path.exists(self.path):
                errors.append(
                    f"指定路径不存在|Specified path not found: {self.path}"
                )

        # 验证配置文件模式
        if self.conf and not os.path.exists(self.conf):
            errors.append(
                f"配置文件不存在|Configuration file not found: {self.conf}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
