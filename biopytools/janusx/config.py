"""
JanusX配置模块|JanusX Configuration Module

包含JanusX各模块的配置类|Contains configuration classes for JanusX modules
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
from ..common.paths import get_tool_path


@dataclass
class JanusXGWASConfig:
    """JanusX GWAS配置类|JanusX GWAS Configuration Class"""

    # 必需参数|Required parameters
    genotype: str  # VCF或PLINK文件路径|VCF or PLINK file path
    pheno: str     # 表型文件路径|Phenotype file path

    # 可选参数|Optional parameters
    genotype_type: str = "vcf"  # 基因型类型: vcf 或 bfile|Genotype type: vcf or bfile
    models: List[str] = field(default_factory=list)  # GWAS模型: lm, lmm, fastlmm, farmcpu
    maf: float = 0.05  # 最小等位基因频率阈值|Minor allele frequency threshold
    geno: float = 0.0  # 缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)
    grm: str = "1"  # 亲缘关系矩阵方法|GRM method: 1, 2, or path to precomputed GRM
    qcov: int = 0  # PCA数量或Q矩阵文件路径|Number of PCs or path to Q matrix file (0=不使用Q矩阵|no Q matrix)
    cov: Optional[str] = None  # 协变量文件|Covariate file path
    ncol: Optional[List[int]] = None  # 表型列索引(零基)|Phenotype column indices (zero-based)
    plot: bool = False  # 是否生成图表|Whether to generate plots
    chunksize: int = 100000  # SNP分块大小|SNP chunk size
    mmap_limit: Optional[int] = None  # 内存映射限制|Memory map limit
    threads: int = 12  # 线程数|Number of threads
    output_dir: str = "./janusx_gwas_output"  # 输出目录|Output directory
    prefix: Optional[str] = None  # 输出文件前缀|Output file prefix

    # JanusX路径配置|JanusX path configuration
    janusx_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认JanusX路径|Set default JanusX path
        if self.janusx_path is None:
            self.janusx_path = get_tool_path('janusx', '~/software/JanusX/JanusX.bin', 'JANUSX_PATH')
        else:
            self.janusx_path = get_tool_path('janusx', self.janusx_path, 'JANUSX_PATH')

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        self.genotype = os.path.abspath(self.genotype) if self.genotype else None
        self.pheno = os.path.abspath(self.pheno) if self.pheno else None
        if self.cov:
            self.cov = os.path.abspath(self.cov)

        # 设置默认模型|Set default models
        if not self.models:
            self.models = ['lmm']

        # 设置默认前缀|Set default prefix
        if self.prefix is None:
            genotype_basename = Path(self.genotype).stem
            self.prefix = genotype_basename

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.genotype:
            errors.append("基因型文件未指定|Genotype file not specified")

        if self.genotype and not os.path.exists(self.genotype):
            errors.append(f"基因型文件不存在|Genotype file not found: {self.genotype}")

        if not self.pheno:
            errors.append("表型文件未指定|Phenotype file not specified")

        if self.pheno and not os.path.exists(self.pheno):
            errors.append(f"表型文件不存在|Phenotype file not found: {self.pheno}")

        # 验证可选文件参数|Validate optional file parameters
        if self.cov and not os.path.exists(self.cov):
            errors.append(f"协变量文件不存在|Covariate file not found: {self.cov}")

        # 验证JanusX路径|Validate JanusX path
        if self.janusx_path and not os.path.exists(self.janusx_path):
            errors.append(f"JanusX路径不存在|JanusX path not found: {self.janusx_path}")

        # 验证基因型类型|Validate genotype type
        if self.genotype_type not in ['vcf', 'bfile']:
            errors.append(f"基因型类型必须是vcf或bfile|Genotype type must be vcf or bfile: {self.genotype_type}")

        # 验证模型|Validate models
        valid_models = ['lm', 'lmm', 'fastlmm', 'farmcpu']
        for model in self.models:
            if model not in valid_models:
                errors.append(f"无效的模型|Invalid model: {model}")

        # 验证参数范围|Validate parameter ranges
        if self.maf < 0 or self.maf > 1:
            errors.append(f"MAF参数必须在0-1之间|MAF must be between 0 and 1: {self.maf}")

        if self.geno < 0 or self.geno > 1:
            errors.append(f"geno参数必须在0-1之间|geno must be between 0 and 1: {self.geno}")

        if self.chunksize <= 0:
            errors.append(f"chunksize参数必须大于0|chunksize must be positive: {self.chunksize}")

        if self.threads < -1:
            errors.append(f"线程数必须>= -1|Thread count must be >= -1: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class JanusXGSConfig:
    """JanusX基因组选择配置类|JanusX Genomic Selection Configuration Class"""

    # 必需参数|Required parameters
    genotype: str  # VCF或PLINK文件路径|VCF or PLINK file path
    pheno: str     # 表型文件路径|Phenotype file path

    # 可选参数|Optional parameters
    genotype_type: str = "vcf"  # 基因型类型: vcf 或 bfile|Genotype type: vcf or bfile
    models: List[str] = field(default_factory=list)  # GS模型: GBLUP, rrBLUP, BayesA, BayesB, BayesCpi
    maf: float = 0.05  # 最小等位基因频率阈值|Minor allele frequency threshold
    geno: float = 0.0  # 缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)
    pcd: bool = False  # 启用PCA降维|Enable PCA-based dimensionality reduction
    ncol: Optional[int] = None  # 表型列索引(零基)|Phenotype column index (zero-based)
    cv: Optional[int] = None  # 交叉验证折数|K-fold cross-validation
    plot: bool = False  # 是否生成图表|Whether to generate plots
    output_dir: str = "./janusx_gs_output"  # 输出目录|Output directory
    prefix: Optional[str] = None  # 输出文件前缀|Output file prefix

    # JanusX路径配置|JanusX path configuration
    janusx_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认JanusX路径|Set default JanusX path
        if self.janusx_path is None:
            self.janusx_path = get_tool_path('janusx', '~/software/JanusX/JanusX.bin', 'JANUSX_PATH')
        else:
            self.janusx_path = get_tool_path('janusx', self.janusx_path, 'JANUSX_PATH')

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        self.genotype = os.path.abspath(self.genotype) if self.genotype else None
        self.pheno = os.path.abspath(self.pheno) if self.pheno else None

        # 设置默认模型|Set default models
        if not self.models:
            self.models = ['GBLUP']

        # 设置默认前缀|Set default prefix
        if self.prefix is None:
            genotype_basename = Path(self.genotype).stem
            self.prefix = genotype_basename

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.genotype:
            errors.append("基因型文件未指定|Genotype file not specified")

        if self.genotype and not os.path.exists(self.genotype):
            errors.append(f"基因型文件不存在|Genotype file not found: {self.genotype}")

        if not self.pheno:
            errors.append("表型文件未指定|Phenotype file not specified")

        if self.pheno and not os.path.exists(self.pheno):
            errors.append(f"表型文件不存在|Phenotype file not found: {self.pheno}")

        # 验证JanusX路径|Validate JanusX path
        if self.janusx_path and not os.path.exists(self.janusx_path):
            errors.append(f"JanusX路径不存在|JanusX path not found: {self.janusx_path}")

        # 验证基因型类型|Validate genotype type
        if self.genotype_type not in ['vcf', 'bfile']:
            errors.append(f"基因型类型必须是vcf或bfile|Genotype type must be vcf or bfile: {self.genotype_type}")

        # 验证模型|Validate models
        valid_models = ['GBLUP', 'rrBLUP', 'BayesA', 'BayesB', 'BayesCpi']
        for model in self.models:
            if model not in valid_models:
                errors.append(f"无效的模型|Invalid model: {model}")

        # 验证参数范围|Validate parameter ranges
        if self.maf < 0 or self.maf > 1:
            errors.append(f"MAF参数必须在0-1之间|MAF must be between 0 and 1: {self.maf}")

        if self.geno < 0 or self.geno > 1:
            errors.append(f"geno参数必须在0-1之间|geno must be between 0 and 1: {self.geno}")

        if self.ncol is not None and self.ncol < 0:
            errors.append(f"ncol参数必须>= 0|ncol must be >= 0: {self.ncol}")

        if self.cv is not None and self.cv <= 0:
            errors.append(f"cv参数必须大于0|cv must be positive: {self.cv}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class JanusXPCAConfig:
    """JanusX PCA配置类|JanusX PCA Configuration Class"""

    # 必需参数|Required parameters
    genotype: str  # VCF或PLINK文件路径|VCF or PLINK file path

    # 可选参数|Optional parameters
    genotype_type: str = "vcf"  # 基因型类型: vcf 或 bfile|Genotype type: vcf or bfile
    grm: Optional[str] = None  # 预计算的GRM前缀|Precomputed GRM prefix
    pcfile: Optional[str] = None  # 预计算的PCA文件|Precomputed PCA file
    dim: int = 3  # 输出主成分数量|Number of PCs to output
    plot: bool = False  # 生成2D散点图|Generate 2D scatter plot
    plot3d: bool = False  # 生成3D旋转GIF|Generate 3D rotating GIF
    group: Optional[str] = None  # 分组文件路径|Group file path
    color: int = 1  # 调色板索引(0-6)|Color palette index (0-6)
    output_dir: str = "./janusx_pca_output"  # 输出目录|Output directory
    prefix: Optional[str] = None  # 输出文件前缀|Output file prefix

    # JanusX路径配置|JanusX path configuration
    janusx_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认JanusX路径|Set default JanusX path
        if self.janusx_path is None:
            self.janusx_path = get_tool_path('janusx', '~/software/JanusX/JanusX.bin', 'JANUSX_PATH')
        else:
            self.janusx_path = get_tool_path('janusx', self.janusx_path, 'JANUSX_PATH')

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        if self.genotype:
            self.genotype = os.path.abspath(self.genotype)
        if self.grm:
            self.grm = os.path.abspath(self.grm)
        if self.pcfile:
            self.pcfile = os.path.abspath(self.pcfile)
        if self.group:
            self.group = os.path.abspath(self.group)

        # 设置默认前缀|Set default prefix
        if self.prefix is None:
            genotype_basename = Path(self.genotype).stem
            self.prefix = genotype_basename

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.genotype:
            errors.append("基因型文件未指定|Genotype file not specified")

        if self.genotype and not os.path.exists(self.genotype):
            errors.append(f"基因型文件不存在|Genotype file not found: {self.genotype}")

        # 验证可选文件参数|Validate optional file parameters
        if self.grm and not os.path.exists(self.grm):
            errors.append(f"GRM文件不存在|GRM file not found: {self.grm}")

        if self.pcfile and not os.path.exists(self.pcfile):
            errors.append(f"PCA文件不存在|PCA file not found: {self.pcfile}")

        if self.group and not os.path.exists(self.group):
            errors.append(f"分组文件不存在|Group file not found: {self.group}")

        # 验证JanusX路径|Validate JanusX path
        if self.janusx_path and not os.path.exists(self.janusx_path):
            errors.append(f"JanusX路径不存在|JanusX path not found: {self.janusx_path}")

        # 验证基因型类型|Validate genotype type
        if self.genotype_type not in ['vcf', 'bfile']:
            errors.append(f"基因型类型必须是vcf或bfile|Genotype type must be vcf or bfile: {self.genotype_type}")

        # 验证参数范围|Validate parameter ranges
        if self.dim <= 0:
            errors.append(f"dim参数必须大于0|dim must be positive: {self.dim}")

        if self.color < 0 or self.color > 6:
            errors.append(f"color参数必须在0-6之间|color must be between 0 and 6: {self.color}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class JanusXPostGWASConfig:
    """JanusX PostGWAS配置类|JanusX PostGWAS Configuration Class"""

    # 必需参数|Required parameters
    files: List[str]  # GWAS结果文件列表|List of GWAS result files

    # 可选参数|Optional parameters
    chr_col: str = "#CHROM"  # 染色体列名|Chromosome column name
    pos_col: str = "POS"  # 位置列名|Position column name
    pvalue_col: str = "p"  # P值列名|P-value column name
    threshold: Optional[float] = None  # 显著性阈值|Significance threshold
    noplot: bool = False  # 禁用绘图|Disable plotting
    color: int = 0  # 颜色方案索引(0-6, -1为auto)|Color scheme index (0-6, -1 for auto)
    highlight: Optional[str] = None  # 高亮区域BED文件|Highlight regions BED file
    format: str = "png"  # 输出格式: pdf, png, svg, tif|Output format
    anno: Optional[str] = None  # 注释文件路径|Annotation file path
    anno_broaden: Optional[int] = None  # 注释窗口|Annotation window (kb)
    desc_item: str = "description"  # GFF描述键|GFF description key
    output_dir: str = "./janusx_postgwas_output"  # 输出目录|Output directory
    prefix: str = "JanusX"  # 输出文件前缀|Output file prefix
    threads: int = -1  # 线程数|Number of threads (-1 for all cores)

    # JanusX路径配置|JanusX path configuration
    janusx_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认JanusX路径|Set default JanusX path
        if self.janusx_path is None:
            self.janusx_path = get_tool_path('janusx', '~/software/JanusX/JanusX.bin', 'JANUSX_PATH')
        else:
            self.janusx_path = get_tool_path('janusx', self.janusx_path, 'JANUSX_PATH')

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        self.files = [os.path.abspath(f) for f in self.files]
        if self.highlight:
            self.highlight = os.path.abspath(self.highlight)
        if self.anno:
            self.anno = os.path.abspath(self.anno)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.files:
            errors.append("未指定GWAS结果文件|GWAS result files not specified")

        for file in self.files:
            if not os.path.exists(file):
                errors.append(f"GWAS结果文件不存在|GWAS result file not found: {file}")

        # 验证可选文件参数|Validate optional file parameters
        if self.highlight and not os.path.exists(self.highlight):
            errors.append(f"高亮文件不存在|Highlight file not found: {self.highlight}")

        if self.anno and not os.path.exists(self.anno):
            errors.append(f"注释文件不存在|Annotation file not found: {self.anno}")

        # 验证JanusX路径|Validate JanusX path
        if self.janusx_path and not os.path.exists(self.janusx_path):
            errors.append(f"JanusX路径不存在|JanusX path not found: {self.janusx_path}")

        # 验证输出格式|Validate output format
        valid_formats = ['pdf', 'png', 'svg', 'tif']
        if self.format not in valid_formats:
            errors.append(f"输出格式必须是{valid_formats}|Output format must be one of {valid_formats}: {self.format}")

        # 验证参数范围|Validate parameter ranges
        if self.threshold is not None and (self.threshold <= 0 or self.threshold >= 1):
            errors.append(f"threshold参数必须在0-1之间|threshold must be between 0 and 1: {self.threshold}")

        if self.color < -1 or self.color > 6:
            errors.append(f"color参数必须在-1到6之间|color must be between -1 and 6: {self.color}")

        if self.anno_broaden is not None and self.anno_broaden < 0:
            errors.append(f"anno_broaden参数必须>= 0|anno_broaden must be >= 0: {self.anno_broaden}")

        if self.threads < -1:
            errors.append(f"线程数必须>= -1|Thread count must be >= -1: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
