"""
微观共线性分析配置类|Microsynteny Analysis Configuration
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path


@dataclass
class MicrosyntenyConfig:
    """微观共线性分析配置类|Microsynteny Analysis Configuration

    Attributes:
        genome_folder: 基因组文件夹路径|Genome folder path
        gene_list: 目标基因列表文件|Target gene list file
        output_dir: 输出目录|Output directory
        jcvi_path: JCVI Python环境路径（用于共线性分析）|JCVI Python environment path (for synteny analysis)
        threads: 线程数|Number of threads
        extend_genes: 延伸基因数|Number of genes to extend on each side
        cscore: 共线性分数阈值|Synteny score threshold
        step: 运行指定步骤|Run specific step only

    Note:
        步骤说明|Step description:
        - Step 1: GFF/BED预处理|GFF/BED preprocessing
        - Step 2: 共线性分析（使用JCVI）|Synteny analysis (using JCVI)
        - Step 3: 提取共线性blocks|Extract synteny blocks
        - Step 4: 可视化（使用pyCirclize）|Visualization (using pyCirclize)
    """

    # 必需参数|Required parameters
    genome_folder: str
    gene_list: str

    # 可选参数|Optional parameters
    output_dir: str = "./microsynteny_output"
    jcvi_path: str = field(
        default_factory=lambda: get_tool_path(
            'jcvi',
            '~/miniforge3/envs/jcvi_v.1.5.7',
            'JCVI_PATH'
        )
    )
    pycirclize_path: str = field(
        default_factory=lambda: get_tool_path(
            'pycirclize',
            '~/miniforge3/envs/pycirclize_v.1.10.1',
            'PYCIRCCLIZE_PATH'
        )
    )
    threads: int = 12
    extend_genes: int = 30
    cscore: float = 0.99
    step: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing

        展开所有路径并创建输出目录结构|Expand all paths and create output directory structure
        """
        # ⚠️ 关键：展开所有包含~的路径|CRITICAL: Expand all paths with ~
        self.genome_folder = expand_path(self.genome_folder)
        self.gene_list = expand_path(self.gene_list)
        self.output_dir = expand_path(self.output_dir)
        self.jcvi_path = expand_path(self.jcvi_path)
        self.pycirclize_path = expand_path(self.pycirclize_path)

        # 创建输出目录结构|Create output directory structure
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建子目录|Create subdirectories
        self.preprocess_dir = self.output_path / "1_preprocess"
        self.synteny_dir = self.output_path / "2_synteny"
        self.blocks_dir = self.output_path / "3_blocks"
        self.plot_dir = self.output_path / "4_plot"
        self.logs_dir = self.output_path / "logs"

        for dir_path in [self.preprocess_dir, self.synteny_dir,
                         self.blocks_dir, self.plot_dir, self.logs_dir]:
            dir_path.mkdir(exist_ok=True)

        # 日志文件路径|Log file path
        self.log_file = self.logs_dir / "microsynteny.log"

        # 解析物种和基因列表|Parse species and gene list
        self.species_genes = self._parse_gene_list()

    def _find_genome_file(self, species_id: str, extensions: List[str]) -> Optional[str]:
        """查找基因组文件（支持多种后缀）|Find genome file (support multiple extensions)

        Args:
            species_id: 物种ID|Species ID
            extensions: 文件后缀列表|File extension list (e.g., ['.fa', '.fasta'])

        Returns:
            Optional[str]: 找到的文件路径|Found file path, or None if not found
        """
        for ext in extensions:
            file_path = os.path.join(self.genome_folder, f"{species_id}{ext}")
            if os.path.exists(file_path):
                return file_path
        return None

    def _parse_gene_list(self) -> dict:
        """解析基因列表文件|Parse gene list file

        Returns:
            dict: {species_id: [gene_ids]}
        """
        species_genes = {}

        if not os.path.exists(self.gene_list):
            raise ValueError(f"基因列表文件不存在|Gene list file not found: {self.gene_list}")

        with open(self.gene_list, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t' if '\t' in line else None)
                if len(parts) < 2:
                    continue

                species_id = parts[0].strip()
                gene_id = parts[1].strip()

                if species_id not in species_genes:
                    species_genes[species_id] = []
                species_genes[species_id].append(gene_id)

        if not species_genes:
            raise ValueError(f"基因列表文件为空或格式错误|Gene list file is empty or malformed: {self.gene_list}")

        return species_genes

    def validate(self):
        """验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 如果参数验证失败|If validation fails
        """
        errors = []

        # 检查基因组文件夹|Check genome folder
        if not os.path.exists(self.genome_folder):
            errors.append(f"基因组文件夹不存在|Genome folder not found: {self.genome_folder}")
        elif not os.path.isdir(self.genome_folder):
            errors.append(f"基因组路径不是文件夹|Genome path is not a directory: {self.genome_folder}")
        else:
            # 检查是否有数据文件|Check if data files exist
            # 支持多种文件后缀|Support multiple file extensions:
            # 序列文件: .fa, .fasta|Sequence files: .fa, .fasta
            # GFF文件: .gff, .gff3|GFF files: .gff, .gff3
            species_ids = set(self.species_genes.keys())
            missing_files = []

            # 存储找到的文件路径|Store found file paths
            self.genome_files = {}

            for sp_id in species_ids:
                # 查找序列文件|Find sequence file
                fa_extensions = ['.fa', '.fasta']
                fa_file = self._find_genome_file(sp_id, fa_extensions)

                # 查找GFF文件|Find GFF file
                gff_extensions = ['.gff', '.gff3']
                gff_file = self._find_genome_file(sp_id, gff_extensions)

                if not fa_file:
                    missing_files.append(f"{sp_id}.fa/.fasta")
                else:
                    self.genome_files[f"{sp_id}_fa"] = fa_file

                if not gff_file:
                    missing_files.append(f"{sp_id}.gff/.gff3")
                else:
                    self.genome_files[f"{sp_id}_gff"] = gff_file

            if missing_files:
                errors.append(f"缺少基因组文件|Missing genome files: {', '.join(missing_files)}")

        # 检查JCVI环境|Check JCVI environment
        if not os.path.exists(self.jcvi_path):
            errors.append(f"JCVI环境不存在|JCVI environment not found: {self.jcvi_path}")

        jcvi_python = os.path.join(self.jcvi_path, "bin", "python")
        if not os.path.exists(jcvi_python):
            errors.append(f"JCVI Python不存在|JCVI Python not found: {jcvi_python}")

        # 检查LAST工具|Check LAST tool (required by JCVI)
        import shutil
        jcvi_bin = os.path.join(self.jcvi_path, "bin")
        lastdb_in_jcvi = os.path.join(jcvi_bin, "lastdb")

        # 优先检查JCVI环境的bin目录|Check JCVI bin directory first
        if os.path.exists(lastdb_in_jcvi):
            pass  # LAST在JCVI环境中|LAST in JCVI environment
        elif not shutil.which("lastdb"):
            errors.append(
                "未找到LAST工具|LAST tool not found.\n"
                f"在JCVI环境中安装|Install in JCVI environment:\n"
                f"  conda activate jcvi_v.1.5.7\n"
                f"  conda install -c bioconda last\n"
                "或访问|or visit: https://gitlab.com/mcfrith/last"
            )

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查cscore范围|Check cscore range
        if not (0 <= self.cscore <= 1):
            errors.append(f"C-score必须在0-1之间|C-score must be between 0 and 1: {self.cscore}")

        # 检查step参数|Check step parameter
        valid_steps = ['1', '2', '3', '4', None]
        if self.step not in valid_steps:
            errors.append(f"无效的步骤参数|Invalid step parameter: {self.step}")

        # 检查pyCirclize环境（如果运行Step 4或全流程）
        if self.step in ['4', None]:
            if not os.path.exists(self.pycirclize_path):
                errors.append(f"pyCirclize环境不存在|pyCirclize environment not found: {self.pycirclize_path}")

            pycirclize_python = os.path.join(self.pycirclize_path, "bin", "python")
            if not os.path.exists(pycirclize_python):
                errors.append(f"pyCirclize Python不存在|pyCirclize Python not found: {pycirclize_python}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_species_list(self) -> List[str]:
        """获取物种列表|Get species list

        Returns:
            List[str]: 物种ID列表|List of species IDs
        """
        return list(self.species_genes.keys())

    def get_target_genes_for_species_pair(self, sp1: str, sp2: str) -> set:
        """获取物种对的目标基因集合|Get target genes set for species pair

        Args:
            sp1: 物种1|Species 1
            sp2: 物种2|Species 2

        Returns:
            set: 两个物种的目标基因集合|Set of target genes from both species
        """
        target_genes = set()

        # 添加物种1的目标基因
        if sp1 in self.species_genes:
            target_genes.update(self.species_genes[sp1])

        # 添加物种2的目标基因
        if sp2 in self.species_genes:
            target_genes.update(self.species_genes[sp2])

        return target_genes

    def get_done_file(self, step: str) -> Path:
        """获取步骤完成标记文件路径|Get step done marker file path

        Args:
            step: 步骤编号|Step number (1/2/3/4)

        Returns:
            Path: 完成标记文件路径|Done marker file path
        """
        return self.output_path / f".step_{step}.done"

    def is_step_done(self, step: str) -> bool:
        """检查步骤是否已完成|Check if step is completed

        Args:
            step: 步骤编号|Step number (1/2/3/4)

        Returns:
            bool: 是否已完成|Whether completed
        """
        return self.get_done_file(step).exists()

    def mark_step_done(self, step: str):
        """标记步骤已完成|Mark step as completed

        Args:
            step: 步骤编号|Step number (1/2/3/4)
        """
        self.get_done_file(step).touch()
