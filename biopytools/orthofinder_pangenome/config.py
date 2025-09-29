# """
# OrthoFinder泛基因组分析配置管理模块 | OrthoFinder Pangenome Analysis Configuration Management Module
# """

# import os
# from dataclasses import dataclass
# from pathlib import Path
# from typing import Optional

# @dataclass
# class PangenomeConfig:
#     """泛基因组分析配置类 | Pangenome Analysis Configuration Class"""
    
#     # 输入输出 | Input/Output
#     input_dir: str
#     output_dir: str = './orthofinder_pangenome_output'
#     project_name: Optional[str] = None
    
#     # OrthoFinder基础参数 | OrthoFinder basic parameters
#     threads: int = 88
#     search_program: str = 'blast'  # blast, diamond, diamond_ultra_sens, mmseqs
#     mcl_inflation: float = 1.2
#     sequence_type: str = 'protein'  # protein, dna
    
#     # 泛基因组分类参数 | Pangenome classification parameters
#     softcore_missing_threshold: int = 2        # Softcore基因缺失阈值 | Softcore missing threshold
#     dispensable_missing_threshold: int = 2     # Dispensable基因缺失阈值 | Dispensable missing threshold
    
#     # 稀释分析参数 | Rarefaction analysis parameters
#     enable_rarefaction: bool = True            # 启用稀释曲线分析 | Enable rarefaction analysis
#     rarefaction_iterations: int = 100          # 每个样本量的迭代次数 | Iterations per sample size

#     # 单拷贝基因分析参数 | Single copy gene analysis parameters
#     enable_single_copy_analysis: bool = True    # 启用单拷贝基因分析 | Enable single copy gene analysis
#     extract_sequences: bool = True              # 提取序列 | Extract sequences
#     single_copy_output_format: str = 'both'     # 输出格式: 'by_orthogroup', 'by_genome', 'both'
    
#     # 分析模式 | Analysis mode
#     basic_analysis_only: bool = True
#     generate_trees: bool = False
#     msa_program: str = 'mafft'  # mafft, muscle
#     tree_program: str = 'fasttree'  # fasttree, raxml, iqtree
    
#     # 断点续跑参数 | Resume parameters
#     resume_from_existing: bool = True          # 默认启用断点续跑 | Default enable resume
#     skip_orthofinder: bool = False             # 跳过OrthoFinder步骤 | Skip OrthoFinder step
#     force_overwrite: bool = False              # 强制覆盖已有结果 | Force overwrite existing results
    
#     # 可视化参数 | Visualization parameters
#     generate_plots: bool = True                # 生成图表 | Generate plots
#     plot_format: str = 'png'
#     figure_dpi: int = 300
    
#     # 工具路径 | Tool paths
#     orthofinder_path: str = 'orthofinder'
    
#     # 内部属性 | Internal attributes
#     base_name: str = 'pangenome_analysis'
    
#     def __post_init__(self):
#         """初始化后处理 | Post-initialization processing"""
#         self.output_path = Path(self.output_dir)
#         self.output_path.mkdir(parents=True, exist_ok=True)
        
#         # 标准化路径 | Normalize paths
#         self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
#         self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
#         # 设置项目名称 | Set project name
#         if not self.project_name:
#             self.project_name = f"pangenome_{Path(self.input_dir).name}"
    
#     def validate(self):
#         """验证配置参数 | Validate configuration parameters"""
#         errors = []
        
#         # 检查输入目录 | Check input directory
#         if not os.path.exists(self.input_dir):
#             errors.append(f"输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
#         # 检查输入文件 | Check input files
#         input_path = Path(self.input_dir)
#         fasta_files = list(input_path.glob('*.fa*'))
#         if not fasta_files:
#             errors.append(f"输入目录中未找到FASTA文件 | No FASTA files found in input directory: {self.input_dir}")
        
#         # 检查参数范围 | Check parameter ranges
#         if self.threads <= 0:
#             errors.append(f"线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
#         if not 0.1 <= self.mcl_inflation <= 10.0:
#             errors.append(f"MCL inflation参数超出合理范围 | MCL inflation parameter out of range: {self.mcl_inflation}")
        
#         if self.softcore_missing_threshold < 0:
#             errors.append(f"Softcore缺失阈值必须为非负整数 | Softcore missing threshold must be non-negative: {self.softcore_missing_threshold}")
        
#         if self.dispensable_missing_threshold < 0:
#             errors.append(f"Dispensable缺失阈值必须为非负整数 | Dispensable missing threshold must be non-negative: {self.dispensable_missing_threshold}")
        
#         if self.rarefaction_iterations <= 0:
#             errors.append(f"稀释分析迭代次数必须为正整数 | Rarefaction iterations must be positive: {self.rarefaction_iterations}")
        
#         # 检查搜索程序选项 | Check search program options
#         valid_search_programs = ['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']
#         if self.search_program not in valid_search_programs:
#             errors.append(f"不支持的搜索程序 | Unsupported search program: {self.search_program}")
        
#         if errors:
#             raise ValueError("\n".join(errors))
        
#         return True

"""
⚙️ OrthoFinder泛基因组分析配置管理模块 | OrthoFinder Pangenome Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PangenomeConfig:
    """⚙️ 泛基因组分析配置类 | Pangenome Analysis Configuration Class"""
    
    # 📁 输入输出 | Input/Output
    input_dir: str
    output_dir: str = './orthofinder_pangenome_output'
    project_name: Optional[str] = None
    
    # 🔧 OrthoFinder基础参数 | OrthoFinder basic parameters
    threads: int = 88
    search_program: str = 'blast'  # blast, diamond, diamond_ultra_sens, mmseqs
    mcl_inflation: float = 1.2
    sequence_type: str = 'protein'  # protein, dna
    
    # 🎯 泛基因组分类参数 | Pangenome classification parameters
    softcore_missing_threshold: int = 2        # 🟠 Softcore基因缺失阈值 | Softcore missing threshold
    dispensable_missing_threshold: int = 2     # 🟡 Dispensable基因缺失阈值 | Dispensable missing threshold
    
    # 📈 稀释分析参数 | Rarefaction analysis parameters
    enable_rarefaction: bool = True            # 📈 启用稀释曲线分析 | Enable rarefaction analysis
    rarefaction_iterations: int = 100          # 🔢 每个样本量的迭代次数 | Iterations per sample size

    # 🔵 单拷贝基因分析参数 | Single copy gene analysis parameters
    enable_single_copy_analysis: bool = True    # 🔵 启用单拷贝基因分析 | Enable single copy gene analysis
    extract_sequences: bool = True              # 🧬 提取序列 | Extract sequences
    single_copy_output_format: str = 'both'     # 📄 输出格式: 'by_orthogroup', 'by_genome', 'both'
    
    # 🔬 分析模式 | Analysis mode
    basic_analysis_only: bool = True
    generate_trees: bool = False
    msa_program: str = 'mafft'  # mafft, muscle
    tree_program: str = 'fasttree'  # fasttree, raxml, iqtree
    
    # 🔄 断点续跑参数 | Resume parameters
    resume_from_existing: bool = True          # 🔄 默认启用断点续跑 | Default enable resume
    skip_orthofinder: bool = False             # ⏭️ 跳过OrthoFinder步骤 | Skip OrthoFinder step
    force_overwrite: bool = False              # 🔄 强制覆盖已有结果 | Force overwrite existing results
    
    # 📊 可视化参数 | Visualization parameters
    generate_plots: bool = True                # 📊 生成图表 | Generate plots
    plot_format: str = 'png'
    figure_dpi: int = 300
    
    # 🔧 工具路径 | Tool paths
    orthofinder_path: str = 'orthofinder'
    
    # 🏷️ 内部属性 | Internal attributes
    base_name: str = 'pangenome_analysis'
    
    def __post_init__(self):
        """🚀 初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 📁 标准化路径 | Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 🏷️ 设置项目名称 | Set project name
        if not self.project_name:
            self.project_name = f"pangenome_{Path(self.input_dir).name}"
    
    def validate(self):
        """✅ 验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 🔍 检查输入目录 | Check input directory
        if not os.path.exists(self.input_dir):
            errors.append(f"❌ 输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
        # 📄 检查输入文件 | Check input files
        input_path = Path(self.input_dir)
        fasta_files = list(input_path.glob('*.fa*'))
        if not fasta_files:
            errors.append(f"❌ 输入目录中未找到FASTA文件 | No FASTA files found in input directory: {self.input_dir}")
        
        # 🔢 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if not 0.1 <= self.mcl_inflation <= 10.0:
            errors.append(f"❌ MCL inflation参数超出合理范围 | MCL inflation parameter out of range: {self.mcl_inflation}")
        
        if self.softcore_missing_threshold < 0:
            errors.append(f"❌ Softcore缺失阈值必须为非负整数 | Softcore missing threshold must be non-negative: {self.softcore_missing_threshold}")
        
        if self.dispensable_missing_threshold < 0:
            errors.append(f"❌ Dispensable缺失阈值必须为非负整数 | Dispensable missing threshold must be non-negative: {self.dispensable_missing_threshold}")
        
        if self.rarefaction_iterations <= 0:
            errors.append(f"❌ 稀释分析迭代次数必须为正整数 | Rarefaction iterations must be positive: {self.rarefaction_iterations}")
        
        # 🔍 检查搜索程序选项 | Check search program options
        valid_search_programs = ['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']
        if self.search_program not in valid_search_programs:
            errors.append(f"❌ 不支持的搜索程序 | Unsupported search program: {self.search_program}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True