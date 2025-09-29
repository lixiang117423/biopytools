# """
# OrthoFinder泛基因组分析命令 | OrthoFinder Pangenome Analysis Command
# 高级优化版本：解决--help响应速度问题
# """

# import click
# import sys
# import os


# def _lazy_import_orthofinder_main():
#     """懒加载orthofinder main函数 | Lazy load orthofinder main function"""
#     try:
#         from ...orthofinder_pangenome.main import main as orthofinder_main
#         return orthofinder_main
#     except ImportError as e:
#         click.echo(f"导入错误 | Import Error: {e}", err=True)
#         sys.exit(1)


# def _is_help_request():
#     """检查是否是帮助请求 | Check if this is a help request"""
#     help_flags = {'-h', '--help'}
#     return any(arg in help_flags for arg in sys.argv)


# def _validate_directory_exists(dir_path):
#     """验证目录是否存在（仅在非帮助模式下）| Validate directory existence (only in non-help mode)"""
#     if not _is_help_request() and not os.path.exists(dir_path):
#         raise click.BadParameter(f"目录不存在 | Directory does not exist: {dir_path}")
#     if not _is_help_request() and not os.path.isdir(dir_path):
#         raise click.BadParameter(f"路径不是目录 | Path is not a directory: {dir_path}")
#     return dir_path


# @click.command(short_help='基于OrthoFinder的泛基因组分析工具：同源基因群聚类和基因家族分类',
#                context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# @click.option('--input', '-i',
#               required=True,
#               callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
#               help='输入蛋白质序列文件目录 | Input protein sequence files directory')
# @click.option('--output', '-o',
#               required=True,
#               type=click.Path(),
#               help='输出目录 | Output directory')
# @click.option('--project-name', '-n',
#               type=str,
#               help='项目名称 | Project name')
# @click.option('--soft-threshold',
#               default='all-1',
#               help='Soft基因家族阈值 | Soft gene family threshold (all-1, all-2, specific, or number) (default: all-1)')
# @click.option('--threads', '-t',
#               type=int,
#               default=88,
#               help='线程数 | Number of threads (default: 88)')
# @click.option('--search', '--search-program', '-s',
#               type=click.Choice(['blast', 'blastp', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']),
#               default='blastp',
#               help='序列搜索程序 | Sequence search program (default: blastp)')
# @click.option('--mcl-inflation',
#               type=float,
#               default=1.2,
#               help='MCL inflation参数 | MCL inflation parameter (default: 1.2)')
# @click.option('--dna', '-d',
#               is_flag=True,
#               help='输入序列为DNA序列 | Input sequences are DNA')
# @click.option('--basic-only',
#               is_flag=True,
#               default=True,
#               help='仅进行基础分析 | Perform basic analysis only (default: enabled)')
# @click.option('--generate-trees',
#               is_flag=True,
#               help='生成系统发育树 | Generate phylogenetic trees')
# @click.option('--msa-program',
#               type=click.Choice(['mafft', 'muscle']),
#               default='mafft',
#               help='多序列比对程序 | Multiple sequence alignment program (default: mafft)')
# @click.option('--tree-program',
#               type=click.Choice(['fasttree', 'fasttree_fastest', 'raxml', 'raxml-ng', 'iqtree']),
#               default='fasttree',
#               help='系统发育树构建程序 | Phylogenetic tree inference program (default: fasttree)')
# @click.option('--orthofinder-path',
#               default='orthofinder',
#               help='OrthoFinder程序路径 | OrthoFinder program path (default: orthofinder)')
# @click.option('--resume',
#               is_flag=True,
#               default=True,
#               help='使用已有结果续跑 | Resume from existing results (default: enabled)')
# @click.option('--force',
#               is_flag=True,
#               help='强制重新分析覆盖已有结果 | Force reanalysis overwriting existing results')
# @click.option('--skip-orthofinder',
#               is_flag=True,
#               help='跳过OrthoFinder步骤直接进行分类 | Skip OrthoFinder step and go directly to classification')
# def orthofinder(input, output, project_name, soft_threshold, threads, search, 
#                          mcl_inflation, dna, basic_only, generate_trees, msa_program, 
#                          tree_program, orthofinder_path, resume, force, skip_orthofinder):
#     """
#     基于OrthoFinder的泛基因组分析工具 | OrthoFinder-based Pangenome Analysis Tool
    
#     一个完整的泛基因组分析流程工具，使用OrthoFinder进行同源基因群聚类，
#     然后根据基因在不同基因组中的分布模式将基因家族分类为核心基因、
#     软核心基因、附属基因和特异基因，适用于比较基因组学和进化分析研究。
    
#     功能特点 | Features:
#     - 全自动OrthoFinder同源基因群分析
#     - 智能泛基因组基因家族分类算法
#     - 多种序列搜索算法支持
#     - 灵活的基因家族阈值设置
#     - 详细的统计报告和可视化
#     - 大规模基因组数据处理优化
#     - 可选的系统发育分析功能
    
#     分析流程 | Analysis Pipeline:
#     1. 验证输入蛋白质序列文件完整性
#     2. 运行OrthoFinder同源基因群聚类分析
#     3. 解析同源基因群结果和基因计数数据
#     4. 基于分布模式进行泛基因组分类
#     5. 生成综合统计报告和详细数据表
    
#     泛基因组分类标准 | Pangenome Classification Criteria:
#     - 核心基因 (Core genes): 存在于所有基因组中
#     - 软核心基因 (Soft-core genes): 存在于大多数基因组中
#     - 附属基因 (Accessory genes): 存在于部分基因组中
#     - 特异基因 (Specific genes): 只存在于少数基因组中
    
#     应用场景 | Use Cases:
#     - 细菌和古菌泛基因组分析
#     - 物种进化和基因组可塑性研究
#     - 核心功能基因识别
#     - 基因组比较和系统发育分析
#     - 病原菌毒力因子挖掘
#     - 工业微生物功能基因发现
    
#     示例 | Examples:
    
#     \b
#     # 基本泛基因组分析
#     biopytools orthofinder -i protein_sequences/ -o pangenome_results
    
#     \b
#     # 指定项目名称和自定义soft阈值
#     biopytools orthofinder -i data/proteomes/ -o results/ \\
#         --project-name "Ecoli_pangenome" --soft-threshold all-2
    
#     \b
#     # 使用Diamond搜索算法和高线程数
#     biopytools orthofinder -i sequences/ -o output/ \\
#         -s diamond -t 64 --soft-threshold 5
    
#     \b
#     # 数值型soft阈值设置
#     biopytools orthofinder -i proteomes/ -o results/ \\
#         --soft-threshold 8 --search diamond_ultra_sens
    
#     \b
#     # 包含系统发育树分析的完整流程
#     biopytools orthofinder -i sequences/ -o comprehensive_analysis/ \\
#         --generate-trees --msa-program muscle --tree-program iqtree \\
#         --search diamond -t 96
    
#     \b
#     # DNA序列分析模式
#     biopytools orthofinder -i dna_sequences/ -o dna_results/ \\
#         --dna --search blast_nucl --mcl-inflation 1.5
    
#     \b
#     # 高性能服务器优化配置
#     biopytools orthofinder -i large_dataset/ -o hpc_results/ \\
#         -t 128 --search mmseqs --mcl-inflation 1.1 \\
#         --soft-threshold all-3 --project-name "HighThroughput_Analysis"
    
#     输入文件要求 | Input File Requirements:
    
#     蛋白质序列文件:
#     - 标准FASTA格式 (.fa, .faa, .fasta)
#     - 每个基因组一个独立文件
#     - 文件名将作为基因组标识符
#     - 建议使用描述性文件名
#     - 序列ID需要在同一基因组内唯一
    
#     DNA序列文件 (使用--dna时):
#     - 标准FASTA格式核酸序列
#     - 通常为CDS序列或基因序列
#     - 需要正确的开放阅读框
#     - 建议预先进行基因预测
    
#     文件组织示例:
#     protein_sequences/
#     ├── genome1.faa
#     ├── genome2.faa  
#     ├── genome3.faa
#     └── genome4.faa
    
#     序列格式示例:
#     >gene001 hypothetical protein
#     MKILVFASLLSLLAAGVQAAPEAQPVIKLEEATGKGQRWLWAGLASRLVDAKMQTIHSASLVS
#     >gene002 DNA polymerase
#     MAKTFEKVKLAASQAGEEAATQSNAQLQMKLVLKATTQADQAIKAKTQASLAALNQADEQLAS
    
#     参数详细说明 | Parameter Details:
    
#     必需参数:
#     --input: 包含所有基因组蛋白质序列文件的目录
#     --output: 结果输出目录，将自动创建
    
#     泛基因组分类参数:
#     --soft-threshold: 软核心基因阈值设置
#         - "all-1": 所有基因组减1 (N-1)
#         - "all-2": 所有基因组减2 (N-2)  
#         - "specific": 基于数据自动确定
#         - 数值: 直接指定基因组数量阈值
    
#     OrthoFinder参数:
#     --threads: 并行计算线程数，建议设置为CPU核心数
#     --search: 序列相似性搜索算法
#         - blastp: 标准BLAST蛋白质搜索（默认）
#         - diamond: 快速Diamond搜索
#         - diamond_ultra_sens: Diamond高敏感模式
#         - mmseqs: MMseqs2快速搜索
#         - blast_nucl: 核酸序列BLAST（DNA模式）
#     --mcl-inflation: MCL聚类算法膨胀参数，影响聚类粒度
    
#     系统发育分析参数:
#     --generate-trees: 启用系统发育树构建
#     --msa-program: 多序列比对工具选择
#     --tree-program: 系统发育树构建方法
    
#     序列类型:
#     --dna: 指定输入为DNA序列而非蛋白质序列
#     --basic-only: 仅进行基础OrthoFinder分析，跳过高级功能
    
#     输出文件说明 | Output Files:
    
#     核心结果文件:
#     - pangenome_classification.txt: 基因家族分类结果
#     - gene_families_table.txt: 详细的基因家族数据表
#     - comprehensive_report.txt: 综合分析报告
#     - orthofinder_results/: OrthoFinder原始输出目录
    
#     详细输出结构:
#     output_directory/
#     ├── pangenome_classification.txt
#     ├── gene_families_table.txt  
#     ├── comprehensive_report.txt
#     ├── orthofinder_results/
#     │   ├── Orthogroups/
#     │   │   ├── Orthogroups.txt
#     │   │   └── Orthogroups.GeneCount.txt
#     │   ├── Phylogenetic_Hierarchical_Orthogroups/
#     │   └── Species_Tree/
#     └── analysis.log
    
#     性能和系统要求 | Performance & System Requirements:
    
#     依赖软件:
#     - OrthoFinder (v2.3.0+): 同源基因群分析核心工具
#     - BLAST+ 或 Diamond: 序列相似性搜索
#     - MCL: 马尔可夫聚类算法
#     - FastTree/RAxML/IQ-TREE: 系统发育树构建（可选）
#     - MAFFT/MUSCLE: 多序列比对（可选）
    
#     系统建议:
#     - RAM: 至少8GB，大数据集推荐32GB+
#     - CPU: 多核处理器，推荐16核+
#     - 存储: 至少输入数据大小的5倍自由空间
#     - 临时空间: 大量临时文件，建议使用SSD
    
#     数据规模估算:
#     - 小规模(<50基因组): ~4GB内存，~20GB存储
#     - 中规模(50-200基因组): ~16GB内存，~100GB存储  
#     - 大规模(200-500基因组): ~64GB内存，~500GB存储
#     - 超大规模(>500基因组): 考虑分批处理
    
#     算法复杂度:
#     - 时间复杂度: O(n²m) (n=基因组数, m=平均基因数)
#     - 空间复杂度: O(nm)
#     - Diamond模式可显著提升速度
    
#     故障排除 | Troubleshooting:
    
#     常见问题:
#     1. "OrthoFinder not found": 检查OrthoFinder安装和PATH
#     2. "Memory error": 增加系统内存或减少线程数
#     3. "Empty orthogroups": 检查输入序列质量和格式
#     4. "MCL inflation error": 调整--mcl-inflation参数
#     5. "Tree generation failed": 检查MSA和树构建软件
    
#     优化建议:
#     - 大数据集优先使用Diamond搜索
#     - 根据内存限制调整线程数
#     - 预先验证输入文件格式
#     - 使用SSD存储提高I/O性能
#     - 监控系统资源使用情况
    
#     最佳实践 | Best Practices:
    
#     1. 数据准备:
#        - 使用高质量的基因组注释
#        - 确保蛋白质序列完整性
#        - 统一序列命名规范
#        - 移除冗余和低质量序列
    
#     2. 参数选择:
#        - 小数据集使用BLAST获得最佳精度
#        - 大数据集使用Diamond平衡速度与精度
#        - 根据研究目标调整soft阈值
#        - 考虑基因组进化距离选择参数
    
#     3. 结果解释:
#        - 核心基因通常为基础代谢功能
#        - 附属基因可能与环境适应相关
#        - 特异基因可能包含新功能或水平转移
#        - 结合功能注释进行生物学解释
    
#     4. 质量控制:
#        - 检查基因组完整性和质量
#        - 验证同源基因群的合理性
#        - 比较不同参数设置的结果
#        - 使用已知基因验证分类准确性
    
#     引用和参考 | Citation & References:
    
#     如果在学术研究中使用此工具，请引用相关文献:
#     - OrthoFinder: Emms & Kelly (2019) Genome Biology
#     - Diamond: Buchfink et al. (2015) Nature Methods
#     - MCL: Van Dongen (2000) PhD Thesis
#     - 泛基因组概念: Tettelin et al. (2005) PNAS
#     """
    
#     # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
#     orthofinder_main = _lazy_import_orthofinder_main()
    
#     # 构建参数列表传递给原始main函数 | Build argument list for original main function
#     args = ['orthofinder.py']
    
#     # 必需参数 | Required parameters
#     args.extend(['-i', input])
#     args.extend(['-o', output])
    
#     # 可选参数（只有在非默认值时才添加）| Optional parameters (add only when non-default)
#     if project_name:
#         args.extend(['-n', project_name])
    
#     if soft_threshold != 'all-1':
#         args.extend(['--soft-threshold', soft_threshold])
    
#     if threads != 88:
#         args.extend(['-t', str(threads)])
    
#     if search != 'blastp':
#         args.extend(['-s', search])
    
#     if mcl_inflation != 1.2:
#         args.extend(['--mcl-inflation', str(mcl_inflation)])
    
#     if msa_program != 'mafft':
#         args.extend(['--msa-program', msa_program])
    
#     if tree_program != 'fasttree':
#         args.extend(['--tree-program', tree_program])
    
#     if orthofinder_path != 'orthofinder':
#         args.extend(['--orthofinder-path', orthofinder_path])
    
#     # 断点续跑参数 | Resume parameters  
#     if not resume:  # 默认是True，当设置为False时添加参数
#         # 注意：这里可能需要根据main.py中的实际参数名调整
#         pass  # 如果main.py中有--no-resume参数的话
    
#     if force:
#         args.append('--force')
    
#     if skip_orthofinder:
#         args.append('--skip-orthofinder')
    
#     # 布尔选项 | Boolean options
#     if dna:
#         args.append('-d')
    
#     if not basic_only:  # 默认是True，所以当设置为False时需要处理
#         # 这里需要根据实际main.py中的参数处理方式来调整
#         pass
    
#     if generate_trees:
#         args.append('--generate-trees')
    
#     # 保存并恢复sys.argv | Save and restore sys.argv
#     original_argv = sys.argv
#     sys.argv = args
    
#     try:
#         # 调用原始的main函数 | Call original main function
#         orthofinder_main()
#     except SystemExit as e:
#         # 处理程序正常退出 | Handle normal program exit
#         if e.code != 0:
#             sys.exit(e.code)
#     except KeyboardInterrupt:
#         click.echo("\nOrthoFinder泛基因组分析被用户中断 | OrthoFinder pangenome analysis interrupted by user", err=True)
#         sys.exit(1)
#     except Exception as e:
#         click.echo(f"OrthoFinder泛基因组分析失败 | OrthoFinder pangenome analysis failed: {e}", err=True)
#         sys.exit(1)
#     finally:
#         sys.argv = original_argv

"""
🧬 OrthoFinder泛基因组分析命令 | OrthoFinder Pangenome Analysis Command
⚡ 高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_orthofinder_main():
    """⚡ 懒加载orthofinder main函数 | Lazy load orthofinder main function"""
    try:
        from ...orthofinder_pangenome.main import main as orthofinder_main
        return orthofinder_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """❓ 检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_directory_exists(dir_path):
    """🔍 验证目录是否存在（仅在非帮助模式下）| Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(dir_path):
        raise click.BadParameter(f"❌ 目录不存在 | Directory does not exist: {dir_path}")
    if not _is_help_request() and not os.path.isdir(dir_path):
        raise click.BadParameter(f"❌ 路径不是目录 | Path is not a directory: {dir_path}")
    return dir_path


@click.command(short_help='🧬 基于OrthoFinder的泛基因组分析工具：同源基因群聚类和基因家族分类',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='📁 输入蛋白质序列文件目录 | Input protein sequence files directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📁 输出目录 | Output directory')
@click.option('--project-name', '-n',
              type=str,
              help='🏷️  项目名称 | Project name')
@click.option('--soft-threshold',
              default='all-1',
              help='🟠 Soft基因家族阈值 | Soft gene family threshold (all-1, all-2, specific, or number) (default: all-1)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--search', '--search-program', '-s',
              type=click.Choice(['blast', 'blastp', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']),
              default='blastp',
              help='🔍 序列搜索程序 | Sequence search program (default: blastp)')
@click.option('--mcl-inflation',
              type=float,
              default=1.2,
              help='🔧 MCL inflation参数 | MCL inflation parameter (default: 1.2)')
@click.option('--dna', '-d',
              is_flag=True,
              help='🧬 输入序列为DNA序列 | Input sequences are DNA')
@click.option('--basic-only',
              is_flag=True,
              default=True,
              help='🔧 仅进行基础分析 | Perform basic analysis only (default: enabled)')
@click.option('--generate-trees',
              is_flag=True,
              help='🌳 生成系统发育树 | Generate phylogenetic trees')
@click.option('--msa-program',
              type=click.Choice(['mafft', 'muscle']),
              default='mafft',
              help='🔤 多序列比对程序 | Multiple sequence alignment program (default: mafft)')
@click.option('--tree-program',
              type=click.Choice(['fasttree', 'fasttree_fastest', 'raxml', 'raxml-ng', 'iqtree']),
              default='fasttree',
              help='🌳 系统发育树构建程序 | Phylogenetic tree inference program (default: fasttree)')
@click.option('--orthofinder-path',
              default='orthofinder',
              help='🔧 OrthoFinder程序路径 | OrthoFinder program path (default: orthofinder)')
@click.option('--resume',
              is_flag=True,
              default=True,
              help='🔄 使用已有结果续跑 | Resume from existing results (default: enabled)')
@click.option('--force',
              is_flag=True,
              help='🔄 强制重新分析覆盖已有结果 | Force reanalysis overwriting existing results')
@click.option('--skip-orthofinder',
              is_flag=True,
              help='⏭️  跳过OrthoFinder步骤直接进行分类 | Skip OrthoFinder step and go directly to classification')
def orthofinder(input, output, project_name, soft_threshold, threads, search, 
                         mcl_inflation, dna, basic_only, generate_trees, msa_program, 
                         tree_program, orthofinder_path, resume, force, skip_orthofinder):
    """
    🧬 基于OrthoFinder的泛基因组分析工具 | OrthoFinder-based Pangenome Analysis Tool
    
    📋 一个完整的泛基因组分析流程工具，使用OrthoFinder进行同源基因群聚类，
    然后根据基因在不同基因组中的分布模式将基因家族分类为核心基因、
    软核心基因、附属基因和特异基因，适用于比较基因组学和进化分析研究。
    
    ✨ 功能特点 | Features:
    - 🤖 全自动OrthoFinder同源基因群分析
    - 🎯 智能泛基因组基因家族分类算法
    - 🔍 多种序列搜索算法支持
    - ⚙️ 灵活的基因家族阈值设置
    - 📊 详细的统计报告和可视化
    - 🚀 大规模基因组数据处理优化
    - 🌳 可选的系统发育分析功能
    
    🔄 分析流程 | Analysis Pipeline:
    1. ✅ 验证输入蛋白质序列文件完整性
    2. 🧬 运行OrthoFinder同源基因群聚类分析
    3. 📥 解析同源基因群结果和基因计数数据
    4. 🎯 基于分布模式进行泛基因组分类
    5. 📝 生成综合统计报告和详细数据表
    
    🎯 泛基因组分类标准 | Pangenome Classification Criteria:
    - 🔴 核心基因 (Core genes): 存在于所有基因组中
    - 🟠 软核心基因 (Soft-core genes): 存在于大多数基因组中
    - 🟡 附属基因 (Accessory genes): 存在于部分基因组中
    - 🟢 特异基因 (Specific genes): 只存在于少数基因组中
    
    🔬 应用场景 | Use Cases:
    - 🦠 细菌和古菌泛基因组分析
    - 🧬 物种进化和基因组可塑性研究
    - 🔴 核心功能基因识别
    - 📊 基因组比较和系统发育分析
    - 🦠 病原菌毒力因子挖掘
    - 🏭 工业微生物功能基因发现
    
    💡 示例 | Examples:
    
    \b
    # 🔰 基本泛基因组分析
    biopytools orthofinder -i protein_sequences/ -o pangenome_results
    
    \b
    # 🏷️ 指定项目名称和自定义soft阈值
    biopytools orthofinder -i data/proteomes/ -o results/ \\
        --project-name "Ecoli_pangenome" --soft-threshold all-2
    
    \b
    # ⚡ 使用Diamond搜索算法和高线程数
    biopytools orthofinder -i sequences/ -o output/ \\
        -s diamond -t 64 --soft-threshold 5
    
    \b
    # 🔢 数值型soft阈值设置
    biopytools orthofinder -i proteomes/ -o results/ \\
        --soft-threshold 8 --search diamond_ultra_sens
    
    \b
    # 🌳 包含系统发育树分析的完整流程
    biopytools orthofinder -i sequences/ -o comprehensive_analysis/ \\
        --generate-trees --msa-program muscle --tree-program iqtree \\
        --search diamond -t 96
    
    \b
    # 🧬 DNA序列分析模式
    biopytools orthofinder -i dna_sequences/ -o dna_results/ \\
        --dna --search blast_nucl --mcl-inflation 1.5
    
    \b
    # 🖥️ 高性能服务器优化配置
    biopytools orthofinder -i large_dataset/ -o hpc_results/ \\
        -t 128 --search mmseqs --mcl-inflation 1.1 \\
        --soft-threshold all-3 --project-name "HighThroughput_Analysis"
    
    📁 输入文件要求 | Input File Requirements:
    
    🧬 蛋白质序列文件:
    - 📄 标准FASTA格式 (.fa, .faa, .fasta)
    - 📦 每个基因组一个独立文件
    - 🏷️ 文件名将作为基因组标识符
    - 💡 建议使用描述性文件名
    - 🆔 序列ID需要在同一基因组内唯一
    
    🧬 DNA序列文件 (使用--dna时):
    - 📄 标准FASTA格式核酸序列
    - 🧬 通常为CDS序列或基因序列
    - ✅ 需要正确的开放阅读框
    - 💡 建议预先进行基因预测
    
    📁 文件组织示例:
    protein_sequences/
    ├── genome1.faa
    ├── genome2.faa  
    ├── genome3.faa
    └── genome4.faa
    
    📄 序列格式示例:
    >gene001 hypothetical protein
    MKILVFASLLSLLAAGVQAAPEAQPVIKLEEATGKGQRWLWAGLASRLVDAKMQTIHSASLVS
    >gene002 DNA polymerase
    MAKTFEKVKLAASQAGEEAATQSNAQLQMKLVLKATTQADQAIKAKTQASLAALNQADEQLAS
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    ⭐ 必需参数:
    --input: 📁 包含所有基因组蛋白质序列文件的目录
    --output: 📁 结果输出目录，将自动创建
    
    🎯 泛基因组分类参数:
    --soft-threshold: 🟠 软核心基因阈值设置
        - "all-1": 所有基因组减1 (N-1)
        - "all-2": 所有基因组减2 (N-2)  
        - "specific": 基于数据自动确定
        - 数值: 直接指定基因组数量阈值
    
    🔧 OrthoFinder参数:
    --threads: 🧵 并行计算线程数，建议设置为CPU核心数
    --search: 🔍 序列相似性搜索算法
        - blastp: 标准BLAST蛋白质搜索（默认）
        - diamond: 快速Diamond搜索
        - diamond_ultra_sens: Diamond高敏感模式
        - mmseqs: MMseqs2快速搜索
        - blast_nucl: 核酸序列BLAST（DNA模式）
    --mcl-inflation: 🔧 MCL聚类算法膨胀参数，影响聚类粒度
    
    🌳 系统发育分析参数:
    --generate-trees: 🌳 启用系统发育树构建
    --msa-program: 🔤 多序列比对工具选择
    --tree-program: 🌳 系统发育树构建方法
    
    🧬 序列类型:
    --dna: 🧬 指定输入为DNA序列而非蛋白质序列
    --basic-only: 🔧 仅进行基础OrthoFinder分析，跳过高级功能
    
    📄 输出文件说明 | Output Files:
    
    📋 核心结果文件:
    - 📊 pangenome_classification.txt: 基因家族分类结果
    - 📋 gene_families_table.txt: 详细的基因家族数据表
    - 📝 comprehensive_report.txt: 综合分析报告
    - 📁 orthofinder_results/: OrthoFinder原始输出目录
    
    📁 详细输出结构:
    output_directory/
    ├── 📊 pangenome_classification.txt
    ├── 📋 gene_families_table.txt  
    ├── 📝 comprehensive_report.txt
    ├── 📁 orthofinder_results/
    │   ├── 📁 Orthogroups/
    │   │   ├── 📄 Orthogroups.txt
    │   │   └── 📊 Orthogroups.GeneCount.txt
    │   ├── 📁 Phylogenetic_Hierarchical_Orthogroups/
    │   └── 🌳 Species_Tree/
    └── 📋 analysis.log
    
    🖥️ 性能和系统要求 | Performance & System Requirements:
    
    🔧 依赖软件:
    - 🧬 OrthoFinder (v2.3.0+): 同源基因群分析核心工具
    - 🔍 BLAST+ 或 Diamond: 序列相似性搜索
    - 🔗 MCL: 马尔可夫聚类算法
    - 🌳 FastTree/RAxML/IQ-TREE: 系统发育树构建（可选）
    - 🔤 MAFFT/MUSCLE: 多序列比对（可选）
    
    💻 系统建议:
    - 🧠 RAM: 至少8GB，大数据集推荐32GB+
    - 🔢 CPU: 多核处理器，推荐16核+
    - 💾 存储: 至少输入数据大小的5倍自由空间
    - ⚡ 临时空间: 大量临时文件，建议使用SSD
    
    📊 数据规模估算:
    - 🔰 小规模(<50基因组): ~4GB内存，~20GB存储
    - 🔶 中规模(50-200基因组): ~16GB内存，~100GB存储  
    - 🔴 大规模(200-500基因组): ~64GB内存，~500GB存储
    - 🚀 超大规模(>500基因组): 考虑分批处理
    
    ⚙️ 算法复杂度:
    - ⏰ 时间复杂度: O(n²m) (n=基因组数, m=平均基因数)
    - 💾 空间复杂度: O(nm)
    - ⚡ Diamond模式可显著提升速度
    
    🔧 故障排除 | Troubleshooting:
    
    ❓ 常见问题:
    1. ❌ "OrthoFinder not found": 检查OrthoFinder安装和PATH
    2. 💾 "Memory error": 增加系统内存或减少线程数
    3. 📭 "Empty orthogroups": 检查输入序列质量和格式
    4. 🔧 "MCL inflation error": 调整--mcl-inflation参数
    5. 🌳 "Tree generation failed": 检查MSA和树构建软件
    
    💡 优化建议:
    - ⚡ 大数据集优先使用Diamond搜索
    - 🧵 根据内存限制调整线程数
    - ✅ 预先验证输入文件格式
    - 💿 使用SSD存储提高I/O性能
    - 📊 监控系统资源使用情况
    
    🏆 最佳实践 | Best Practices:
    
    1. 📋 数据准备:
       - ✅ 使用高质量的基因组注释
       - 🔍 确保蛋白质序列完整性
       - 🏷️ 统一序列命名规范
       - 🧹 移除冗余和低质量序列
    
    2. ⚙️ 参数选择:
       - 🎯 小数据集使用BLAST获得最佳精度
       - ⚡ 大数据集使用Diamond平衡速度与精度
       - 🎯 根据研究目标调整soft阈值
       - 🧬 考虑基因组进化距离选择参数
    
    3. 📊 结果解释:
       - 🔴 核心基因通常为基础代谢功能
       - 🟡 附属基因可能与环境适应相关
       - 🟢 特异基因可能包含新功能或水平转移
       - 📚 结合功能注释进行生物学解释
    
    4. ✅ 质量控制:
       - 🔍 检查基因组完整性和质量
       - ✅ 验证同源基因群的合理性
       - 🔄 比较不同参数设置的结果
       - 🧪 使用已知基因验证分类准确性
    
    📚 引用和参考 | Citation & References:
    
    📖 如果在学术研究中使用此工具，请引用相关文献:
    - 🧬 OrthoFinder: Emms & Kelly (2019) Genome Biology
    - ⚡ Diamond: Buchfink et al. (2015) Nature Methods
    - 🔗 MCL: Van Dongen (2000) PhD Thesis
    - 🧬 泛基因组概念: Tettelin et al. (2005) PNAS
    """
    
    # ⚡ 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    orthofinder_main = _lazy_import_orthofinder_main()
    
    # 🔧 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['orthofinder.py']
    
    # ⭐ 必需参数 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    
    # ⚙️ 可选参数（只有在非默认值时才添加）| Optional parameters (add only when non-default)
    if project_name:
        args.extend(['-n', project_name])
    
    if soft_threshold != 'all-1':
        args.extend(['--soft-threshold', soft_threshold])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if search != 'blastp':
        args.extend(['-s', search])
    
    if mcl_inflation != 1.2:
        args.extend(['--mcl-inflation', str(mcl_inflation)])
    
    if msa_program != 'mafft':
        args.extend(['--msa-program', msa_program])
    
    if tree_program != 'fasttree':
        args.extend(['--tree-program', tree_program])
    
    if orthofinder_path != 'orthofinder':
        args.extend(['--orthofinder-path', orthofinder_path])
    
    # 🔄 断点续跑参数 | Resume parameters  
    if not resume:  # 默认是True，当设置为False时添加参数
        # 注意：这里可能需要根据main.py中的实际参数名调整
        pass  # 如果main.py中有--no-resume参数的话
    
    if force:
        args.append('--force')
    
    if skip_orthofinder:
        args.append('--skip-orthofinder')
    
    # ✅ 布尔选项 | Boolean options
    if dna:
        args.append('-d')
    
    if not basic_only:  # 默认是True，所以当设置为False时需要处理
        # 这里需要根据实际main.py中的参数处理方式来调整
        pass
    
    if generate_trees:
        args.append('--generate-trees')
    
    # 💾 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 🚀 调用原始的main函数 | Call original main function
        orthofinder_main()
    except SystemExit as e:
        # ✅ 处理程序正常退出 | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n⛔ OrthoFinder泛基因组分析被用户中断 | OrthoFinder pangenome analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ OrthoFinder泛基因组分析失败 | OrthoFinder pangenome analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv