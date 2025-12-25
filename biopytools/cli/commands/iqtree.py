"""
🌲 IQ-TREE系统发育树分析命令 | IQ-TREE Phylogenetic Tree Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_iqtree_main():
    """懒加载IQ-TREE main函数 | Lazy load IQ-TREE main function"""
    try:
        from ...iqtree.main import main as iqtree_main
        return iqtree_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"❌ 文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='🌲 IQ-TREE系统发育树分析：最大似然法构建高精度进化树',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 输入比对文件 | Input alignment file (FASTA/PHYLIP/NEXUS)')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📂 输出目录 | Output directory')
@click.option('--prefix', '-p',
              required=True,
              help='🏷️  输出文件前缀 | Output file prefix')

# ===== 核心参数 | Core Parameters =====
@click.option('--model', '-m',
              help='🧬 进化模型 | Evolutionary model (auto-select if not specified)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')

# ===== Bootstrap参数 | Bootstrap Parameters =====
@click.option('--bootstrap', '-b',
              type=int,
              default=1000,
              help='🔄 Bootstrap重复次数 | Bootstrap replicates (default: 1000)')
@click.option('--boot-type',
              type=click.Choice(['ufboot', 'standard']),
              default='ufboot',
              help='📊 Bootstrap类型 | Bootstrap type (default: ufboot)')
@click.option('--save-boot-trees',
              is_flag=True,
              help='💾 保存所有bootstrap树 | Save all bootstrap trees to file')

# ===== 树选项 | Tree Options =====
@click.option('--outgroup',
              help='🎯 外群名称 | Outgroup taxon names (comma-separated)')
@click.option('--constraint',
              type=click.Path(exists=True),
              help='🔗 约束树文件 | Constraint tree file')

# ===== 分区分析 | Partition Analysis =====
@click.option('--partition',
              type=click.Path(exists=True),
              help='📊 分区文件 | Partition file')
@click.option('--partition-mode',
              type=click.Choice(['edge-linked', 'edge-equal', 'edge-unlinked']),
              default='edge-linked',
              help='🔀 分区模式 | Partition mode (default: edge-linked)')

# ===== 高级功能 | Advanced Features =====
@click.option('--concordance',
              type=click.Path(exists=True),
              help='🔬 一致性因子分析基因树文件 | Concordance factor gene tree file')
@click.option('--ancestral',
              is_flag=True,
              help='🧬 启用祖先状态重建 | Enable ancestral state reconstruction')

# ===== 其他参数 | Other Parameters =====
@click.option('--seed',
              type=int,
              help='🎲 随机种子 | Random seed')
@click.option('--runs',
              type=int,
              default=1,
              help='🔁 独立运行次数 | Number of independent runs (default: 1)')
@click.option('--redo',
              is_flag=True,
              help='🔄 重新运行分析 | Redo analysis')
@click.option('--iqtree-path',
              default='iqtree',
              help='🔧 IQ-TREE程序路径 | IQ-TREE program path (default: iqtree)')
def iqtree(input, output, prefix, model, threads, bootstrap, boot_type, save_boot_trees,
           outgroup, constraint, partition, partition_mode,
           concordance, ancestral, seed, runs, redo, iqtree_path):
    """
    🌲 IQ-TREE系统发育树分析工具 | IQ-TREE Phylogenetic Tree Analysis Tool
    
    基于最大似然法构建高精度系统发育树，支持多种进化模型、Bootstrap评估、
    分区分析、一致性因子分析和祖先状态重建等高级功能。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🔬 核心功能:
       • 最大似然法建树: 高精度进化树构建
       • 模型自动选择: ModelFinder自动选择最佳模型
       • 超快Bootstrap: UFBoot提供快速可靠的支持度评估
       • 分区分析: 支持多基因/多位点联合分析
       • 一致性因子: 评估基因树与物种树的一致性
       • 祖先重建: 推断祖先序列和性状状态
    
    \b
    🚀 性能优化:
       • 多线程并行计算
       • 智能算法优化
       • 内存高效处理
       • 支持大规模数据集
    
    \b
    🎯 质量保证:
       • 自动模型检验
       • Bootstrap置信度评估
       • SH-aLRT快速似然比检验
       • 多重统计检验
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    比对文件格式:
       • FASTA: 最常用格式，*.fasta, *.fa
       • PHYLIP: 系统发育标准格式，*.phy
       • NEXUS: 包含额外信息，*.nex
       • CLUSTAL: ClustalW输出，*.aln
       • MAF: 多序列比对格式
    
    \b
    比对质量要求:
       • 序列已经完成多序列比对
       • 移除低质量区域或gap过多的列
       • 序列名称唯一且不含特殊字符
       • 建议使用MAFFT或MUSCLE进行比对
    
    \b
    分区文件格式(可选):
       # NEXUS格式
       #NEXUS
       begin sets;
           charset gene1 = 1-500;
           charset gene2 = 501-1000;
           charset gene3 = 1001-1500;
       end;
       
       # RAxML格式
       DNA, gene1 = 1-500
       DNA, gene2 = 501-1000
       DNA, gene3 = 1001-1500
    
    📊 输出文件说明 | Output Files:
    
    \b
    主要输出文件:
       • *.treefile - 最佳进化树(Newick格式)
       • *.iqtree - 详细分析报告
       • *.log - 运行日志
       • *.model.gz - 模型选择结果
       • *.bionj - NJ初始树
    
    \b
    Bootstrap相关:
       • *.ufboot - UFBoot支持度值
       • *.contree - 共识树
       • *.splits.nex - Split网络
    
    \b
    高级分析输出:
       • *.cf.tree - 一致性因子标注树
       • *.cf.stat - 一致性统计
       • *.state - 祖先状态重建结果
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本建树分析
    biopytools iqtree -i alignment.fasta -o tree_results -p mytree
    
    \b
    # 🧬 指定进化模型
    biopytools iqtree -i protein.phy -o results -p phylo \\
        -m LG+G8+F -t 32
    
    \b
    # 🔄 使用标准Bootstrap
    biopytools iqtree -i nucleotide.fa -o output -p tree \\
        -b 1000 --boot-type standard
    
    \b
    # 📊 分区分析(多基因)
    biopytools iqtree -i concat.fasta -o partition_out -p multigene \\
        --partition partitions.nex --partition-mode edge-linked
    
    \b
    # 🎯 指定外群
    biopytools iqtree -i aligned.fa -o rooted_tree -p myphy \\
        --outgroup "species1,species2"
    
    \b
    # 🔬 一致性因子分析
    biopytools iqtree -i species.aln -o concordance -p cf_tree \\
        --concordance gene_trees.nwk -t 64
    
    \b
    # 🧬 完整分析流程(包含祖先重建)
    biopytools iqtree -i full_alignment.fasta -o complete_analysis -p full \\
        --partition parts.txt --concordance genes.tre --ancestral \\
        -b 1000 --seed 12345
    
    \b
    # 🚀 高性能配置
    biopytools iqtree -i large_dataset.phy -o fast_run -p speedtest \\
        -t 128 --runs 10 -m GTR+I+G
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 物种进化关系研究
    • 基因家族进化分析
    • 分子钟定年分析
    • 群体遗传学研究
    • 病原体溯源分析
    • 保护生物学研究
    • 比较基因组学
    
    🧬 进化模型选择 | Model Selection:
    
    \b
    核苷酸模型:
       • JC/F81: 简单模型
       • HKY/TN: 考虑转换/颠换差异
       • GTR: 通用时间可逆模型(最常用)
       • +I: 不变位点
       • +G: Gamma速率异质性
       • +F: 经验碱基频率
    
    \b
    蛋白质模型:
       • LG: 一般蛋白质进化
       • WAG: 全局蛋白质替换
       • JTT: Jones-Taylor-Thornton
       • Blosum62: 基于BLOSUM矩阵
       • +C: 位点异质性分类
    
    \b
    推荐组合:
       • DNA: GTR+I+G (通用)
       • 蛋白质: LG+G8+F (通用)
       • 密码子: GY+F3X4 (编码区)
       • 不确定时: 使用ModelFinder自动选择
    
    ⚡ 性能要求 | Performance Requirements:
    
    \b
    硬件建议:
       • CPU: 16核以上推荐(支持多线程)
       • 内存: 8GB+(大数据集需32GB+)
       • 存储: 输入文件的10-20倍空间
       • 计算时间: 取决于序列数和长度
    
    \b
    性能估算:
       • 100序列×1000bp: ~10分钟(32线程)
       • 500序列×5000bp: ~2小时(64线程)
       • 1000序列×10000bp: ~10小时(128线程)
       • Bootstrap会增加计算时间
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1. 似然值不收敛 → 检查比对质量，尝试不同模型
       2. Bootstrap值过低 → 数据信号弱，增加位点数
       3. 内存不足 → 减少线程数或使用分区分析
       4. 树拓扑不稳定 → 增加Bootstrap次数，使用约束树
    
    \b
    优化建议:
       • 移除gap过多的位点(>50% gap)
       • 使用信息量高的位点
       • 合理设置分区策略
       • 长枝吸引可尝试移除问题序列
       • 多次独立运行检验一致性
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据准备
       • 使用高质量的序列比对(MAFFT/MUSCLE)
       • 去除gap过多或异常的序列
       • 检查序列方向一致性
       • 验证没有重复序列
    
    \b
    2️⃣ 模型选择
       • 小数据集: 让ModelFinder自动选择
       • 大数据集: 先用子集测试模型
       • 分区数据: 考虑每个分区独立模型
       • 蛋白质: 优先考虑LG或WAG模型
    
    \b
    3️⃣ Bootstrap评估
       • UFBoot: 通用选择，快速且准确
       • 标准Bootstrap: 保守评估需求
       • 建议1000次以上获得稳定支持度
       • 可结合SH-aLRT检验
    
    \b
    4️⃣ 结果验证
       • 检查似然值收敛性
       • 评估Bootstrap支持度分布
       • 对比不同运行结果的一致性
       • 使用TreeViewer可视化检查
    
    \b
    5️⃣ 发表标准
       • 报告使用的进化模型
       • 说明Bootstrap设置
       • 提供关键节点支持度
       • 存档比对文件和树文件
    
    📚 相关资源 | Related Resources:
    
    \b
    • IQ-TREE官网: http://www.iqtree.org/
    • 在线文档: http://www.iqtree.org/doc/
    • 进化模型: http://www.iqtree.org/doc/Substitution-Models
    • 教程视频: http://www.iqtree.org/workshop/
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 确保输入序列已正确比对
    • 模型选择会显著影响结果
    • Bootstrap值<70%支持度较弱
    • 大数据集建议使用分区分析
    • 多次运行验证拓扑稳定性
    • 长枝吸引可能导致错误分组
    • 保存所有输出文件便于重现分析
    """
    
    # 🚀 懒加载
    iqtree_main = _lazy_import_iqtree_main()
    
    # 构建参数列表
    args = ['iqtree.py']
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-p', prefix])
    
    # 核心参数
    if model:
        args.extend(['-m', model])
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # Bootstrap参数
    if bootstrap != 1000:
        args.extend(['-b', str(bootstrap)])
    if boot_type != 'ufboot':
        args.extend(['--boot-type', boot_type])
    if save_boot_trees:
        args.append('--save-boot-trees')
    
    # 树选项
    if outgroup:
        args.extend(['--outgroup', outgroup])
    if constraint:
        args.extend(['--constraint', constraint])
    
    # 分区分析
    if partition:
        args.extend(['--partition', partition])
    if partition_mode != 'edge-linked':
        args.extend(['--partition-mode', partition_mode])
    
    # 高级功能
    if concordance:
        args.extend(['--concordance', concordance])
    if ancestral:
        args.append('--ancestral')
    
    # 其他参数
    if seed:
        args.extend(['--seed', str(seed)])
    if runs != 1:
        args.extend(['--runs', str(runs)])
    if redo:
        args.append('--redo')
    if iqtree_path != 'iqtree':
        args.extend(['--iqtree-path', iqtree_path])
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        iqtree_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 IQ-TREE分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 IQ-TREE分析执行失败 | Analysis execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv