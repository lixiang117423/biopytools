"""
系统发育树构建命令 | Phylogenetic Tree Builder Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_phylo_main():
    """懒加载系统发育树main函数 | Lazy load phylo main function"""
    try:
        from ...mafft_fasttree.main import main as phylo_main
        return phylo_main
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


@click.command(short_help='🌳 系统发育树构建工具：基于MAFFT和FastTree的进化树分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧬 输入序列文件 (FASTA格式) | Input sequence file (FASTA format)')
@click.option('--output', '-o',
              default='./phylo_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./phylo_output)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--seq-type',
              type=click.Choice(['protein', 'nucleotide']),
              help='🔬 序列类型 (不指定则自动检测) | Sequence type (auto-detect if not specified)')
@click.option('--mafft-params',
              default='--auto',
              help='🧩 MAFFT额外参数 | Additional MAFFT parameters (default: --auto)')
@click.option('--fasttree-params',
              default='',
              help='🌳 FastTree额外参数 | Additional FastTree parameters')
@click.option('--mafft-path',
              default='mafft',
              help='🔧 MAFFT软件路径 | MAFFT software path (default: mafft)')
@click.option('--fasttree-path',
              default='fasttree',
              help='🔧 FastTree软件路径 | FastTree software path (default: fasttree)')
def mafft_fasttree(input, output, threads, seq_type, mafft_params, fasttree_params, 
               mafft_path, fasttree_path):
    """
    🌳 系统发育树构建工具 | Phylogenetic Tree Builder Tool
    
    专业的系统发育树构建流程工具，整合MAFFT多序列比对和FastTree进化树构建，
    支持蛋白质和核酸序列分析，广泛应用于分子进化、系统发育和比较基因组学研究。
    
    ✨ 功能特点 | Features:
    - 🧬 自动序列类型检测
    - 🧩 高效的MAFFT多序列比对
    - 🌳 快速的FastTree进化树构建
    - 🔄 智能序列ID处理
    - 📊 完整的分析日志记录
    - ⚡ 多线程并行加速
    - 🎯 灵活的参数自定义
    
    🔄 分析流程 | Analysis Pipeline:
    1. 📖 检测和验证输入序列类型
    2. 🧹 清理序列和处理序列ID
    3. 🧩 MAFFT多序列比对分析
    4. 🌳 FastTree系统发育树构建
    5. 💾 生成分析结果和日志
    
    🎯 应用场景 | Use Cases:
    - 🌳 分子系统发育分析
    - 🧬 进化关系研究
    - 📊 比较基因组学分析
    - 🔍 序列同源性评估
    - 🧪 蛋白质家族分类
    - 🌿 物种进化树构建
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本系统发育树构建
    biopytools mafft-fasttree -i sequences.fa -o phylo_results
    
    \b
    # 🚀 高性能多线程分析
    biopytools mafft-fasttree -i proteins.fa -o results -t 64
    
    \b
    # ⚙️ 自定义MAFFT和FastTree参数
    biopytools mafft-fasttree -i sequences.fa -o results \\
        --mafft-params "--maxiterate 1000" --fasttree-params "-gamma"
    
    \b
    # 🔬 指定序列类型
    biopytools mafft-fasttree -i nucleotides.fa -o results --seq-type nucleotide
    
    \b
    # 🧪 蛋白质序列系统发育分析
    biopytools mafft-fasttree -i proteins.fasta -o protein_tree \\
        --seq-type protein --threads 96 --mafft-params "--globalpair --maxiterate 1000"
    
    \b
    # 🌿 核酸序列进化树构建
    biopytools mafft-fasttree -i dna_sequences.fa -o dna_tree \\
        --seq-type nucleotide --fasttree-params "-gtr -gamma"
    
    📁 输入文件要求 | Input File Requirements:
    
    🧬 序列文件格式:
    - 📄 标准FASTA格式 (.fa, .fasta, .fas)
    - 🔤 蛋白质序列：标准20个氨基酸字符
    - 🧬 核酸序列：ATCGN字符
    - 🏷️ 序列ID应唯一且有意义
    - 📏 序列数量：至少3条以上
    
    📝 FASTA格式示例:
    >Species_1_gene_X
    MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
    >Species_2_gene_X
    MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    📋 必需参数:
    --input (-i): 🧬 输入FASTA序列文件
        - 包含待分析的序列数据
        - 支持蛋白质和核酸序列
    
    📁 输出控制:
    --output (-o): 📁 结果输出目录
        - 自动创建必要的文件结构
        - 包含比对结果和进化树文件
    
    ⚡ 性能参数:
    --threads (-t): 🧵 并行处理线程数
        - 建议设置为CPU核心数
        - 影响MAFFT比对速度
    
    🔬 序列类型:
    --seq-type: 🧬 指定序列数据类型
        - protein: 蛋白质序列
        - nucleotide: 核酸序列
        - 不指定则自动检测
    
    🧩 MAFFT参数:
    --mafft-params: ⚙️ MAFFT额外参数
        - --auto: 自动选择算法（默认）
        - --globalpair: 全局配对比对
        - --localpair: 局部配对比对
        - --maxiterate N: 最大迭代次数
        - --retree N: 引导树重建次数
    
    🌳 FastTree参数:
    --fasttree-params: 🔧 FastTree额外参数
        - -gamma: 使用Gamma分布速率
        - -gtr: 使用GTR核酸替换模型
        - -lg: 使用LG蛋白质替换模型
        - -wag: 使用WAG蛋白质替换模型
        - -boot N: Bootstrap重抽样次数
    
    🛠️ 工具路径:
    --mafft-path: 🔧 MAFFT可执行文件路径
    --fasttree-path: 🔧 FastTree可执行文件路径
    
    📁 输出文件说明 | Output Files:
    
    📝 主要输出文件:
    - 🧹 cleaned_sequences.fa: 清理后的序列文件
    - 🧩 mafft_alignment.fa: MAFFT比对结果
    - 🌳 phylogenetic_tree.nwk: Newick格式进化树
    - 🗺️ id_mapping.txt: 序列ID映射表
    - 📋 phylo_analysis.log: 详细分析日志
    
    📊 输出目录结构:
    output_directory/
    ├── 🧹 cleaned_sequences.fa
    ├── 🧩 mafft_alignment.fa
    ├── 🌳 phylogenetic_tree.nwk
    ├── 🗺️ id_mapping.txt
    └── 📋 phylo_analysis.log
    
    🌳 进化树文件格式:
    - Newick格式标准系统发育树
    - 可用于FigTree、iTOL等可视化工具
    - 包含分支长度和拓扑结构信息
    
    ⚡ 性能和系统要求 | Performance & System Requirements:
    
    🔧 依赖软件:
    - 🧩 MAFFT (v7.0+): 多序列比对工具
    - 🌳 FastTree (v2.1+): 系统发育树构建工具
    - 🐍 Python (3.7+): 核心运行环境
    
    💻 系统建议:
    - 💾 RAM: 至少2GB，大数据集建议16GB+
    - 🗄️ 存储: 输出目录需要输入文件3倍空间
    - 🧵 CPU: 多核处理器提升MAFFT速度
    
    📊 处理能力估算:
    - 📄 小数据集(<100序列): 几分钟
    - 📋 中等数据集(100-1000序列): 10-60分钟
    - 📁 大数据集(1000-10000序列): 1-10小时
    - 🗂️ 超大数据集(>10000序列): 考虑分批或使用更快算法
    
    💾 内存使用估算:
    - 🔵 基础内存: ~500MB
    - 📈 序列相关: 序列数量×序列长度×2MB
    - 🧵 线程相关: 每线程额外~100MB
    
    🛠️ 故障排除 | Troubleshooting:
    
    ⚠️ 常见问题:
    1. ❌ "MAFFT not found": 检查MAFFT安装和PATH
    2. ❌ "FastTree not found": 检查FastTree安装和PATH
    3. 💾 "内存不足": 减少线程数或序列数量
    4. 📄 "序列格式错误": 验证FASTA格式规范性
    5. 🌳 "树构建失败": 检查比对质量和参数设置
    
    📋 数据质量检查:
    - 🔍 验证FASTA文件格式完整性
    - 📏 检查序列长度一致性
    - 🏷️ 确认序列ID唯一性
    - 🧬 验证序列字符合法性
    
    💡 优化建议:
    - 🎯 根据序列数量调整MAFFT算法
    - 🧵 合理设置线程数平衡速度和内存
    - 📊 大数据集考虑使用--auto以外的快速算法
    - 🌳 根据研究目标选择合适的替换模型
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 使用高质量的序列数据
       - 🔍 预先验证序列格式和完整性
       - 📏 确保序列具有足够的变异性
       - 🏷️ 使用描述性的序列标识符
    
    2️⃣ ⚙️ 参数选择:
       - 🧩 小数据集(<100)使用globalpair获得最佳精度
       - 📊 中等数据集使用--auto平衡速度和精度
       - 🚀 大数据集使用FFT-NS-2等快速算法
       - 🌳 根据序列类型选择合适的替换模型
    
    3️⃣ 🔍 质量控制:
       - 📈 检查比对结果的质量
       - 🌳 验证进化树拓扑结构合理性
       - 🎲 使用Bootstrap评估分支支持度
       - 📚 与已知系统发育关系进行比较
    
    4️⃣ 📦 结果分析:
       - 🎨 使用FigTree等工具可视化进化树
       - 📊 分析分支长度和进化距离
       - 🔍 识别单系群和进化关系
       - 📝 记录分析参数便于重现
    
    🔍 结果验证建议 | Result Validation:
    
    🤖 自动验证:
    - 📊 检查比对覆盖度和gap比例
    - 🌳 验证进化树拓扑结构完整性
    - 📏 评估分支长度分布合理性
    - 🏷️ 确认所有序列都被包含
    
    👥 人工验证:
    - 🎲 检查关键分支的支持度
    - 🔍 验证已知进化关系的正确性
    - 📚 与文献报道的系统发育一致性
    - 🧬 评估外群位置的合理性
    
    📚 相关工具和应用 | Related Tools & Applications:
    
    🔧 上游工具:
    - 📊 序列收集: NCBI, UniProt, Ensembl
    - 🧬 序列编辑: AliView, SeaView, MEGA
    - 🔍 序列过滤: CD-HIT, USEARCH
    
    🔬 下游分析:
    - 🎨 树可视化: FigTree, iTOL, ETE Toolkit
    - 📊 树比较: TreeCmp, RAxML
    - 🧪 进化分析: PAML, HyPhy, BEAST
    - 🌳 树编辑: Dendroscope, TreeGraph
    
    📖 替换模型参考:
    - 🧬 核酸: GTR, HKY, JC69
    - 🧪 蛋白质: LG, WAG, JTT, Blosum62
    - ⚡ 速率异质性: Gamma分布, 不变位点
    
    📚 引用和参考 | Citation & References:
    
    如果在学术研究中使用此工具，请引用:
    - 🧩 MAFFT: Katoh & Standley (2013) Mol Biol Evol
    - 🌳 FastTree: Price et al. (2010) PLoS ONE
    """
    
    # 🚀 懒加载
    phylo_main = _lazy_import_phylo_main()
    
    # 构建参数列表
    args = ['mafft_fasttree.py']
    args.extend(['-i', input])
    
    if output != './phylo_output':
        args.extend(['-o', output])
    if threads != 88:
        args.extend(['-t', str(threads)])
    if seq_type:
        args.extend(['--seq-type', seq_type])
    if mafft_params != '--auto':
        args.extend(['--mafft-params', mafft_params])
    if fasttree_params:
        args.extend(['--fasttree-params', fasttree_params])
    if mafft_path != 'mafft':
        args.extend(['--mafft-path', mafft_path])
    if fasttree_path != 'fasttree':
        args.extend(['--fasttree-path', fasttree_path])
    
    original_argv = sys.argv
    sys.argv = args
    
    try:
        phylo_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 系统发育树构建被用户中断", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 系统发育树构建失败: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv