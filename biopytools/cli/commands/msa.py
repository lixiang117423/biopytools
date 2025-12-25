"""
🧬 多序列比对分析命令 | Multiple Sequence Alignment Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_msa_main():
    """懒加载MSA main函数 | Lazy load MSA main function"""
    try:
        from ...msa_toolkit.main import main as msa_main
        return msa_main
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
    short_help='🧬 多序列比对工具：支持MAFFT/Clustal Omega/MUSCLE/T-Coffee',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📁 输入序列文件 | Input sequence file (FASTA format)')
@click.option('--output', '-o',
              required=True,
              help='📝 输出文件前缀 | Output file prefix')

# ===== 比对方法 | Alignment Method =====
@click.option('--method', '-m',
              type=click.Choice(['mafft', 'clustalo', 'muscle', 't_coffee']),
              default='mafft',
              help='🔧 比对方法 | Alignment method (default: mafft)')

# ===== 通用参数 | Common Parameters =====
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--format', '-f',
              type=click.Choice(['fasta', 'clustal', 'phylip', 'nexus']),
              default='fasta',
              help='💾 输出格式 | Output format (default: fasta)')

# ===== MAFFT参数 | MAFFT Parameters =====
@click.option('--mafft-strategy',
              type=click.Choice(['auto', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi']),
              default='auto',
              help='🎯 MAFFT比对策略 | MAFFT alignment strategy (default: auto)')
@click.option('--mafft-maxiterate',
              type=int,
              default=1000,
              help='🔄 MAFFT最大迭代次数 | MAFFT max iterations (default: 1000)')

# ===== Clustal Omega参数 | Clustal Omega Parameters =====
@click.option('--clustalo-iterations',
              type=int,
              default=0,
              help='🔄 Clustal Omega迭代次数 | Clustal Omega iterations (default: 0)')

# ===== MUSCLE参数 | MUSCLE Parameters =====
@click.option('--muscle-maxiters',
              type=int,
              default=16,
              help='🔄 MUSCLE最大迭代次数 | MUSCLE max iterations (default: 16)')

# ===== 工具路径 | Tool Paths =====
@click.option('--mafft-path',
              default='mafft',
              help='📍 MAFFT程序路径 | MAFFT program path (default: mafft)')
@click.option('--clustalo-path',
              default='clustalo',
              help='📍 Clustal Omega程序路径 | Clustal Omega program path (default: clustalo)')
@click.option('--muscle-path',
              default='muscle',
              help='📍 MUSCLE程序路径 | MUSCLE program path (default: muscle)')
@click.option('--tcoffee-path',
              default='t_coffee',
              help='📍 T-Coffee程序路径 | T-Coffee program path (default: t_coffee)')
def msa(input, output, method, threads, format,
        mafft_strategy, mafft_maxiterate,
        clustalo_iterations, muscle_maxiters,
        mafft_path, clustalo_path, muscle_path, tcoffee_path):
    """
    🧬 多序列比对分析工具 | Multiple Sequence Alignment Analysis Tool
    
    整合多种主流比对算法，为蛋白质和核酸序列提供高质量的多序列比对服务，
    支持MAFFT、Clustal Omega、MUSCLE和T-Coffee等多种比对方法。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🔬 多种比对算法:
       • MAFFT: 快速准确的比对算法(推荐)
       • Clustal Omega: 经典的高质量比对
       • MUSCLE: 快速高效的比对方法
       • T-Coffee: 高精度的一致性比对
    
    \b
    🚀 性能优化:
       • 多线程并行处理
       • 自动选择最佳策略
       • 支持大规模数据集
       • 渐进式比对算法
    
    \b
    🎯 高质量输出:
       • 多种输出格式支持
       • 详细的比对统计
       • 质量评估报告
       • 可视化比对结果
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    序列文件格式:
       • 必须是FASTA格式
       • 每条序列有唯一的ID
       • 支持蛋白质和核酸序列
       • 自动检测序列类型
    
    \b
    序列命名规范:
       >sequence_id1 [description]
       ATCGATCGATCG...
       >sequence_id2 [description]
       ATCGATCGATCG...
    
    \b
    数据质量要求:
       • 序列ID唯一且无特殊字符
       • 建议序列数: 3-10000条
       • 序列长度: 50bp-50000bp
       • 避免过短或过长的序列
    
    📊 输出文件说明 | Output Files:
    
    \b
    主要输出文件:
       • *.aln.fasta - 比对后的序列(FASTA格式)
       • *.aln.clustal - Clustal格式比对结果
       • *.stats.txt - 比对统计信息
       • *.log - 详细运行日志
    
    \b
    统计信息包含:
       • 序列数量和平均长度
       • 比对长度和gap比例
       • 保守位点数量
       • 一致性分数
       • 质量评估指标
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本用法 - 使用MAFFT
    biopytools msa -i sequences.fasta -o alignment
    
    \b
    # 🔧 选择比对方法
    biopytools msa -i seqs.fa -o result -m clustalo -t 32
    
    \b
    # 🎯 MAFFT高精度比对
    biopytools msa -i input.fasta -o output -m mafft \\
        --mafft-strategy linsi --threads 64
    
    \b
    # 🔄 Clustal Omega迭代优化
    biopytools msa -i sequences.fa -o aligned -m clustalo \\
        --clustalo-iterations 5 -t 32
    
    \b
    # 💪 MUSCLE快速比对
    biopytools msa -i seqs.fasta -o muscle_result -m muscle \\
        --muscle-maxiters 32 -t 16
    
    \b
    # 📊 指定输出格式
    biopytools msa -i input.fa -o output -f phylip -m mafft
    
    \b
    # 🚀 高性能大数据集比对
    biopytools msa -i large_dataset.fasta -o big_align \\
        -m mafft --mafft-strategy fftns -t 128
    
    \b
    # 🎨 T-Coffee精确比对
    biopytools msa -i proteins.fa -o tcoffee_out -m t_coffee -t 16
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 系统发育树构建前的序列比对
    • 蛋白质家族分析
    • 保守结构域识别
    • 功能位点预测
    • 进化分析研究
    • 基因组注释辅助
    • motif识别和分析
    
    🔧 比对方法选择指南 | Method Selection Guide:
    
    \b
    MAFFT (推荐 | Recommended):
       ✅ 优点: 速度快、精度高、稳定性好
       📊 适用: 大多数场景，特别是大数据集
       🎯 策略选择:
          • auto: 自动选择(推荐)
          • linsi: 最高精度(小数据集<200序列)
          • ginsi: 高精度，适合全局相似序列
          • einsi: 高精度，适合包含大插入/删除
          • fftns: 快速(大数据集>10000序列)
          • fftnsi: 快速迭代版本
    
    \b
    Clustal Omega:
       ✅ 优点: 经典算法、结果可靠、文献广泛
       📊 适用: 中等规模数据集(100-5000序列)
       🔄 迭代次数: 0(快速)到5(高质量)
    
    \b
    MUSCLE:
       ✅ 优点: 速度快、内存效率高
       📊 适用: 中小规模数据集(<1000序列)
       🔄 迭代次数: 16(标准)到32(高质量)
    
    \b
    T-Coffee:
       ✅ 优点: 最高精度、一致性好
       ⚠️ 缺点: 速度较慢、内存消耗大
       📊 适用: 小数据集(<100序列)，需要最高质量
    
    \b
    快速选择指南:
       • <100序列 → MAFFT linsi 或 T-Coffee
       • 100-1000序列 → MAFFT auto 或 Clustal Omega
       • 1000-10000序列 → MAFFT auto
       • >10000序列 → MAFFT fftns
       • 需要最高精度 → MAFFT linsi 或 T-Coffee
       • 需要快速结果 → MAFFT fftns 或 MUSCLE
    
    ⚡ 性能要求 | Performance Requirements:
    
    \b
    硬件建议:
       • CPU: 8核以上推荐
       • 内存: 
         - 小数据集(<1000序列): 4GB
         - 中等数据集(1000-5000): 16GB
         - 大数据集(>5000): 32GB+
       • 存储: 输入文件的5-10倍空间
    
    \b
    时间估算(MAFFT auto, 32线程):
       • 100序列×500bp: ~1分钟
       • 500序列×1000bp: ~5分钟
       • 1000序列×2000bp: ~20分钟
       • 5000序列×1000bp: ~1小时
       • 10000序列×500bp: ~2小时
    
    \b
    软件依赖:
       • Python 3.7+
       • Biopython
       • 至少一种比对工具:
         - MAFFT 7.0+ (推荐)
         - Clustal Omega 1.2+
         - MUSCLE 3.8+
         - T-Coffee 11.0+
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1. 内存不足 → 减少线程数或使用更快的策略
       2. 速度太慢 → 选择更快的算法或策略
       3. 结果质量差 → 使用更精确的策略或迭代
       4. 程序找不到 → 指定完整的工具路径
    
    \b
    优化建议:
       • 预先过滤低质量序列
       • 移除冗余序列(如>95%相似度)
       • 考虑分批次比对大数据集
       • 使用SSD提升I/O性能
       • 合理设置线程数(不超过物理核心)
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据准备
       • 检查序列质量和完整性
       • 统一序列方向(DNA/RNA)
       • 移除过短或过长的异常序列
       • 确保序列ID唯一且规范
    
    \b
    2️⃣ 方法选择
       • 根据数据规模选择合适算法
       • 考虑精度和速度的平衡
       • 对于关键分析，尝试多种方法对比
       • 参考文献中常用的方法
    
    \b
    3️⃣ 参数优化
       • 小数据集: 使用高精度设置
       • 大数据集: 优先考虑速度
       • 迭代次数: 根据质量需求调整
       • 线程数: 设置为物理核心数
    
    \b
    4️⃣ 质量控制
       • 检查比对的gap比例(<30%)
       • 评估保守区域的一致性
       • 可视化检查比对结果
       • 使用工具评估比对质量(如TrimAl)
    
    \b
    5️⃣ 结果处理
       • 移除gap过多的列(可选)
       • 保存多种格式便于后续分析
       • 记录使用的参数便于重现
       • 可视化检查比对质量
    
    📚 相关资源 | Related Resources:
    
    \b
    • MAFFT: https://mafft.cbrc.jp/alignment/software/
    • Clustal Omega: http://www.clustal.org/omega/
    • MUSCLE: https://www.drive5.com/muscle/
    • T-Coffee: http://www.tcoffee.org/
    • Biopython: https://biopython.org/
    
    \b
    可视化工具:
    • Jalview: 功能强大的比对编辑器
    • AliView: 轻量级比对查看器
    • MEGA: 综合进化分析软件
    • Geneious: 商业生物信息学平台
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 比对前确保序列已去除引物和接头
    • 超大数据集(>10000序列)建议使用MAFFT fftns
    • T-Coffee适合小数据集，大数据集会非常慢
    • 定期更新比对软件获得最佳性能
    • 保存原始未比对序列便于重新分析
    • 比对质量直接影响下游系统发育分析
    • 对于发表文章，明确报告使用的方法和参数
    """
    
    # 🚀 懒加载
    msa_main = _lazy_import_msa_main()
    
    # 构建参数列表
    args = ['msa.py']
    args.extend(['-i', input])
    args.extend(['-o', output])
    
    # 比对方法
    if method != 'mafft':
        args.extend(['-m', method])
    
    # 通用参数
    if threads != 88:
        args.extend(['-t', str(threads)])
    if format != 'fasta':
        args.extend(['-f', format])
    
    # MAFFT参数
    if mafft_strategy != 'auto':
        args.extend(['--mafft-strategy', mafft_strategy])
    if mafft_maxiterate != 1000:
        args.extend(['--mafft-maxiterate', str(mafft_maxiterate)])
    
    # Clustal Omega参数
    if clustalo_iterations != 0:
        args.extend(['--clustalo-iterations', str(clustalo_iterations)])
    
    # MUSCLE参数
    if muscle_maxiters != 16:
        args.extend(['--muscle-maxiters', str(muscle_maxiters)])
    
    # 工具路径
    if mafft_path != 'mafft':
        args.extend(['--mafft-path', mafft_path])
    if clustalo_path != 'clustalo':
        args.extend(['--clustalo-path', clustalo_path])
    if muscle_path != 'muscle':
        args.extend(['--muscle-path', muscle_path])
    if tcoffee_path != 't_coffee':
        args.extend(['--tcoffee-path', tcoffee_path])
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        msa_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 多序列比对被用户中断 | MSA interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 多序列比对执行失败 | MSA execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv