"""
🧬 Filter ANNOVAR命令 | Filter ANNOVAR Command
基因区域变异提取工具 - 高级优化版本
"""

import click
import sys
import os


def _lazy_import_filter_annovar_main():
    """懒加载Filter ANNOVAR main函数 | Lazy load filter annovar main function"""
    try:
        from ...filter_annovar.main import main as filter_annovar_main
        return filter_annovar_main
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
    short_help='🧬 Filter ANNOVAR：基因区域变异提取工具',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 基因选择参数 | Gene Selection Parameters =====
@click.option('--gene-id', '-g',
              help='🎯 单个基因ID | Single gene ID')
@click.option('--gene-list', '-G',
              type=click.Path(exists=True),
              help='📋 基因ID列表文件 | Gene ID list file (one gene per line)')

# ===== 必需输入文件 | Required Input Files =====
@click.option('--gff',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 GFF注释文件路径 | GFF annotation file path')
@click.option('--exonic',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 ANNOVAR外显子注释文件 | ANNOVAR exonic annotation file (.exonic_variant_function)')
@click.option('--all',
              'all_variant',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 ANNOVAR所有变异注释文件 | ANNOVAR all variants file (.variant_function)')

# ===== 输出控制 | Output Control =====
@click.option('--output', '-o',
              default='./filter_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./filter_output)')
@click.option('--format', '-f',
              type=click.Choice(['excel', 'txt']),
              default='excel',
              help='📊 输出格式 | Output format (default: excel)')

# ===== 分析参数 | Analysis Parameters =====
@click.option('--extend', '-e',
              type=int,
              default=5000,
              help='📏 上下游扩展范围(bp) | Upstream/downstream extension range in bp (default: 5000)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
def filter_annovar(gene_id, gene_list, gff, exonic, all_variant, output, format, extend, threads):
    """
    🧬 Filter ANNOVAR - 基因区域变异提取工具 | Gene Region Variant Extraction Tool
    
    从ANNOVAR注释结果中提取特定基因及其上下游区域的变异信息，
    支持单基因和批量基因分析，自动整合外显子和全基因组变异数据。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🎯 灵活的基因选择:
       • 单基因分析模式
       • 批量基因分析模式
       • 自动基因ID识别
       • 支持多线程并行处理
    
    \b
    📊 完整的变异信息:
       • 外显子区域变异(exonic)
       • 内含子区域变异(intronic)
       • 上游区域变异(upstream)
       • 下游区域变异(downstream)
       • UTR区域变异(UTR5/UTR3)
       • 非编码区域变异(ncRNA)
    
    \b
    📄 多种输出格式:
       • Excel格式(推荐，易于查看)
       • 文本格式(便于后续处理)
       • 自动整理分类
       • 清晰的统计信息
    
    🔄 分析流程 | Analysis Pipeline:
    
    \b
    步骤1: GFF解析 | GFF Parsing
       • 读取GFF注释文件
       • 提取基因位置信息
       • 计算扩展区域范围
       • 建立基因位置索引
    
    \b
    步骤2: 变异提取 | Variant Extraction
       • 读取ANNOVAR外显子注释
       • 读取ANNOVAR所有变异注释
       • 根据基因位置提取变异
       • 包含上下游扩展区域
    
    \b
    步骤3: 变异分类 | Variant Classification
       • 按变异类型分类
       • 按变异位置分类
       • 按功能影响分类
       • 生成分类统计
    
    \b
    步骤4: 结果输出 | Result Output
       • 生成Excel或文本文件
       • 包含详细注释信息
       • 提供统计摘要
       • 保存日志文件
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    GFF文件格式:
       • 标准GFF3格式
       • 包含基因注释信息
       • 包含CDS、mRNA等特征
       • 示例:
         chr1  maker  gene  1000  5000  .  +  .  ID=gene-001
    
    \b
    ANNOVAR外显子文件(.exonic_variant_function):
       • ANNOVAR annotate_variation.pl输出
       • 包含外显子变异注释
       • 格式: 
         line1  nonsynonymous SNV  gene-001:...  chr1  1234  ...
    
    \b
    ANNOVAR所有变异文件(.variant_function):
       • ANNOVAR annotate_variation.pl输出
       • 包含所有区域变异注释
       • 格式:
         exonic  gene-001  chr1  1234  1234  A  T  ...
         intronic  gene-002  chr1  5678  5678  G  C  ...
    
    \b
    基因列表文件(可选):
       • 每行一个基因ID
       • 支持注释行(#开头)
       • 示例:
         gene-001
         gene-002
         # gene-003 (跳过)
         gene-004
    
    📊 输出文件说明 | Output Files:
    
    \b
    Excel格式输出(默认):
       gene-ID_variants.xlsx
       • Sheet1: 外显子变异
       • Sheet2: 内含子变异
       • Sheet3: 上游变异
       • Sheet4: 下游变异
       • Sheet5: 统计摘要
    
    \b
    文本格式输出:
       gene-ID_variants.txt
       • 制表符分隔
       • 包含所有变异信息
       • 按位置排序
    
    \b
    统计信息:
       • 总变异数量
       • 各类型变异数量
       • 各区域变异数量
       • 功能影响分类统计
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 单个基因分析（最常用）
    biopytools filter-annovar -g gene-PHYSODRAFT_288440 \\
        --gff genome.gff \\
        --exonic variants.exonic_variant_function \\
        --all variants.variant_function -o output_dir
    
    \b
    # 📋 批量基因分析
    biopytools filter-annovar -G genes.txt \\
        --gff genome.gff \\
        --exonic variants.exonic_variant_function \\
        --all variants.variant_function -o results
    
    \b
    # 📝 输出文本格式
    biopytools filter-annovar -g gene-ID \\
        --gff genome.gff --exonic exonic.txt --all all.txt \\
        -o output -f txt
    
    \b
    # 📏 自定义扩展范围（扩大到10kb）
    biopytools filter-annovar -g gene-ID \\
        --gff genome.gff --exonic exonic.txt --all all.txt \\
        -e 10000
    
    \b
    # 🚀 多线程批量处理
    biopytools filter-annovar -G gene_list.txt \\
        --gff genome.gff --exonic exonic.txt --all all.txt \\
        -o results -t 64
    
    \b
    # 🎨 完整配置示例
    biopytools filter-annovar -G my_genes.txt \\
        --gff annotation.gff \\
        --exonic sample.exonic_variant_function \\
        --all sample.variant_function \\
        -o gene_variants -f excel -e 8000 -t 44
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 🧬 候选基因变异分析
    • 🔬 功能基因组学研究
    • 🏥 临床相关基因检测
    • 🌱 重要农艺性状基因分析
    • 📊 群体遗传学研究
    • 🧪 基因编辑靶点选择
    • 💊 药物靶标变异分析
    
    ⚙️ 参数详解 | Parameter Details:
    
    \b
    🎯 基因选择:
       --gene-id (-g):
       • 分析单个基因
       • 适合快速查询
       • 基因ID必须与GFF匹配
       • 示例: gene-001, AT1G01010
       
       --gene-list (-G):
       • 批量分析多个基因
       • 文件每行一个基因ID
       • 支持#注释行
       • 自动并行处理
    
    \b
    📏 扩展范围 (--extend):
       • 基因上下游扩展距离
       • 单位: 碱基对(bp)
       • 默认: 5000bp (5kb)
       • 建议值:
         - 启动子分析: 2000-3000bp
         - 调控元件分析: 5000-10000bp
         - 大范围筛选: 10000-20000bp
       
       扩展示例:
       基因: chr1:10000-15000
       扩展5kb后: chr1:5000-20000
    
    \b
    📊 输出格式:
       excel (推荐):
       ✅ 多Sheet组织清晰
       ✅ 易于查看和筛选
       ✅ 包含格式化表格
       ✅ 支持Excel分析工具
       
       txt:
       📝 制表符分隔
       ⚡ 便于后续程序处理
       💾 文件体积更小
       🔧 易于命令行操作
    
    \b
    🧵 线程数:
       • 仅在批量分析时有效
       • 单基因分析不使用多线程
       • 建议: CPU核心数×0.8-1.0
       • 过多线程可能导致I/O瓶颈
    
    📈 变异类型说明 | Variant Type Description:
    
    \b
    外显子变异类型:
       • nonsynonymous SNV: 非同义SNP
       • synonymous SNV: 同义SNP
       • stopgain: 提前终止密码子
       • stoploss: 终止密码子缺失
       • frameshift: 移码突变
       • nonframeshift: 非移码插入/缺失
    
    \b
    基因区域类型:
       • exonic: 外显子区域
       • intronic: 内含子区域
       • upstream: 上游区域
       • downstream: 下游区域
       • UTR5: 5'非翻译区
       • UTR3: 3'非翻译区
       • splicing: 剪接位点
       • ncRNA: 非编码RNA
    
    🔍 结果解读 | Result Interpretation:
    
    \b
    高影响变异:
       • stopgain: 可能导致蛋白截短
       • stoploss: 可能延长蛋白
       • frameshift: 可能完全改变蛋白
       • splicing: 可能影响剪接
    
    \b
    中等影响变异:
       • nonsynonymous SNV: 改变氨基酸
       • nonframeshift indel: 氨基酸插入/缺失
    
    \b
    低影响变异:
       • synonymous SNV: 不改变氨基酸
       • intronic: 通常影响较小
       • UTR: 可能影响表达调控
    
    \b
    调控区域变异:
       • upstream: 可能影响启动子
       • downstream: 可能影响终止子
       • UTR5/UTR3: 可能影响翻译效率
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1️⃣ 基因ID未找到
          → 检查基因ID是否正确
          → 确认GFF文件中有该基因
          → 注意ID格式和大小写
       
       2️⃣ 无变异结果
          → 检查该基因区域是否有变异
          → 尝试增大--extend扩展范围
          → 确认ANNOVAR文件格式正确
       
       3️⃣ 文件格式错误
          → 确认使用ANNOVAR标准输出
          → 检查文件完整性
          → 验证GFF格式规范性
       
       4️⃣ 内存不足
          → 减少并行线程数
          → 分批处理基因列表
          → 使用文本格式减少内存
    
    \b
    调试建议:
       • 先用单基因测试
       • 检查日志文件
       • 验证输入文件格式
       • 使用小数据集测试
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据准备:
       • 确保GFF和ANNOVAR文件来自同一基因组
       • 验证染色体命名一致性
       • 检查文件编码格式(UTF-8)
       • 确认文件完整无损坏
    
    \b
    2️⃣ 参数选择:
       # 启动子区域分析
       -e 3000
       
       # 增强子区域分析
       -e 10000
       
       # 大范围调控元件
       -e 20000
    
    \b
    3️⃣ 批量处理策略:
       # 小批量(<100基因)
       -G genes.txt -t 32
       
       # 大批量(>100基因)
       # 分割基因列表
       split -l 100 genes.txt batch_
       # 分批处理
       for batch in batch_*; do
           biopytools filter-annovar -G $batch ...
       done
    
    \b
    4️⃣ 结果筛选:
       Excel输出后可以:
       • 按变异类型筛选
       • 按功能影响排序
       • 关注高影响变异
       • 结合其他数据库验证
    
    \b
    5️⃣ 下游分析:
       • 使用SnpEff/VEP进一步注释
       • 与dbSNP/ClinVar交叉验证
       • 结合表达数据分析
       • 进行实验验证
    
    📚 相关资源 | Related Resources:
    
    \b
    • ANNOVAR官网: https://annovar.openbioinformatics.org/
    • GFF3格式规范
    • 变异注释最佳实践
    • 功能变异预测工具
    
    \b
    推荐工具:
    • SnpEff: 变异效应预测
    • VEP: 变异注释
    • IGV: 可视化验证
    • SIFT/PolyPhen: 功能预测
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • GFF文件必须包含完整的基因注释
    • ANNOVAR文件必须是标准格式
    • 基因ID必须在GFF中存在
    • 扩展范围不要过大导致重叠
    • Excel格式更适合人工审查
    • 文本格式更适合程序处理
    • 批量分析建议使用多线程
    • 保留日志文件便于追踪问题
    
    🔬 高级应用 | Advanced Applications:
    
    \b
    候选基因优先级排序:
    1. 提取所有候选基因变异
    2. 按功能影响分类
    3. 结合表达数据
    4. 优先分析高影响变异
    
    \b
    调控区域变异分析:
    1. 设置较大扩展范围(-e 10000)
    2. 关注upstream区域变异
    3. 结合转录因子结合位点
    4. 分析潜在调控影响
    
    \b
    比较基因组学应用:
    1. 提取多个物种同源基因变异
    2. 比较变异模式
    3. 识别保守/快速进化区域
    4. 推断功能重要性
    """
    
    # 验证基因选择参数
    if not gene_id and not gene_list:
        raise click.UsageError(
            "❌ 必须指定 --gene-id 或 --gene-list 之一 | "
            "Must specify either --gene-id or --gene-list"
        )
    
    if gene_id and gene_list:
        raise click.UsageError(
            "❌ --gene-id 和 --gene-list 不能同时使用 | "
            "Cannot use both --gene-id and --gene-list"
        )
    
    # 🚀 懒加载
    filter_annovar_main = _lazy_import_filter_annovar_main()
    
    # 构建参数列表
    args = ['filter_annovar.py']
    
    # 基因选择参数
    if gene_id:
        args.extend(['-g', gene_id])
    if gene_list:
        args.extend(['-G', gene_list])
    
    # 必需输入文件
    args.extend(['--gff', gff])
    args.extend(['--exonic', exonic])
    args.extend(['--all', all_variant])
    
    # 输出控制
    if output != './filter_output':
        args.extend(['-o', output])
    if format != 'excel':
        args.extend(['-f', format])
    
    # 分析参数
    if extend != 5000:
        args.extend(['-e', str(extend)])
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        filter_annovar_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 分析执行失败 | Analysis execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv