"""
INDEL PAV分析命令 | INDEL PAV Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_pav_main():
    """懒加载PAV main函数 | Lazy load PAV main function"""
    try:
        from ...indel_pav.main import main as pav_main
        return pav_main
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


@click.command(short_help='🧬 INDEL PAV分析工具：基于VCF文件的插入缺失存在/缺失变异检测',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📁 输入VCF文件路径 (支持压缩和未压缩) | Input VCF file path (supports compressed and uncompressed)')
@click.option('--output', '-o',
              default='./indel_pav.txt',
              type=click.Path(),
              help='📁 输出文件路径 | Output file path (default: ./indel_pav.txt)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🚀 线程数 | Number of threads (default: 88)')
@click.option('--min-length',
              type=int,
              default=1,
              help='📏 最小INDEL长度 (bp) | Minimum INDEL length (bp) (default: 1)')
@click.option('--max-length',
              type=int,
              help='📏 最大INDEL长度 (bp) | Maximum INDEL length (bp)')
@click.option('--min-quality', '-q',
              type=float,
              default=20.0,
              help='⭐ 最小质量分数 | Minimum quality score (default: 20.0)')
@click.option('--min-depth', '-d',
              type=int,
              default=5,
              help='📊 最小深度 | Minimum depth (default: 5)')
@click.option('--max-missing',
              type=float,
              default=0.8,
              help='❓ 最大缺失率 (0-1) | Maximum missing rate (0-1) (default: 0.8)')
@click.option('--include-complex',
              is_flag=True,
              help='🔀 包含复杂变异 | Include complex variants')
@click.option('--compress',
              is_flag=True,
              help='🗜️  压缩输出文件 | Compress output file')
@click.option('--bcftools-path',
              default='bcftools',
              help='🔧 BCFtools软件路径 | BCFtools software path (default: bcftools)')
def indelpav(vcf, output, threads, min_length, max_length, min_quality, min_depth, 
              max_missing, include_complex, compress, bcftools_path):
    """
    🧬 INDEL PAV分析工具 | INDEL PAV Analysis Tool
    
    一个专业的插入缺失(INDEL)存在/缺失变异(PAV)分析工具，用于从VCF文件中
    识别和分析结构变异模式，支持群体遗传学研究中的基因组多样性评估和
    适应性进化分析。
    
    ✨ 功能特点 | Features:
    - 🔍 高效的VCF文件INDEL变异识别
    - 📊 灵活的质量和深度过滤参数
    - 🧬 支持简单和复杂结构变异分析
    - 📈 详细的PAV统计和分布报告
    - ⚡ 多线程并行处理优化
    - 🗜️ 支持压缩输入输出格式
    - 📋 标准化的结果输出格式
    
    🔄 分析流程 | Analysis Pipeline:
    1. 📖 解析VCF文件并提取INDEL变异
    2. 🔍 根据长度和质量参数过滤变异
    3. 🧬 计算样本间的PAV矩阵
    4. 📊 生成统计报告和汇总信息
    5. 📝 输出标准化的PAV结果文件
    
    📝 PAV分析概念 | PAV Analysis Concept:
    PAV (Presence/Absence Variation) 是指基因组中某些序列片段在不同个体中
    存在或缺失的变异模式，是重要的结构变异类型，与适应性进化、疾病易感性
    和农作物性状等密切相关。
    
    🎯 应用场景 | Use Cases:
    - 🔬 群体遗传学和进化分析
    - 🌾 农作物品种改良和性状挖掘
    - 🏥 人类疾病关联性研究
    - 🧬 比较基因组学分析
    - 📊 基因组多样性评估
    - 🔍 适应性进化检测
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本PAV分析
    biopytools indelpav -v variants.vcf -o pav_results.txt
    
    \b
    # ⚙️ 自定义过滤参数
    biopytools indelpav -v data.vcf.gz -o filtered_pav.txt \\
        --min-length 5 --max-length 1000 -q 30 -d 10
    
    \b
    # 🚀 高性能处理大文件
    biopytools indelpav -v large_variants.vcf.gz -o results.txt \\
        -t 64 --compress --max-missing 0.5
    
    \b
    # 🔀 包含复杂变异的全面分析
    biopytools indelpav -v complex.vcf -o comprehensive.txt \\
        --include-complex --min-quality 25 --min-depth 8
    
    \b
    # 📊 严格质控的高质量分析
    biopytools indelpav -v filtered.vcf.gz -o high_quality.txt \\
        -q 40 -d 15 --max-missing 0.2 --min-length 10
    
    \b
    # 🔧 指定外部工具路径
    biopytools indelpav -v variants.vcf -o results.txt \\
        --bcftools-path /usr/local/bin/bcftools
    
    📁 输入文件要求 | Input File Requirements:
    
    🧬 VCF文件格式:
    - 📋 标准VCF 4.0+格式
    - 🗜️ 支持压缩文件 (.vcf.gz)
    - 🔍 包含INDEL变异信息
    - 📊 包含质量和深度信息
    - 🏷️ 正确的样本标识符
    
    ✅ 必需字段:
    - 📍 CHROM, POS: 染色体和位置信息
    - 🆔 ID: 变异标识符 (可为'.')
    - 🧬 REF, ALT: 参考和替代序列
    - ⭐ QUAL: 变异质量分数
    - 🔧 FORMAT: 样本格式字段
    - 📊 样本基因型和深度信息
    
    📝 VCF格式示例:
    ##fileformat=VCFv4.2
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2
    chr1    1000    .       ATCG    A       35.2    PASS    .       GT:DP   1/1:12  0/1:8
    chr1    2000    .       G       GTCA    42.8    PASS    .       GT:DP   0/0:15  1/1:10
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    📋 必需参数:
    --vcf: 📁 输入VCF文件，包含INDEL变异信息
    
    📁 输出控制:
    --output: 📁 PAV结果输出文件路径
    --compress: 🗜️ 启用输出文件压缩 (gzip格式)
    
    🔍 变异过滤参数:
    --min-length: 📏 INDEL最小长度阈值，过滤短变异
    --max-length: 📏 INDEL最大长度阈值，过滤超长变异
    --min-quality: ⭐ 最小质量分数，提高变异可靠性
    --min-depth: 📊 最小测序深度，确保检测置信度
    --max-missing: ❓ 最大样本缺失率，控制数据完整性
    
    🧬 变异类型控制:
    --include-complex: 🔀 包含复杂结构变异
        - 默认只分析简单INDEL
        - 启用后包含多等位基因和复杂重排
    
    ⚡ 性能参数:
    --threads: 🚀 并行处理线程数
        - 建议设置为CPU核心数
        - 影响VCF解析和计算速度
    
    🔧 工具路径:
    --bcftools-path: 🔧 BCFtools可执行文件路径
        - 用于VCF文件处理和索引
        - 支持压缩文件快速访问
    
    📁 输出文件说明 | Output Files:
    
    📝 主要输出:
    - 📄 {output_file}: PAV矩阵结果文件
    - 📊 {output_file}.summary: 统计汇总报告
    - 📋 {output_file}.log: 详细分析日志
    
    📊 PAV矩阵格式:
    第一列：变异ID (CHROM:POS:REF:ALT)
    后续列：各样本的存在/缺失状态
    - 1: 存在 (Present)
    - 0: 缺失 (Absent)  
    - NA: 数据缺失
    
    📈 统计报告内容:
    - 🔢 总变异数量和样本数量
    - 📊 变异长度分布统计
    - 📈 样本间PAV相似性矩阵
    - 🎯 核心/特异/稀有变异分类
    - 📍 染色体分布统计
    
    ⚡ 性能和系统要求 | Performance & System Requirements:
    
    🔧 依赖软件:
    - 🧰 BCFtools (v1.9+): VCF文件处理
    - 🐍 Python 3.7+: 核心运行环境
    - 📚 标准Python库: pandas, numpy
    
    💻 系统建议:
    - 💾 RAM: 至少4GB，大文件建议16GB+
    - 🗄️ 存储: 输出目录需要输入文件2倍空间
    - 🚀 CPU: 多核处理器提升并行处理效率
    
    📊 性能指标:
    - 📄 小文件(<100MB): 通常几分钟内完成
    - 📋 中等文件(100MB-1GB): 10-30分钟
    - 📁 大文件(1GB-10GB): 30分钟-2小时
    - 🗂️ 超大文件(>10GB): 考虑分批处理
    
    💾 内存使用估算:
    - 🔵 基础内存: ~500MB
    - 📈 文件相关: ~VCF文件大小的2-3倍
    - 🧵 多线程: 每线程额外~100MB
    
    🛠️ 故障排除 | Troubleshooting:
    
    ⚠️ 常见问题:
    1. ❌ "BCFtools not found": 检查BCFtools安装和PATH设置
    2. 💾 "内存不足": 减少线程数或增加系统内存
    3. 📁 "VCF格式错误": 验证输入文件格式规范性
    4. 🔒 "权限拒绝": 检查输出目录写入权限
    5. ⚡ "处理速度慢": 调整线程数和过滤参数
    
    📋 数据质量检查:
    - 🔍 验证VCF文件头信息完整性
    - 📊 检查样本基因型缺失率
    - ⭐ 确认质量分数分布合理
    - 📏 验证INDEL长度分布
    
    💡 优化建议:
    - 🎯 合理设置过滤阈值平衡质量和数量
    - 🚀 使用压缩格式减少I/O时间
    - 📊 预先统计VCF文件基本信息
    - 🧵 根据内存限制调整线程数
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 预先验证VCF文件完整性和格式
       - 🔍 检查样本标识符一致性
       - 📊 统计变异类型和质量分布
       - 💾 确保足够的存储空间
    
    2️⃣ ⚙️ 参数调优:
       - 🎯 根据研究目标设置长度阈值
       - ⭐ 基于数据质量选择质量分数
       - 📊 考虑样本量调整缺失率阈值
       - 🧬 评估是否需要包含复杂变异
    
    3️⃣ 🔍 质量控制:
       - 📈 检查PAV矩阵的合理性
       - 🎲 随机验证部分变异的准确性
       - 📊 比较不同参数设置的结果
       - 🔗 与已知变异数据库进行比较
    
    4️⃣ 📦 结果解释:
       - 🧬 结合生物学知识解释PAV模式
       - 📍 关注特定基因组区域的变异
       - 🌳 考虑群体结构对PAV分布的影响
       - 📚 参考相关文献验证发现
    
    🔍 结果验证 | Result Validation:
    
    🤖 自动检查:
    - 📊 PAV矩阵维度和数据类型
    - 🔢 变异和样本数量统计
    - 📈 缺失数据比例检查
    - ✅ 输出文件格式验证
    
    👥 手动验证:
    - 🎲 随机抽查变异的基因型调用
    - 🔍 检查已知变异的PAV状态
    - 📊 验证统计结果的生物学合理性
    - 🌐 与公共数据库变异频率比较
    
    📚 引用和参考 | Citation & References:
    
    📖 相关方法和工具:
    - 🧰 BCFtools: Danecek et al. (2021) GigaScience
    - 📊 PAV分析: Abel et al. (2020) Nature Genetics
    - 🧬 结构变异: Sudmant et al. (2015) Nature
    - 📈 群体遗传学: Excoffier & Lischer (2010) Mol Ecol Resour
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    pav_main = _lazy_import_pav_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['indelpav.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-v', vcf])
    
    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    if output != './indel_pav.txt':
        args.extend(['-o', output])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if min_length != 1:
        args.extend(['--min-length', str(min_length)])
    
    if max_length is not None:
        args.extend(['--max-length', str(max_length)])
    
    if min_quality != 20.0:
        args.extend(['-q', str(min_quality)])
    
    if min_depth != 5:
        args.extend(['-d', str(min_depth)])
    
    if max_missing != 0.8:
        args.extend(['--max-missing', str(max_missing)])
    
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
    
    # 布尔选项 🚩 | Boolean options
    if include_complex:
        args.append('--include-complex')
    
    if compress:
        args.append('--compress')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        pav_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 INDEL PAV分析被用户中断 | INDEL PAV analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 INDEL PAV分析失败 | INDEL PAV analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv