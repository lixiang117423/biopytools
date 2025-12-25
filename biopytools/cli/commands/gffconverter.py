"""
GFF格式转换命令 | GFF Format Conversion Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_gff_main():
    """懒加载GFF main函数 | Lazy load GFF main function"""
    try:
        from ...gff_converter.main import main as gff_main
        return gff_main
    except ImportError as e:
        click.echo(f"导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help='GFF格式转换工具：标准化基因ID格式和物种注释',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入GFF文件路径 | Input GFF file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出GFF文件路径 | Output GFF file path')
@click.option('--species-name', '-s',
              required=True,
              type=str,
              help='物种名称 (如: OV53) | Species name (e.g., OV53)')
@click.option('--species-prefix', '-p',
              required=True,
              type=str,
              help='物种缩写 (如: Ov) | Species prefix (e.g., Ov)')
@click.option('--start-num',
              type=int,
              default=10,
              help='起始编号 | Starting number for gene numbering (default: 10)')
@click.option('--step',
              type=int,
              default=10,
              help='编号步长 | Step size for gene numbering (default: 10)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='线程数 | Number of threads (default: 88)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式 | Verbose output mode')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件 | Keep intermediate files')
@click.option('--show-sample',
              type=int,
              metavar='N',
              help='显示N个转换示例然后退出 | Show N conversion samples and exit')
def gffconverter(input, output, species_name, species_prefix, start_num, step, 
                 threads, verbose, keep_intermediate, show_sample):
    """
    GFF格式转换工具 | GFF Format Conversion Tool
    
    一个专业的GFF/GFF3格式标准化工具，用于统一基因注释文件的ID命名规范，
    支持自定义物种标识符和编号规则，适用于多物种比较基因组学研究和
    基因组数据库标准化处理。
    
    功能特点 | Features:
    - 智能GFF/GFF3文件格式解析和验证
    - 灵活的基因ID重命名和标准化
    - 自定义物种标识符和编号系统
    - 保持原始注释信息完整性
    - 支持大规模基因组数据处理
    - 详细的转换日志和统计报告
    - 多线程并行处理优化
    
    转换流程 | Conversion Pipeline:
    1. 解析输入GFF文件并验证格式
    2. 提取原始基因特征和属性信息
    3. 根据指定规则生成新的基因ID
    4. 更新特征属性和交叉引用关系
    5. 格式化输出并保持GFF标准兼容
    
    ID命名规范 | ID Naming Convention:
    基因ID格式: {species_prefix}g{padded_number}
    示例: 物种前缀"Ov" + 起始编号10 + 步长10 = Ovg010, Ovg020, Ovg030...
    
    应用场景 | Use Cases:
    - 多物种基因组比较分析预处理
    - 基因组数据库标准化和整合
    - 注释文件格式转换和清理
    - 基因ID冲突解决和重编号
    - 系统发育和比较基因组学研究
    - 生物信息学流程标准化
    
    示例 | Examples:
    
    \b
    # 基本格式转换
    biopytools renamegff -i input.gff -o output.gff \\
        -s OV53 -p Ov
    
    \b
    # 自定义编号规则
    biopytools renamegff -i annotation.gff -o standardized.gff \\
        -s "Ecoli_K12" -p Ec --start-num 100 --step 5
    
    \b
    # 显示转换示例（预览模式）
    biopytools renamegff -i data.gff3 -o output.gff \\
        -s Species01 -p Sp --show-sample 10
    
    \b
    # 详细模式转换并保留中间文件
    biopytools renamegff -i complex.gff -o processed.gff \\
        -s OV53 -p Ov --verbose --keep-intermediate
    
    \b
    # 高性能批处理模式
    biopytools renamegff -i large_genome.gff -o converted.gff \\
        -s LargeGenome -p Lg -t 64 --start-num 1 --step 1
    
    \b
    # 多物种项目标准化示例
    biopytools renamegff -i species_A.gff -o std_species_A.gff \\
        -s SpeciesA -p SpA --start-num 1000 --step 10
    biopytools renamegff -i species_B.gff -o std_species_B.gff \\
        -s SpeciesB -p SpB --start-num 2000 --step 10
    
    输入文件要求 | Input File Requirements:
    
    GFF/GFF3文件格式:
    - 标准GFF3格式 (推荐) 或GFF2格式
    - 包含完整的基因特征注释
    - 正确的九列制表符分隔格式
    - 有效的属性字段 (第9列)
    - 支持压缩文件 (.gz)
    
    必需特征类型:
    - gene: 基因特征记录
    - mRNA/transcript: 转录本信息 (可选)
    - CDS/exon: 编码序列和外显子 (可选)
    
    标准GFF3格式示例:
    ##gff-version 3
    chr1    source  gene    1000    5000    .    +    .    ID=gene-001;Name=hypothetical_protein
    chr1    source  mRNA    1000    5000    .    +    .    ID=mRNA-001;Parent=gene-001
    chr1    source  CDS     1200    4800    .    +    0    ID=CDS-001;Parent=mRNA-001
    
    参数详细说明 | Parameter Details:
    
    必需参数:
    --input: 输入GFF文件路径，支持.gff、.gff3、.gtf格式
    --output: 输出文件路径，将生成标准化的GFF3文件
    --species-name: 完整物种名称，用于文件头和文档记录
    --species-prefix: 物种缩写，用于生成基因ID前缀
    
    编号控制参数:
    --start-num: 基因编号起始数值，默认从10开始
    --step: 编号递增步长，默认每次增加10
        - 建议使用10的倍数便于后续插入新基因
        - 小步长(1-5)适用于密集编号需求
        - 大步长(50-100)适用于稀疏编号或预留空间
    
    性能参数:
    --threads: 并行处理线程数，影响大文件处理速度
    
    输出控制:
    --verbose: 启用详细日志输出，显示转换详情
    --keep-intermediate: 保留处理过程中的临时文件
    --show-sample: 预览模式，显示指定数量的转换示例
    
    输出文件说明 | Output Files:
    
    主要输出:
    - {output_file}: 标准化后的GFF3文件
    - {output_file}.log: 详细的转换日志
    - conversion_stats.txt: 转换统计报告 (verbose模式)
    
    中间文件 (--keep-intermediate时保留):
    - parsed_features.tmp: 解析后的特征数据
    - id_mappings.txt: 原始ID到新ID的映射表
    - validation_report.txt: 格式验证报告
    
    输出格式特点:
    - 符合GFF3标准规范
    - 保持原始生物学信息完整
    - 统一的基因ID命名格式
    - 正确的层次结构关系
    - 兼容主流基因组浏览器
    
    ID生成规则 | ID Generation Rules:
    
    基因ID格式: {prefix}g{number}
    - prefix: 用户指定的物种前缀
    - g: 固定的基因标识符
    - number: 零填充的序号
    
    编号示例:
    prefix="Ov", start_num=10, step=10:
    - 第1个基因: Ovg010
    - 第2个基因: Ovg020  
    - 第3个基因: Ovg030
    
    prefix="Ec", start_num=1, step=1:
    - 第1个基因: Ecg001
    - 第2个基因: Ecg002
    - 第3个基因: Ecg003
    
    相关特征ID:
    - mRNA: {gene_id}.t1, {gene_id}.t2 (多转录本)
    - CDS: {mrna_id}.cds
    - protein: {gene_id}.p1, {gene_id}.p2
    
    性能和系统要求 | Performance & System Requirements:
    
    依赖要求:
    - Python 3.7+
    - 标准库模块：pathlib, argparse, logging
    - 无外部软件依赖，纯Python实现
    
    系统建议:
    - RAM: 至少2GB，大文件(>1GB)建议8GB+
    - 存储: 输出目录需要输入文件2倍空间
    - CPU: 多核处理器可提升大文件处理速度
    
    性能指标:
    - 小文件(<10MB): 通常几秒内完成
    - 中等文件(10MB-100MB): 1-10分钟
    - 大文件(100MB-1GB): 10-60分钟
    - 超大文件(>1GB): 考虑分批处理
    
    内存使用:
    - 基础内存: ~100MB
    - 文件大小相关: ~2-3倍于输入文件大小
    - 多线程模式: 每线程额外~50MB
    
    故障排除 | Troubleshooting:
    
    常见问题:
    1. "文件不存在": 检查输入文件路径拼写
    2. "格式错误": 验证GFF文件格式规范性
    3. "内存不足": 减少线程数或处理较小文件块
    4. "权限拒绝": 检查输出目录写入权限
    5. "编码错误": 确保文件使用UTF-8编码
    
    格式要求:
    - 文件必须是制表符分隔
    - 第9列属性必须包含ID字段
    - 基因特征是处理的基础单位
    - 避免重复的ID值
    
    解决方案:
    - 使用--show-sample预览转换效果
    - 启用--verbose查看详细错误信息
    - 检查输入文件前几行格式
    - 验证物种前缀不包含特殊字符
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 预先验证GFF文件完整性
       - 💾 备份原始注释文件
       - 🗄️ 确保足够的磁盘空间
       - 🏷️ 使用描述性的物种标识符
    
    2️⃣ ⚙️ 参数选择:
       - 🎯 选择有意义的物种前缀
       - 🔢 起始编号避免与现有ID冲突
       - 📏 步长设置考虑未来扩展需求
       - 🧵 根据文件大小调整线程数
    
    3️⃣ 🔍 质量控制:
       - 👀 使用--show-sample预览结果
       - 📊 比较转换前后的基因数量
       - ✅ 验证关键基因的ID映射正确性
       - 📋 检查输出文件的GFF3格式合规性
    
    4️⃣ 📦 批量处理:
       - 📝 建立统一的命名规范
       - 🎯 为不同物种分配不同的编号范围
       - 📋 记录转换参数便于重现
       - ✅ 建立质量检查流程
    
    🔍 输出验证 | Output Validation:
    
    🤖 自动检查项目:
    - 📋 GFF3格式标准合规性
    - 🆔 基因ID唯一性验证
    - 🔗 父子关系完整性检查
    - 📍 坐标范围合理性验证
    - 🏷️ 属性字段完整性确认
    
    👥 手动验证建议:
    - 🎲 随机抽查转换后的基因ID格式
    - 🧬 确认关键基因的功能注释保持不变
    - 🔗 验证转录本-基因关系正确性
    - 🔤 检查特殊字符处理情况
    - 🛠️ 测试输出文件在下游工具中的兼容性
    
    📚 引用和参考 | Citation & References:
    
    📖 相关标准和规范:
    - 📋 GFF3 Format Specification: Sequence Ontology Project
    - 🧬 Gene Ontology: http://geneontology.org/
    - ⚖️ FAIR Data Principles: Wilkinson et al. (2016)
    """
    
    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    gff_main = _lazy_import_gff_main()
    
    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['gffconverter.py']
    
    # 必需参数 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-s', species_name])
    args.extend(['-p', species_prefix])
    
    # 可选参数（只有在非默认值时才添加）| Optional parameters (add only when non-default)
    if start_num != 10:
        args.extend(['--start-num', str(start_num)])
    
    if step != 10:
        args.extend(['--step', str(step)])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # 布尔选项 | Boolean options
    if verbose:
        args.append('--verbose')
    
    if keep_intermediate:
        args.append('--keep-intermediate')
    
    # 特殊参数 | Special parameters
    if show_sample is not None:
        args.extend(['--show-sample', str(show_sample)])
    
    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 | Call original main function
        gff_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 GFF格式转换被用户中断 | GFF format conversion interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 GFF格式转换失败 | GFF format conversion failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv