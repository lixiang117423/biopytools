"""
BUSCO质量评估分析命令 | BUSCO Quality Assessment Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_busco_main():
    """懒加载BUSCO main函数 | Lazy load BUSCO main function"""
    try:
        from ...busco.main import main as busco_main
        return busco_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径是否存在（仅在非帮助模式下）| Validate path existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(f"❌ 路径不存在 | Path does not exist: {path}")
    return path


@click.command(short_help='🧬 BUSCO质量评估分析工具：基因组/转录组/蛋白质组完整性和质量评估',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='🔤 输入文件或目录路径 | Input file or directory path')
@click.option('--lineage', '-l',
              required=True,
              type=str,
              help='🗄️  BUSCO数据库/谱系名称 | BUSCO database/lineage name')
@click.option('--output', '-o',
              default='./busco_output',
              type=click.Path(),
              help='📂 输出目录 | Output directory (default: ./busco_output)')
@click.option('--mode', '-m',
              type=click.Choice(['genome', 'geno', 'transcriptome', 'tran', 'proteins', 'prot']),
              default='genome',
              help='🔬 BUSCO分析模式 | BUSCO analysis mode (default: genome)')
@click.option('--cpu', '-c',
              type=int,
              default=88,
              help='⚡ CPU线程数 | Number of CPU threads (default: 88)')
@click.option('--sample-suffix',
              default='*.fa',
              help='🏷️  样本名称提取后缀模式 | Sample name extraction suffix pattern (default: *.fa)')
@click.option('--output-format',
              type=click.Choice(['txt', 'csv', 'xlsx']),
              default='txt',
              help='📄 输出文件格式 | Output file format (default: txt)')
@click.option('--force', '-f',
              is_flag=True,
              help='💪 强制重写现有文件 | Force rewriting of existing files')
@click.option('--augustus',
              is_flag=True,
              help='🔮 使用Augustus基因预测器 | Use Augustus gene predictor')
@click.option('--augustus-parameters',
              type=str,
              help='🎛️  Augustus附加参数 | Additional Augustus parameters')
@click.option('--augustus-species',
              type=str,
              help='🐛 Augustus物种名称 | Augustus species name')
@click.option('--auto-lineage',
              is_flag=True,
              help='🤖 自动选择谱系 | Automatically select lineage')
@click.option('--auto-lineage-euk',
              is_flag=True,
              help='🌿 自动选择真核生物谱系 | Automatically select eukaryote lineage')
@click.option('--auto-lineage-prok',
              is_flag=True,
              help='🦠 自动选择原核生物谱系 | Automatically select prokaryote lineage')
@click.option('--contig-break',
              type=int,
              default=10,
              help='🧩 Contig断点N数量 | Number of Ns for contig break (default: 10)')
@click.option('--datasets-version',
              default='odb12',
              help='📊 数据集版本 | Dataset version (default: odb12)')
@click.option('--download-path',
              type=str,
              help='💾 数据集下载路径 | Dataset download path')
@click.option('--evalue', '-e',
              type=float,
              default=1e-3,
              help='🎯 BLAST E值阈值 | BLAST E-value threshold (default: 1e-3)')
@click.option('--limit',
              type=int,
              default=3,
              help='🎯 候选区域数量限制 | Candidate region limit (default: 3)')
@click.option('--long',
              is_flag=True,
              help='⏱️  启用Augustus长模式优化 | Enable Augustus long mode optimization')
@click.option('--metaeuk',
              is_flag=True,
              help='🧠 使用Metaeuk基因预测器 | Use Metaeuk gene predictor')
@click.option('--metaeuk-parameters',
              type=str,
              help='🎛️  Metaeuk附加参数 | Additional Metaeuk parameters')
@click.option('--metaeuk-rerun-parameters',
              type=str,
              help='🔄 Metaeuk重运行参数 | Metaeuk rerun parameters')
@click.option('--miniprot',
              is_flag=True,
              help='🧬 使用Miniprot基因预测器 | Use Miniprot gene predictor')
@click.option('--skip-bbtools',
              is_flag=True,
              help='⏩ 跳过BBTools统计 | Skip BBTools statistics')
@click.option('--offline',
              is_flag=True,
              help='📡 离线模式 | Offline mode')
@click.option('--restart', '-r',
              is_flag=True,
              help='🔄 重启未完成的分析 | Restart incomplete analysis')
@click.option('--quiet', '-q',
              is_flag=True,
              help='🔇 静默模式 | Quiet mode')
@click.option('--scaffold-composition',
              is_flag=True,
              help='📊 生成scaffold组成文件 | Generate scaffold composition file')
@click.option('--tar',
              is_flag=True,
              help='📦 压缩子目录 | Compress subdirectories')
@click.option('--busco-path',
              default='busco',
              help='🔧 BUSCO软件路径 | BUSCO software path (default: busco)')
def busco(input, lineage, output, mode, cpu, sample_suffix, output_format, 
                  force, augustus, augustus_parameters, augustus_species, auto_lineage,
                  auto_lineage_euk, auto_lineage_prok, contig_break, datasets_version,
                  download_path, evalue, limit, long, metaeuk, metaeuk_parameters,
                  metaeuk_rerun_parameters, miniprot, skip_bbtools, offline, restart,
                  quiet, scaffold_composition, tar, busco_path):
    """
    🧬 BUSCO质量评估分析工具 | BUSCO Quality Assessment Analysis Tool
    
    一个专业的基因组、转录组和蛋白质组完整性评估工具，基于BUSCO (Benchmarking 
    Universal Single-Copy Orthologs) 方法，通过检测保守单拷贝直系同源基因来
    评估生物序列数据的完整性和质量，广泛应用于基因组学和转录组学研究。
    
    ✨ 功能特点 | Features:
    - 🎯 精确的基因组/转录组完整性评估
    - 🗄️ 支持多种BUSCO数据库谱系
    - 🔬 多种分析模式适配不同数据类型
    - 📊 详细的统计报告和可视化
    - ⚡ 批量样本并行处理优化
    - 🤖 智能谱系自动选择功能
    - 📦 多种输出格式支持
    
    🔄 分析流程 | Analysis Pipeline:
    1. 📖 解析输入序列文件或目录
    2. 🗄️ 加载指定的BUSCO数据库谱系
    3. 🔍 执行同源基因搜索和识别
    4. 📊 计算完整性统计指标
    5. 📝 生成详细的质量评估报告
    
    📊 BUSCO评估指标 | BUSCO Assessment Metrics:
    - 🟢 Complete (C): 完整基因数量和比例
    - 🔵 Complete Single-copy (S): 单拷贝完整基因
    - 🟡 Complete Duplicated (D): 重复拷贝完整基因
    - 🟠 Fragmented (F): 片段化基因数量
    - 🔴 Missing (M): 缺失基因数量
    
    🎯 应用场景 | Use Cases:
    - 🧬 基因组组装质量评估
    - 📊 转录组组装完整性检验
    - 🧪 蛋白质组预测准确性验证
    - 🔬 不同组装方法效果比较
    - 📈 测序深度充分性评估
    - 🌳 系统发育分析数据预处理
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本基因组分析
    biopytools busco -i genome.fa -l eukaryota_odb12
    
    \b
    # 🔬 转录组模式分析
    biopytools busco -i transcripts.fa -l fungi_odb12 \\
        -m transcriptome -o busco_transcriptome
    
    \b
    # 🧪 蛋白质组批量分析
    biopytools busco -i proteins_dir/ -l vertebrata_odb12 \\
        -m proteins --sample-suffix "*.pep.fa" -c 64
    
    \b
    # 🤖 自动谱系选择
    biopytools busco -i unknown_genome.fa \\
        --auto-lineage-euk -o auto_busco
    
    \b
    # 🔮 使用Augustus基因预测
    biopytools busco -i draft_genome.fa -l metazoa_odb12 \\
        --augustus --augustus-species human --long
    
    \b
    # 📊 多格式输出和压缩
    biopytools busco -i genomes_dir/ -l bacteria_odb12 \\
        --output-format xlsx --tar --force
    
    \b
    # 🧠 使用Metaeuk进行宏基因组分析
    biopytools busco -i metagenome.fa -l bacteria_odb12 \\
        --metaeuk --evalue 1e-5 --limit 5
    
    📁 输入数据要求 | Input Data Requirements:
    
    🧬 基因组模式 (genome):
    - 📄 FASTA格式的基因组序列文件
    - 🧩 支持完整基因组或草图基因组
    - 📦 可以是单个文件或包含多个文件的目录
    - 🔤 建议使用标准的染色体/contig命名
    
    📊 转录组模式 (transcriptome):
    - 📄 FASTA格式的转录本序列文件
    - 🧬 de novo组装或参考基因组组装均可
    - 📏 建议去除过短的转录本(<200bp)
    - 🔍 支持Trinity、SPAdes等组装工具输出
    
    🧪 蛋白质组模式 (proteins):
    - 📄 FASTA格式的蛋白质序列文件
    - 🔬 基因预测或实验验证的蛋白质序列
    - ⭐ 要求完整的开放阅读框
    - 🏷️ 序列标识符应唯一且有意义
    
    📝 文件格式示例:
    >chromosome_1
    ATCGATCGATCGATCG...
    >chromosome_2
    GCTAGCTAGCTAGCTA...
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    📋 必需参数:
    --input: 📁 输入文件或目录路径
        - 单文件：直接指定FASTA文件路径
        - 目录：批量处理目录内所有匹配文件
    --lineage: 🗄️ BUSCO数据库谱系名称
        - 常用谱系：eukaryota_odb12, bacteria_odb12, archaea_odb12
        - 专门谱系：vertebrata_odb12, fungi_odb12, plantae_odb12
        - 可使用 busco --list-datasets 查看所有可用谱系
    
    📂 输出控制:
    --output: 📂 结果输出目录
    --output-format: 📄 输出格式选择
        - txt: 制表符分隔文本文件 (默认)
        - csv: 逗号分隔文件
        - xlsx: Excel电子表格格式
    
    🔬 分析模式:
    --mode: 分析类型选择
        - genome/geno: 基因组完整性评估
        - transcriptome/tran: 转录组完整性评估
        - proteins/prot: 蛋白质组完整性评估
    
    ⚡ 性能参数:
    --cpu: 🚀 并行处理线程数
        - 建议设置为CPU核心数
        - 影响搜索和比对速度
    
    🤖 智能功能:
    --auto-lineage: 🤖 自动选择最佳谱系
    --auto-lineage-euk: 🌿 限定真核生物谱系自动选择
    --auto-lineage-prok: 🦠 限定原核生物谱系自动选择
    
    🔮 基因预测器选项:
    --augustus: 使用Augustus进行基因预测
        - 适用于真核生物基因组
        - 可指定物种模型提高准确性
    --metaeuk: 🧠 使用Metaeuk进行基因预测
        - 适用于宏基因组和环境样本
        - 对基因组碎片化容忍度高
    --miniprot: 🧬 使用Miniprot进行基因预测
        - 基于蛋白质同源性的快速预测
        - 适用于远缘物种分析
    
    🎯 搜索参数:
    --evalue: E值阈值控制搜索敏感性
    --limit: 🎯 每个BUSCO基因的候选区域数量限制
    
    🗄️ 数据管理:
    --datasets-version: 📊 BUSCO数据库版本
    --download-path: 💾 数据库下载存储路径
    --offline: 📡 离线模式，不检查数据库更新
    
    🔧 流程控制:
    --force: 💪 覆盖已存在的输出文件
    --restart: 🔄 从中断点继续未完成的分析
    --quiet: 🔇 减少日志输出的静默模式
    
    📦 输出增强:
    --scaffold-composition: 📊 生成scaffold组成统计
    --tar: 📦 自动压缩输出子目录
    
    📁 输出文件说明 | Output Files:
    
    📝 主要结果文件:
    - 📊 busco_summary.txt: 汇总统计表格
    - 📋 busco_results.{txt|csv|xlsx}: 详细结果数据
    - 📈 busco_report.txt: 分析过程报告
    - 📦 individual_results/: 各样本详细结果
    
    📊 BUSCO原始输出:
    - 🗂️ run_{sample}/: 每个样本的完整BUSCO输出
    - 📋 short_summary*.txt: BUSCO标准汇总
    - 📄 full_table*.tsv: 详细基因检测表
    - 📊 busco_sequences/: 检测到的序列文件
    
    📈 统计可视化 (如果支持):
    - 📊 busco_plot.png: 完整性统计柱状图
    - 📈 completeness_comparison.png: 样本间比较图
    - 📉 missing_genes_analysis.png: 缺失基因分析
    
    ⚡ 性能和系统要求 | Performance & System Requirements:
    
    🔧 依赖软件:
    - 🧰 BUSCO (v5.0+): 核心评估程序
    - 🔍 BLAST+ 或 DIAMOND: 序列搜索引擎
    - 🔮 Augustus: 真核基因预测 (可选)
    - 🧠 MetaEuk: 宏基因组基因预测 (可选)
    - 🧬 Miniprot: 蛋白质比对基因预测 (可选)
    - 🛠️ BBTools: 序列统计工具 (可选)
    
    💻 系统建议:
    - 💾 RAM: 至少8GB，大基因组建议32GB+
    - 🗄️ 存储: 至少输入文件大小的5倍空间
    - 🚀 CPU: 多核处理器显著提升速度
    - 🌐 网络: 首次运行需要下载数据库
    
    📊 数据规模处理能力:
    - 📄 小基因组(<100MB): 几分钟到数小时
    - 📋 中等基因组(100MB-1GB): 数小时
    - 📁 大基因组(1GB-10GB): 几小时到一天
    - 🗂️ 批量处理: 支持数百个样本并行
    
    💾 内存使用估算:
    - 🔵 基础内存: ~2GB
    - 📈 数据相关: 基因组大小的2-5倍
    - 🧵 线程相关: 每线程额外~500MB
    
    🛠️ 故障排除 | Troubleshooting:
    
    ⚠️ 常见问题:
    1. ❌ "BUSCO not found": 检查BUSCO安装和PATH设置
    2. 🗄️ "Dataset download failed": 检查网络连接或使用离线模式
    3. 💾 "内存不足": 减少线程数或增加系统内存
    4. 📁 "输出目录权限": 确保有写入权限
    5. 🔧 "Augustus error": 检查Augustus配置和物种模型
    
    📋 数据质量检查:
    - 🔍 验证输入文件FASTA格式正确性
    - 📏 检查序列长度分布合理性
    - 🏷️ 确认序列标识符格式规范
    - 🧬 评估输入数据的生物学意义
    
    💡 优化建议:
    - 🎯 根据物种选择最接近的谱系数据库
    - 🚀 合理设置线程数避免资源竞争
    - 📊 对比多种参数设置的结果
    - 🔄 使用restart功能避免重复计算
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 使用高质量的输入序列
       - 🔍 预先验证文件格式和完整性
       - 📏 移除过短或低质量序列
       - 🏷️ 使用有意义的序列标识符
    
    2️⃣ ⚙️ 参数选择:
       - 🎯 选择最接近的谱系数据库
       - 🔬 根据数据类型选择合适模式
       - 🤖 不确定物种时使用自动谱系选择
       - 🔮 复杂基因组考虑使用Augustus
    
    3️⃣ 🔍 结果解释:
       - 📊 关注Complete基因的比例
       - 🔴 分析Missing基因的模式
       - 🟡 评估Duplicated基因的意义
       - 📈 比较同类研究的BUSCO分数
    
    4️⃣ 📦 质量控制:
       - 🎲 随机验证部分基因的检测结果
       - 🔄 使用不同参数验证结果稳定性
       - 📚 参考相关物种的BUSCO基准
       - 🧬 结合其他质量评估方法
    
    🔍 结果解释指南 | Result Interpretation Guide:
    
    📊 BUSCO分数标准:
    - 🌟 优秀 (>95%): 高质量基因组/转录组
    - ✅ 良好 (90-95%): 可接受的质量水平
    - ⚠️ 一般 (80-90%): 可能需要进一步优化
    - ❌ 较差 (<80%): 建议重新组装或预测
    
    🎯 不同模式期望值:
    - 🧬 基因组模式: 通常>90%完整度
    - 📊 转录组模式: 因表达差异可能较低
    - 🧪 蛋白质组模式: 取决于预测方法质量
    
    📚 引用和参考 | Citation & References:
    
    📖 核心文献:
    - 🧰 BUSCO: Manni et al. (2021) Mol Biol Evol
    - 🔮 Augustus: Stanke et al. (2006) Nucleic Acids Res
    - 🧠 MetaEuk: Levy Karin et al. (2020) Microbiome
    - 🧬 Miniprot: Li (2023) Bioinformatics
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    busco_main = _lazy_import_busco_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['busco_analysis.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input])
    args.extend(['-l', lineage])
    
    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    if output != './busco_output':
        args.extend(['-o', output])
    
    if mode != 'genome':
        args.extend(['-m', mode])
    
    if cpu != 88:
        args.extend(['-c', str(cpu)])
    
    if sample_suffix != '*.fa':
        args.extend(['--sample-suffix', sample_suffix])
    
    if output_format != 'txt':
        args.extend(['--output-format', output_format])
    
    if contig_break != 10:
        args.extend(['--contig-break', str(contig_break)])
    
    if datasets_version != 'odb12':
        args.extend(['--datasets-version', datasets_version])
    
    if evalue != 1e-3:
        args.extend(['-e', str(evalue)])
    
    if limit != 3:
        args.extend(['--limit', str(limit)])
    
    if busco_path != 'busco':
        args.extend(['--busco-path', busco_path])
    
    # 字符串参数 🔤 | String parameters
    if augustus_parameters:
        args.extend(['--augustus-parameters', augustus_parameters])
    
    if augustus_species:
        args.extend(['--augustus-species', augustus_species])
    
    if download_path:
        args.extend(['--download-path', download_path])
    
    if metaeuk_parameters:
        args.extend(['--metaeuk-parameters', metaeuk_parameters])
    
    if metaeuk_rerun_parameters:
        args.extend(['--metaeuk-rerun-parameters', metaeuk_rerun_parameters])
    
    # 布尔选项 🚩 | Boolean options
    if force:
        args.append('-f')
    
    if augustus:
        args.append('--augustus')
    
    if auto_lineage:
        args.append('--auto-lineage')
    
    if auto_lineage_euk:
        args.append('--auto-lineage-euk')
    
    if auto_lineage_prok:
        args.append('--auto-lineage-prok')
    
    if long:
        args.append('--long')
    
    if metaeuk:
        args.append('--metaeuk')
    
    if miniprot:
        args.append('--miniprot')
    
    if skip_bbtools:
        args.append('--skip-bbtools')
    
    if offline:
        args.append('--offline')
    
    if restart:
        args.append('-r')
    
    if quiet:
        args.append('-q')
    
    if scaffold_composition:
        args.append('--scaffold-composition')
    
    if tar:
        args.append('--tar')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        busco_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 BUSCO质量评估分析被用户中断 | BUSCO quality assessment analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 BUSCO质量评估分析失败 | BUSCO quality assessment analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv