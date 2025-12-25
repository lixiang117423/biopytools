"""
🧬 BWA-GATK流程分析命令 | BWA-GATK Pipeline Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_bwa_gatk_main():
    """懒加载BWA-GATK main函数 | Lazy load BWA-GATK main function"""
    try:
        from ...bwa_gatk.main import main as bwa_gatk_main
        return bwa_gatk_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path, path_type="file"):
    """验证路径是否存在（仅在非帮助模式下）| Validate path existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(
            f"❌ {path_type}不存在 | {path_type.capitalize()} does not exist: {path}"
        )
    return path


@click.command(
    short_help='🧬 BWA-GATK流程：高通量基因组变异检测分析流程',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value, "input") if value else None,
              help='📁 输入FASTQ文件或目录 | Input FASTQ file or directory')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value, "reference") if value else None,
              help='🧬 参考基因组FASTA文件 | Reference genome FASTA file')

# ===== 输出控制 | Output Control =====
@click.option('--output', '-o',
              default='./bwa_gatk_output',
              type=click.Path(),
              help='📂 输出目录 | Output directory (default: ./bwa_gatk_output)')
@click.option('--sample-name',
              help='🏷️  样本名称 | Sample name (auto-detect if not specified)')

# ===== 性能参数 | Performance Parameters =====
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--memory',
              default='200G',
              help='💾 最大内存使用 | Maximum memory usage (default: 200G)')

# ===== BWA比对参数 | BWA Alignment Parameters =====
@click.option('--bwa-algorithm',
              type=click.Choice(['mem', 'aln']),
              default='mem',
              help='⚙️  BWA比对算法 | BWA alignment algorithm (default: mem)')
@click.option('--min-seed-length',
              type=int,
              default=19,
              help='🌱 最小种子长度 | Minimum seed length (default: 19)')
@click.option('--band-width',
              type=int,
              default=100,
              help='📏 带宽参数 | Band width (default: 100)')

# ===== GATK变异检测参数 | GATK Variant Calling Parameters =====
@click.option('--min-base-quality',
              type=int,
              default=20,
              help='🎯 最小碱基质量值 | Minimum base quality score (default: 20)')
@click.option('--min-mapping-quality',
              type=int,
              default=20,
              help='🗺️  最小比对质量值 | Minimum mapping quality score (default: 20)')
@click.option('--stand-call-conf',
              type=float,
              default=30.0,
              help='📞 变异检出置信度阈值 | Stand call confidence threshold (default: 30.0)')

# ===== 变异过滤参数 | Variant Filtering Parameters =====
@click.option('--filter-expression',
              help='🔍 变异过滤表达式 | Variant filter expression')
@click.option('--min-depth',
              type=int,
              default=10,
              help='📊 最小测序深度 | Minimum sequencing depth (default: 10)')
@click.option('--max-depth',
              type=int,
              help='📈 最大测序深度 | Maximum sequencing depth')

# ===== 流程控制 | Pipeline Control =====
@click.option('--skip-bwa',
              is_flag=True,
              help='⏭️  跳过BWA比对步骤 | Skip BWA alignment step')
@click.option('--skip-gatk',
              is_flag=True,
              help='⏭️  跳过GATK变异检测步骤 | Skip GATK variant calling step')
@click.option('--only-align',
              is_flag=True,
              help='🎯 仅执行比对步骤 | Only perform alignment step')
@click.option('--only-call',
              is_flag=True,
              help='🎯 仅执行变异检测步骤 | Only perform variant calling step')
@click.option('--resume',
              is_flag=True,
              default=True,
              help='🔄 启用断点续传 | Enable resume from checkpoint (default: enabled)')
@click.option('--force',
              is_flag=True,
              help='🔥 强制重新运行 | Force rerun all steps')

# ===== 高级选项 | Advanced Options =====
@click.option('--known-sites',
              type=click.Path(exists=True),
              help='📋 已知变异位点VCF文件 | Known variant sites VCF file')
@click.option('--bqsr',
              is_flag=True,
              help='🔧 启用碱基质量分数重校正 | Enable Base Quality Score Recalibration')
@click.option('--remove-duplicates',
              is_flag=True,
              default=True,
              help='🗑️  移除PCR重复 | Remove PCR duplicates (default: enabled)')
@click.option('--emit-ref-confidence',
              type=click.Choice(['NONE', 'BP_RESOLUTION', 'GVCF']),
              default='NONE',
              help='📤 参考位点置信度输出模式 | Reference confidence mode (default: NONE)')

# ===== 工具路径 | Tool Paths =====
@click.option('--bwa-path',
              help='🔧 BWA工具路径 | BWA tool path')
@click.option('--gatk-path',
              help='🔧 GATK工具路径 | GATK tool path')
@click.option('--samtools-path',
              help='🔧 Samtools工具路径 | Samtools tool path')

# ===== 其他选项 | Other Options =====
@click.option('--verbose', '-v',
              is_flag=True,
              help='📢 详细输出模式 | Verbose output mode')
@click.option('--quiet', '-q',
              is_flag=True,
              help='🔇 静默模式 | Quiet mode')
@click.option('--keep-intermediate',
              is_flag=True,
              help='💾 保留中间文件 | Keep intermediate files')
def bwa_gatk(input, reference, output, sample_name, threads, memory,
             bwa_algorithm, min_seed_length, band_width,
             min_base_quality, min_mapping_quality, stand_call_conf,
             filter_expression, min_depth, max_depth,
             skip_bwa, skip_gatk, only_align, only_call, resume, force,
             known_sites, bqsr, remove_duplicates, emit_ref_confidence,
             bwa_path, gatk_path, samtools_path,
             verbose, quiet, keep_intermediate):
    """
    🧬 BWA-GATK变异检测流程 | BWA-GATK Variant Calling Pipeline
    
    整合BWA比对和GATK变异检测的完整基因组分析流程，支持从原始测序数据
    到高质量变异检测的自动化处理。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🔬 完整的分析流程:
       • BWA比对: 将测序数据比对到参考基因组
       • 质量控制: 碱基质量重校正(BQSR)
       • 去重处理: 移除PCR重复序列
       • 变异检测: GATK HaplotypeCaller
       • 变异过滤: 多维度质量过滤
    
    \b
    🚀 性能优化:
       • 多线程并行处理
       • 内存使用优化
       • 断点续传支持
       • 智能缓存机制
    
    \b
    🎯 高质量结果:
       • 准确的变异检测
       • 严格的质量控制
       • 详细的统计报告
       • 标准VCF输出
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    FASTQ文件格式:
       • 支持单端或双端测序数据
       • 文件命名: *_1.fastq.gz, *_2.fastq.gz
       • 支持压缩格式(.gz, .bz2)
       • 质量编码: Phred33或Phred64
    
    \b
    参考基因组要求:
       • FASTA格式(.fa, .fasta)
       • 需要预先建立索引
       • 支持自动索引生成
    
    📊 输出文件说明 | Output Files:
    
    \b
    主要输出文件:
       • *.bam - 比对结果文件
       • *.vcf.gz - 变异检测结果
       • *.html - 质量评估报告
       • *_stats.txt - 统计信息
    
    \b
    中间文件(可选保留):
       • *.sorted.bam - 排序后的BAM
       • *.dedup.bam - 去重后的BAM
       • *.recal.table - BQSR校正表
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本用法 - 完整流程
    biopytools bwa-gatk -i reads.fastq.gz -r ref.fasta -o output/
    
    \b
    # 🔬 双端测序数据分析
    biopytools bwa-gatk -i fastq_dir/ -r ref.fa -o results/ \\
        --sample-name sample01 -t 32
    
    \b
    # 🎯 仅执行比对
    biopytools bwa-gatk -i reads.fq.gz -r ref.fa -o align_out/ \\
        --only-align --threads 16
    
    \b
    # 🔧 启用BQSR质量校正
    biopytools bwa-gatk -i data.fq -r ref.fasta -o output/ \\
        --known-sites dbsnp.vcf --bqsr
    
    \b
    # ⚙️ 自定义变异过滤
    biopytools bwa-gatk -i reads.fastq -r genome.fa -o vcf_out/ \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \\
        --min-depth 20 --stand-call-conf 50
    
    \b
    # 🚀 高性能配置
    biopytools bwa-gatk -i input/ -r ref.fa -o output/ \\
        -t 64 --memory 400G --resume
    
    \b
    # 🧬 生成GVCF格式
    biopytools bwa-gatk -i sample.fq -r ref.fa -o gvcf_out/ \\
        --emit-ref-confidence GVCF
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 全基因组重测序(WGS)分析
    • 外显子组测序(WES)分析
    • 目标区域测序变异检测
    • 群体遗传学研究
    • 临床基因检测
    • 育种辅助选择
    
    ⚡ 性能要求 | Performance Requirements:
    
    \b
    硬件建议:
       • CPU: 16核以上推荐
       • 内存: 32GB以上(大基因组需64GB+)
       • 存储: 输入数据量的5-10倍空间
       • 临时空间: /tmp需要足够空间
    
    \b
    软件依赖:
       • BWA 0.7.17+
       • GATK 4.0+
       • Samtools 1.10+
       • Picard Tools
       • Python 3.7+
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1. 内存不足 → 减少--threads或增加--memory
       2. 磁盘空间不足 → 清理中间文件或扩展存储
       3. 参考基因组索引缺失 → 自动生成或手动建立
       4. GATK内存错误 → 调整--memory参数
    
    \b
    优化建议:
       • 使用SSD存储提升I/O性能
       • 合理设置线程数(不超过物理核心数)
       • 已知变异位点可提升BQSR效果
       • 大样本分析建议使用GVCF模式
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据质量检查
       • 运行前检查FASTQ质量(FastQC)
       • 确认参考基因组完整性
       • 验证已知变异位点VCF格式
    
    \b
    2️⃣ 参数优化
       • WGS建议: --min-depth 10 --stand-call-conf 30
       • WES建议: --min-depth 20 --stand-call-conf 50
       • 低覆盖度: 降低阈值但增加后续过滤
    
    \b
    3️⃣ 结果验证
       • 检查BAM文件比对率和覆盖度
       • 评估变异检测的Ti/Tv比值
       • 对比已知变异位点验证准确性
    
    \b
    4️⃣ 批量处理
       • 使用循环处理多个样本
       • 考虑使用任务调度系统
       • 启用--resume避免重复计算
    
    📚 相关资源 | Related Resources:
    
    \b
    • GATK最佳实践: https://gatk.broadinstitute.org/
    • BWA手册: http://bio-bwa.sourceforge.net/
    • VCF格式规范: https://samtools.github.io/hts-specs/
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 确保有足够的磁盘空间(输入数据的5-10倍)
    • 大基因组分析建议使用高内存节点
    • 变异过滤参数需根据具体项目调整
    • BQSR需要高质量的已知变异位点
    • 临床应用需要严格的质量控制
    """
    
    # 🚀 懒加载
    bwa_gatk_main = _lazy_import_bwa_gatk_main()
    
    # 构建参数列表
    args = ['bwa_gatk.py']
    args.extend(['-i', input])
    args.extend(['-r', reference])
    args.extend(['-o', output])
    
    # 样本信息
    if sample_name:
        args.extend(['--sample-name', sample_name])
    
    # 性能参数
    if threads != 88:
        args.extend(['-t', str(threads)])
    if memory != '200G':
        args.extend(['--memory', memory])
    
    # BWA参数
    if bwa_algorithm != 'mem':
        args.extend(['--bwa-algorithm', bwa_algorithm])
    if min_seed_length != 19:
        args.extend(['--min-seed-length', str(min_seed_length)])
    if band_width != 100:
        args.extend(['--band-width', str(band_width)])
    
    # GATK参数
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if min_mapping_quality != 20:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])
    if stand_call_conf != 30.0:
        args.extend(['--stand-call-conf', str(stand_call_conf)])
    
    # 过滤参数
    if filter_expression:
        args.extend(['--filter-expression', filter_expression])
    if min_depth != 10:
        args.extend(['--min-depth', str(min_depth)])
    if max_depth:
        args.extend(['--max-depth', str(max_depth)])
    
    # 流程控制
    if skip_bwa:
        args.append('--skip-bwa')
    if skip_gatk:
        args.append('--skip-gatk')
    if only_align:
        args.append('--only-align')
    if only_call:
        args.append('--only-call')
    if force:
        args.append('--force')
    if not resume:  # resume默认为True
        args.append('--no-resume')
    
    # 高级选项
    if known_sites:
        args.extend(['--known-sites', known_sites])
    if bqsr:
        args.append('--bqsr')
    if not remove_duplicates:  # remove_duplicates默认为True
        args.append('--no-remove-duplicates')
    if emit_ref_confidence != 'NONE':
        args.extend(['--emit-ref-confidence', emit_ref_confidence])
    
    # 工具路径
    if bwa_path:
        args.extend(['--bwa-path', bwa_path])
    if gatk_path:
        args.extend(['--gatk-path', gatk_path])
    if samtools_path:
        args.extend(['--samtools-path', samtools_path])
    
    # 其他选项
    if verbose:
        args.append('-v')
    if quiet:
        args.append('-q')
    if keep_intermediate:
        args.append('--keep-intermediate')
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        bwa_gatk_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 BWA-GATK流程被用户中断 | Pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 BWA-GATK流程执行失败 | Pipeline execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv