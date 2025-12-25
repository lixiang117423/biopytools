# """
# BWA比对分析命令 | BWA Alignment Analysis Command
# 高级优化版本：解决--help响应速度问题
# """

# import click
# import sys
# import os


# def _lazy_import_bwa_main():
#     """懒加载BWA main函数 | Lazy load BWA main function"""
#     try:
#         from ...bwa.main import main as bwa_main
#         return bwa_main
#     except ImportError as e:
#         click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
#         sys.exit(1)


# def _is_help_request():
#     """检查是否是帮助请求 | Check if this is a help request"""
#     help_flags = {'-h', '--help'}
#     return any(arg in help_flags for arg in sys.argv)


# def _validate_path_exists(path):
#     """验证路径是否存在（仅在非帮助模式下）| Validate path existence (only in non-help mode)"""
#     if not _is_help_request() and not os.path.exists(path):
#         raise click.BadParameter(f"❌ 路径不存在 | Path does not exist: {path}")
#     return path


# @click.command(short_help='🧬 BWA全基因组比对分析工具：高效的序列比对和覆盖度分析',
#                context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# @click.option('--genome', '-g',
#               required=True,
#               callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
#               help='🧬 参考基因组文件 | Reference genome file')
# @click.option('--input-dir', '-i',
#               required=True,
#               callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
#               help='📁 输入FASTQ文件目录 | Input FASTQ directory')
# @click.option('--pattern', '-p',
#               required=True,
#               type=str,
#               help='🔍 FASTQ文件匹配模式 (如 _1.clean.fq.gz) | FASTQ file pattern')
# @click.option('--output-dir', '-o',
#               default='./bwa_output',
#               type=click.Path(),
#               help='📂 输出目录 | Output directory (default: ./bwa_output)')
# @click.option('--threads', '-t',
#               type=int,
#               default=88,
#               help='🧵 线程数 | Number of threads (default: 88)')
# # BWA算法参数
# @click.option('--bwa-k',
#               type=int,
#               default=19,
#               help='🔬 最小种子长度 | Minimum seed length (default: 19)')
# @click.option('--bwa-w',
#               type=int,
#               default=100,
#               help='📏 带宽 | Band width (default: 100)')
# @click.option('--bwa-d',
#               type=int,
#               default=100,
#               help='📐 X-dropoff | Off-diagonal X-dropoff (default: 100)')
# @click.option('--bwa-r',
#               type=float,
#               default=1.5,
#               help='🎯 内部种子因子 | Internal seed factor (default: 1.5)')
# @click.option('--bwa-c',
#               type=int,
#               default=500,
#               help='🔢 种子出现次数阈值 | Seed occurrence threshold (default: 500)')
# @click.option('--bwa-D',
#               type=float,
#               default=0.50,
#               help='📉 短链丢弃比例 | Short chain drop fraction (default: 0.50)')
# @click.option('--bwa-W',
#               type=int,
#               default=0,
#               help='📏 最小链长 | Minimum chain length (default: 0)')
# @click.option('--bwa-m',
#               type=int,
#               default=50,
#               help='🔄 配对拯救轮数 | Mate rescue rounds (default: 50)')
# @click.option('--bwa-S',
#               is_flag=True,
#               help='⏭️  跳过配对拯救 | Skip mate rescue')
# @click.option('--bwa-P',
#               is_flag=True,
#               help='⏭️  跳过配对 | Skip pairing')
# # BWA打分参数
# @click.option('--bwa-A',
#               type=int,
#               default=1,
#               help='✅ 匹配得分 | Match score (default: 1)')
# @click.option('--bwa-B',
#               type=int,
#               default=4,
#               help='❌ 错配罚分 | Mismatch penalty (default: 4)')
# @click.option('--bwa-O',
#               default="6,6",
#               help='🔓 gap开放罚分 | Gap open penalty (default: 6,6)')
# @click.option('--bwa-E',
#               default="1,1",
#               help='➡️  gap延伸罚分 | Gap extension penalty (default: 1,1)')
# @click.option('--bwa-L',
#               default="5,5",
#               help='✂️  末端剪切罚分 | Clipping penalty (default: 5,5)')
# @click.option('--bwa-U',
#               type=int,
#               default=17,
#               help='💔 未配对罚分 | Unpaired penalty (default: 17)')
# # BWA输出参数
# @click.option('--bwa-M',
#               is_flag=True,
#               help='🏷️  标记次要比对 | Mark shorter split hits as secondary')
# @click.option('--bwa-T',
#               type=int,
#               default=30,
#               help='⭐ 最小输出得分 | Minimum score to output (default: 30)')
# @click.option('--bwa-a',
#               is_flag=True,
#               help='📋 输出所有比对 | Output all alignments')
# @click.option('--bwa-C',
#               is_flag=True,
#               help='💬 附加FASTQ注释 | Append FASTA/FASTQ comment')
# @click.option('--bwa-V',
#               is_flag=True,
#               help='📋 输出参考序列头 | Output reference FASTA header')
# @click.option('--bwa-Y',
#               is_flag=True,
#               help='✂️  软剪切补充比对 | Soft clipping for supplementary alignments')
# # 后处理参数
# @click.option('--markdup',
#               is_flag=True,
#               help='🔖 标记重复序列 | Mark duplicate reads')
# @click.option('--remove-dup',
#               is_flag=True,
#               help='🗑️  移除重复序列 | Remove duplicate reads')
# # 覆盖度参数
# @click.option('--min-base-quality',
#               type=int,
#               default=0,
#               help='⭐ 最小碱基质量 | Minimum base quality (default: 0)')
# @click.option('--min-mapping-quality',
#               type=int,
#               default=0,
#               help='🎯 最小比对质量 | Minimum mapping quality (default: 0)')
# @click.option('--max-depth',
#               type=int,
#               default=0,
#               help='📊 最大深度限制(0=无限制) | Max depth limit (0=no limit, default: 0)')
# # 滑窗参数
# @click.option('--window-size',
#               type=int,
#               default=1000000,
#               help='🪟  窗口大小(bp) | Window size in bp (default: 1000000)')
# @click.option('--step-size',
#               type=int,
#               default=100000,
#               help='👣 步长(bp) | Step size in bp (default: 100000)')
# # 其他参数
# @click.option('--resume',
#               is_flag=True,
#               help='🔄 断点续传，跳过已完成样品 | Resume, skip completed samples')
# @click.option('--keep-sam',
#               is_flag=True,
#               help='💾 保留SAM文件 | Keep SAM files')
# def bwa(genome, input_dir, pattern, output_dir, threads, bwa_k, bwa_w, bwa_d, bwa_r,
#              bwa_c, bwa_D, bwa_W, bwa_m, bwa_S, bwa_P, bwa_A, bwa_B, bwa_O, bwa_E, bwa_L, bwa_U,
#              bwa_M, bwa_T, bwa_a, bwa_C, bwa_V, bwa_Y, markdup, remove_dup, min_base_quality,
#              min_mapping_quality, max_depth, window_size, step_size, resume, keep_sam):
#     """
#     🧬 BWA全基因组比对分析工具 | BWA Whole Genome Alignment Tool
    
#     专业的BWA-MEM比对工具，提供从FASTQ到BAM的完整分析流程，
#     包括序列比对、重复序列标记、覆盖度分析和统计报告生成。
    
#     ✨ 功能特点 | Features:
#     - 🚀 高效的BWA-MEM算法比对
#     - 📊 全面的比对质量统计
#     - 📈 详细的覆盖度分析
#     - 🔖 智能重复序列处理
#     - 🪟 滑窗覆盖度计算
#     - 🔄 支持断点续传
#     - 📦 批量样品并行处理
    
#     💡 示例 | Examples:
    
#     \b
#     # 🎯 基本比对
#     biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz
    
#     \b
#     # 🔬 完整分析流程
#     biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
#         -o results --markdup -t 64
    
#     \b
#     # ⚙️ 自定义BWA参数
#     biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
#         --bwa-k 25 --bwa-M --bwa-t 88
    
#     \b
#     # 📊 高质量过滤分析
#     biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
#         --min-base-quality 20 --min-mapping-quality 30
    
#     📁 输入要求 | Input Requirements:
#     - 🧬 参考基因组FASTA文件
#     - 📁 配对的FASTQ文件目录
#     - 🔍 正确的文件匹配模式
    
#     📊 输出文件 | Output Files:
#     - 🗂️ BAM比对文件和索引
#     - 📈 覆盖度统计结果
#     - 📋 比对质量报告
#     - 📊 汇总统计信息
#     """
    
#     # 🚀 懒加载
#     bwa_main = _lazy_import_bwa_main()
    
#     # 构建参数列表
#     args = ['bwa_align.py']
#     args.extend(['-g', genome])
#     args.extend(['-i', input_dir])
#     args.extend(['-p', pattern])
    
#     if output_dir != './bwa_output':
#         args.extend(['-o', output_dir])
#     if threads != 88:
#         args.extend(['-t', str(threads)])
    
#     # BWA算法参数
#     if bwa_k != 19:
#         args.extend(['--bwa-k', str(bwa_k)])
#     if bwa_w != 100:
#         args.extend(['--bwa-w', str(bwa_w)])
#     if bwa_d != 100:
#         args.extend(['--bwa-d', str(bwa_d)])
#     if bwa_r != 1.5:
#         args.extend(['--bwa-r', str(bwa_r)])
#     if bwa_c != 500:
#         args.extend(['--bwa-c', str(bwa_c)])
#     if bwa_D != 0.50:
#         args.extend(['--bwa-D', str(bwa_D)])
#     if bwa_W != 0:
#         args.extend(['--bwa-W', str(bwa_W)])
#     if bwa_m != 50:
#         args.extend(['--bwa-m', str(bwa_m)])
#     if bwa_S:
#         args.append('--bwa-S')
#     if bwa_P:
#         args.append('--bwa-P')
    
#     # BWA打分参数
#     if bwa_A != 1:
#         args.extend(['--bwa-A', str(bwa_A)])
#     if bwa_B != 4:
#         args.extend(['--bwa-B', str(bwa_B)])
#     if bwa_O != "6,6":
#         args.extend(['--bwa-O', bwa_O])
#     if bwa_E != "1,1":
#         args.extend(['--bwa-E', bwa_E])
#     if bwa_L != "5,5":
#         args.extend(['--bwa-L', bwa_L])
#     if bwa_U != 17:
#         args.extend(['--bwa-U', str(bwa_U)])
    
#     # BWA输出参数
#     if bwa_M:
#         args.append('--bwa-M')
#     if bwa_T != 30:
#         args.extend(['--bwa-T', str(bwa_T)])
#     if bwa_a:
#         args.append('--bwa-a')
#     if bwa_C:
#         args.append('--bwa-C')
#     if bwa_V:
#         args.append('--bwa-V')
#     if bwa_Y:
#         args.append('--bwa-Y')
    
#     # 后处理参数
#     if markdup:
#         args.append('--markdup')
#     if remove_dup:
#         args.append('--remove-dup')
    
#     # 覆盖度参数
#     if min_base_quality != 0:
#         args.extend(['--min-base-quality', str(min_base_quality)])
#     if min_mapping_quality != 0:
#         args.extend(['--min-mapping-quality', str(min_mapping_quality)])
#     if max_depth != 0:
#         args.extend(['--max-depth', str(max_depth)])
    
#     # 滑窗参数
#     if window_size != 1000000:
#         args.extend(['--window-size', str(window_size)])
#     if step_size != 100000:
#         args.extend(['--step-size', str(step_size)])
    
#     # 其他参数
#     if resume:
#         args.append('--resume')
#     if keep_sam:
#         args.append('--keep-sam')
    
#     original_argv = sys.argv
#     sys.argv = args
    
#     try:
#         bwa_main()
#     except SystemExit as e:
#         if e.code != 0:
#             sys.exit(e.code)
#     except KeyboardInterrupt:
#         click.echo("\n🛑 BWA比对分析被用户中断", err=True)
#         sys.exit(1)
#     except Exception as e:
#         click.echo(f"💥 BWA比对分析失败: {e}", err=True)
#         sys.exit(1)
#     finally:
#         sys.argv = original_argv


"""
BWA比对分析命令 | BWA Alignment Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_bwa_main():
    """懒加载BWA main函数 | Lazy load BWA main function"""
    try:
        from ...bwa.main import main as bwa_main
        return bwa_main
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


@click.command(short_help='🧬 BWA全基因组比对分析工具：高效的序列比对和覆盖度分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='🧬 参考基因组文件 | Reference genome file')
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='📁 输入FASTQ文件目录 | Input FASTQ directory')
@click.option('--pattern', '-p',
              required=True,
              type=str,
              help='🔍 FASTQ文件匹配模式 (如 _1.clean.fq.gz) | FASTQ file pattern')
@click.option('--output-dir', '-o',
              default='./bwa_output',
              type=click.Path(),
              help='📂 输出目录 | Output directory (default: ./bwa_output)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
# BWA算法参数
@click.option('--bwa-k',
              type=int,
              default=19,
              help='🔬 最小种子长度 | Minimum seed length (default: 19)')
@click.option('--bwa-w',
              type=int,
              default=100,
              help='📏 带宽 | Band width (default: 100)')
@click.option('--bwa-d',
              type=int,
              default=100,
              help='📐 X-dropoff | Off-diagonal X-dropoff (default: 100)')
@click.option('--bwa-r',
              type=float,
              default=1.5,
              help='🎯 内部种子因子 | Internal seed factor (default: 1.5)')
@click.option('--bwa-c',
              type=int,
              default=500,
              help='🔢 种子出现次数阈值 | Seed occurrence threshold (default: 500)')
@click.option('--bwa-D', 'bwa_drop_ratio',  # 注意这里！指定参数名
              type=float,
              default=0.50,
              help='📉 短链丢弃比例 | Short chain drop fraction (default: 0.50)')
@click.option('--bwa-W', 'bwa_min_chain',  # 注意这里！指定参数名
              type=int,
              default=0,
              help='📏 最小链长 | Minimum chain length (default: 0)')
@click.option('--bwa-m',
              type=int,
              default=50,
              help='🔄 配对拯救轮数 | Mate rescue rounds (default: 50)')
@click.option('--bwa-S', 'bwa_s',  # 注意这里！小写
              is_flag=True,
              help='⏭️  跳过配对拯救 | Skip mate rescue')
@click.option('--bwa-P', 'bwa_p',  # 注意这里！小写
              is_flag=True,
              help='⏭️  跳过配对 | Skip pairing')
# BWA打分参数
@click.option('--bwa-A', 'bwa_score_match',  # 注意这里！避免和bwa_a冲突
              type=int,
              default=1,
              help='✅ 匹配得分 | Match score (default: 1)')
@click.option('--bwa-B', 'bwa_score_mismatch',  # 注意这里！描述性名称
              type=int,
              default=4,
              help='❌ 错配罚分 | Mismatch penalty (default: 4)')
@click.option('--bwa-O', 'bwa_gap_open',  # 注意这里！
              default="6,6",
              help='🔓 gap开放罚分 | Gap open penalty (default: 6,6)')
@click.option('--bwa-E', 'bwa_gap_ext',  # 注意这里！
              default="1,1",
              help='➡️  gap延伸罚分 | Gap extension penalty (default: 1,1)')
@click.option('--bwa-L', 'bwa_clip',  # 注意这里！
              default="5,5",
              help='✂️  末端剪切罚分 | Clipping penalty (default: 5,5)')
@click.option('--bwa-U', 'bwa_unpaired',  # 注意这里！
              type=int,
              default=17,
              help='💔 未配对罚分 | Unpaired penalty (default: 17)')
# BWA输出参数
@click.option('--bwa-M', 'bwa_mark_secondary',  # 注意这里！避免和bwa_m冲突
              is_flag=True,
              help='🏷️  标记次要比对 | Mark shorter split hits as secondary')
@click.option('--bwa-T', 'bwa_min_score',  # 注意这里！
              type=int,
              default=30,
              help='⭐ 最小输出得分 | Minimum score to output (default: 30)')
@click.option('--bwa-a', 'bwa_all_align',  # 注意这里！避免和bwa_A冲突
              is_flag=True,
              help='📋 输出所有比对 | Output all alignments')
@click.option('--bwa-C', 'bwa_append_comment',  # 注意这里！避免和bwa_c冲突
              is_flag=True,
              help='💬 附加FASTQ注释 | Append FASTA/FASTQ comment')
@click.option('--bwa-V', 'bwa_ref_header',  # 注意这里！
              is_flag=True,
              help='📋 输出参考序列头 | Output reference FASTA header')
@click.option('--bwa-Y', 'bwa_soft_clip',  # 注意这里！避免和bwa_y冲突
              is_flag=True,
              help='✂️  软剪切补充比对 | Soft clipping for supplementary alignments')
# 后处理参数
@click.option('--markdup',
              is_flag=True,
              help='🔖 标记重复序列 | Mark duplicate reads')
@click.option('--remove-dup',
              is_flag=True,
              help='🗑️  移除重复序列 | Remove duplicate reads')
# 覆盖度参数
@click.option('--min-base-quality',
              type=int,
              default=0,
              help='⭐ 最小碱基质量 | Minimum base quality (default: 0)')
@click.option('--min-mapping-quality',
              type=int,
              default=0,
              help='🎯 最小比对质量 | Minimum mapping quality (default: 0)')
@click.option('--max-depth',
              type=int,
              default=0,
              help='📊 最大深度限制(0=无限制) | Max depth limit (0=no limit, default: 0)')
# 滑窗参数
@click.option('--window-size',
              type=int,
              default=1000000,
              help='🪟  窗口大小(bp) | Window size in bp (default: 1000000)')
@click.option('--step-size',
              type=int,
              default=100000,
              help='👣 步长(bp) | Step size in bp (default: 100000)')
# 其他参数
@click.option('--resume',
              is_flag=True,
              help='🔄 断点续传，跳过已完成样品 | Resume, skip completed samples')
@click.option('--keep-sam',
              is_flag=True,
              help='💾 保留SAM文件 | Keep SAM files')
def bwa(genome, input_dir, pattern, output_dir, threads, 
        bwa_k, bwa_w, bwa_d, bwa_r, bwa_c, 
        bwa_drop_ratio, bwa_min_chain, bwa_m, bwa_s, bwa_p,  # 注意这里改了
        bwa_score_match, bwa_score_mismatch, bwa_gap_open, bwa_gap_ext, bwa_clip, bwa_unpaired,  # 注意这里改了
        bwa_mark_secondary, bwa_min_score, bwa_all_align, bwa_append_comment, bwa_ref_header, bwa_soft_clip,  # 注意这里改了
        markdup, remove_dup, min_base_quality, min_mapping_quality, max_depth, 
        window_size, step_size, resume, keep_sam):
    """
    🧬 BWA全基因组比对分析工具 | BWA Whole Genome Alignment Tool
    
    专业的BWA-MEM比对工具，提供从FASTQ到BAM的完整分析流程，
    包括序列比对、重复序列标记、覆盖度分析和统计报告生成。
    
    ✨ 功能特点 | Features:
    - 🚀 高效的BWA-MEM算法比对
    - 📊 全面的比对质量统计
    - 📈 详细的覆盖度分析
    - 🔖 智能重复序列处理
    - 🪟 滑窗覆盖度计算
    - 🔄 支持断点续传
    - 📦 批量样品并行处理
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本比对
    biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz
    
    \b
    # 🔬 完整分析流程
    biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
        -o results --markdup -t 64
    
    \b
    # ⚙️ 自定义BWA参数
    biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
        --bwa-k 25 --bwa-M --bwa-t 88
    
    \b
    # 📊 高质量过滤分析
    biopytools bwa -g genome.fa -i fastq_dir -p _1.clean.fq.gz \\
        --min-base-quality 20 --min-mapping-quality 30
    
    📁 输入要求 | Input Requirements:
    - 🧬 参考基因组FASTA文件
    - 📁 配对的FASTQ文件目录
    - 🔍 正确的文件匹配模式
    
    📊 输出文件 | Output Files:
    - 🗂️ BAM比对文件和索引
    - 📈 覆盖度统计结果
    - 📋 比对质量报告
    - 📊 汇总统计信息
    """
    
    # 🚀 懒加载
    bwa_main = _lazy_import_bwa_main()
    
    # 构建参数列表
    args = ['bwa_align.py']
    args.extend(['-g', genome])
    args.extend(['-i', input_dir])
    args.extend(['-p', pattern])
    
    if output_dir != './bwa_output':
        args.extend(['-o', output_dir])
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # BWA算法参数
    if bwa_k != 19:
        args.extend(['--bwa-k', str(bwa_k)])
    if bwa_w != 100:
        args.extend(['--bwa-w', str(bwa_w)])
    if bwa_d != 100:
        args.extend(['--bwa-d', str(bwa_d)])
    if bwa_r != 1.5:
        args.extend(['--bwa-r', str(bwa_r)])
    if bwa_c != 500:
        args.extend(['--bwa-c', str(bwa_c)])
    if bwa_drop_ratio != 0.50:  # 改了
        args.extend(['--bwa-D', str(bwa_drop_ratio)])
    if bwa_min_chain != 0:  # 改了
        args.extend(['--bwa-W', str(bwa_min_chain)])
    if bwa_m != 50:
        args.extend(['--bwa-m', str(bwa_m)])
    if bwa_s:  # 改了
        args.append('--bwa-S')
    if bwa_p:  # 改了
        args.append('--bwa-P')
    
    # BWA打分参数
    if bwa_score_match != 1:  # 改了
        args.extend(['--bwa-A', str(bwa_score_match)])
    if bwa_score_mismatch != 4:  # 改了
        args.extend(['--bwa-B', str(bwa_score_mismatch)])
    if bwa_gap_open != "6,6":  # 改了
        args.extend(['--bwa-O', bwa_gap_open])
    if bwa_gap_ext != "1,1":  # 改了
        args.extend(['--bwa-E', bwa_gap_ext])
    if bwa_clip != "5,5":  # 改了
        args.extend(['--bwa-L', bwa_clip])
    if bwa_unpaired != 17:  # 改了
        args.extend(['--bwa-U', str(bwa_unpaired)])
    
    # BWA输出参数
    if bwa_mark_secondary:  # 改了
        args.append('--bwa-M')
    if bwa_min_score != 30:  # 改了
        args.extend(['--bwa-T', str(bwa_min_score)])
    if bwa_all_align:  # 改了
        args.append('--bwa-a')
    if bwa_append_comment:  # 改了
        args.append('--bwa-C')
    if bwa_ref_header:  # 改了
        args.append('--bwa-V')
    if bwa_soft_clip:  # 改了
        args.append('--bwa-Y')
    
    # 后处理参数
    if markdup:
        args.append('--markdup')
    if remove_dup:
        args.append('--remove-dup')
    
    # 覆盖度参数
    if min_base_quality != 0:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if min_mapping_quality != 0:
        args.extend(['--min-mapping-quality', str(min_mapping_quality)])
    if max_depth != 0:
        args.extend(['--max-depth', str(max_depth)])
    
    # 滑窗参数
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if step_size != 100000:
        args.extend(['--step-size', str(step_size)])
    
    # 其他参数
    if resume:
        args.append('--resume')
    if keep_sam:
        args.append('--keep-sam')
    
    original_argv = sys.argv
    sys.argv = args
    
    try:
        bwa_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 BWA比对分析被用户中断", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 BWA比对分析失败: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv