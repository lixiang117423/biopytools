# """
# 🚀 Parabricks WGS分析命令 | Parabricks WGS Analysis Command
# 高级优化版本：支持GVCF输出和Joint Calling联合变异检测
# """

# import click
# import sys
# import os


# def _lazy_import_parabricks_main():
#     """懒加载Parabricks main函数 | Lazy load Parabricks main function"""
#     try:
#         from ...parabricks.main import main as parabricks_main
#         return parabricks_main
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


# @click.command(
#     short_help='🚀 基于GPU加速的全基因组测序(WGS)批处理分析工具（支持GVCF和Joint Calling）',
#     context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
# )
# # ===== 必需参数 | Required Parameters =====
# @click.option('--input-dir', '-i',
#               required=True,
#               callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
#               help='📁 输入目录路径 | Input directory path (containing clean FASTQ files)')
# @click.option('--output-dir', '-o',
#               required=True,
#               type=click.Path(),
#               help='📂 输出目录路径 | Output directory path')
# @click.option('--reference', '-r',
#               required=True,
#               callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
#               help='🧬 参考基因组文件路径 | Reference genome file path')

# # ===== 性能参数 | Performance Parameters =====
# @click.option('--threads', '-t',
#               default=88,
#               type=int,
#               help='🧵 线程数 | Number of threads (default: 88)')
# @click.option('--parabricks-path',
#               default='/share/apps/containers/parabricks.sif',
#               type=str,
#               help='🔧 Parabricks程序路径 | Parabricks program path')
# @click.option('--tmp-dir',
#               type=click.Path(),
#               help='💾 临时目录路径 | Temporary directory path (default: tmp under output directory)')

# # ===== GVCF和Joint Calling参数 | GVCF and Joint Calling Parameters =====
# # ===== GVCF和Joint Calling参数 | GVCF and Joint Calling Parameters =====
# @click.option('--gvcf/--no-gvcf',
#               default=True,
#               help='🧬 输出GVCF格式 | Output GVCF format (default: GVCF, use --no-gvcf for VCF)')
# @click.option('--no-joint-calling', 'joint_calling',  # 注意这里：映射到 joint_calling
#               is_flag=True,           # ✅ 改为 True
#               flag_value=False,       # ✅ 添加：当提供 --no-joint-calling 时，joint_calling=False
#               default=True,           # ✅ 添加：默认值为 True
#               help='🔗 禁用Joint Calling | Disable Joint Calling (enabled by default)')
# @click.option('--combined-output',
#               default='combined.g.vcf',
#               help='📜 Joint Calling输出文件名 | Joint Calling output filename (default: combined.g.vcf)')

# # ===== 质量控制参数 | Quality Control Parameters =====
# @click.option('--min-confidence',
#               default=30,
#               type=int,
#               help='🎯 最小置信度阈值 | Minimum confidence threshold (default: 30)')
# @click.option('--min-base-quality',
#               default=20,
#               type=int,
#               help='⭐ 最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
# @click.option('--ploidy',
#               default=2,
#               type=int,
#               help='🧬 倍性 | Ploidy (default: 2)')
# @click.option('--pcr-indel-model',
#               default='CONSERVATIVE',
#               type=str,
#               help='🔬 PCR indel模型 | PCR indel model (default: CONSERVATIVE)')

# # ===== 文件模式参数 | File Pattern Parameters =====
# @click.option('--read1-pattern',
#               default='*_1.clean.fq.gz',
#               type=str,
#               help='🔍 R1文件匹配模式 | R1 file pattern (default: *_1.clean.fq.gz)')
# @click.option('--read2-pattern',
#               default='*_2.clean.fq.gz',
#               type=str,
#               help='🔍 R2文件匹配模式 | R2 file pattern (default: *_2.clean.fq.gz)')
# def parabricks(input_dir, output_dir, reference, threads, parabricks_path, tmp_dir,
#                gvcf, joint_calling, combined_output,
#                min_confidence, min_base_quality, ploidy, pcr_indel_model,
#                read1_pattern, read2_pattern):
#     """
#     🚀 Parabricks WGS批处理分析工具 | Parabricks WGS Batch Analysis Tool
    
#     基于NVIDIA Parabricks的GPU加速全基因组测序(WGS)批处理分析工具，
#     提供从FASTQ到VCF/GVCF的完整分析流程，支持Joint Calling联合变异检测。
    
#     ✨ 功能特点 | Features:
    
#     \b
#     🚀 GPU加速性能:
#        • GPU加速的高性能计算
#        • 比传统CPU工具快10-50倍
#        • 支持多GPU并行处理
#        • 自动负载均衡优化
    
#     \b
#     🧬 GVCF和Joint Calling支持:
#        • GVCF格式输出（推荐用于群体分析）
#        • Joint Calling联合变异检测
#        • 提高低频变异检测准确性
#        • 支持大规模队列分析
    
#     \b
#     📦 完整分析流程:
#        • FASTQ质量评估和预处理
#        • GPU加速序列比对(BWA-MEM)
#        • 读段排序和去重处理
#        • 碱基质量分数重校准
#        • GPU加速变异检测
#        • GVCF/VCF文件生成
#        • Joint Calling联合检测
    
#     🔄 分析流程 | Analysis Pipeline:
    
#     \b
#     标准流程（VCF输出）:
#     1. 📁 FASTQ文件配对和验证
#     2. 🚀 GPU加速序列比对
#     3. 📊 读段排序和去重
#     4. 🎯 碱基质量重校准
#     5. 🧬 变异检测生成VCF
#     6. 📄 质量过滤和统计
    
#     \b
#     GVCF流程（推荐用于群体分析）:
#     1. 📁 FASTQ文件配对和验证
#     2. 🚀 GPU加速序列比对
#     3. 📊 读段排序和去重
#     4. 🎯 碱基质量重校准
#     5. 🧬 变异检测生成GVCF
#     6. 🔗 Joint Calling合并GVCF（可选）
#     7. 📄 质量过滤和统计
    
#     🎯 应用场景 | Use Cases:
    
#     \b
#     • 🧬 大规模群体基因组学研究
#     • 🏥 临床基因组学诊断流程
#     • 🔬 基因组关联研究(GWAS)
#     • 🌱 农业基因组育种项目
#     • 🧪 药物基因组学研究
#     • 📈 进化基因组学分析
#     • 💊 精准医疗数据处理
#     • 👥 队列测序项目
    
#     💡 使用示例 | Usage Examples:
    
#     \b
#     # 🎯 基本WGS批处理分析（默认GVCF输出）
#     biopytools parabricks -i /path/to/clean/data \\
#         -o /path/to/output -r /path/to/reference.fa
    
#     \b
#     # 📄 输出VCF格式（传统模式）
#     biopytools parabricks -i ./clean_data -o ./results \\
#         -r ./reference.fa --no-gvcf
    
#     \b
#     # 🔗 启用Joint Calling（多样本联合检测）
#     biopytools parabricks -i ./cohort_fastq -o ./cohort_results \\
#         -r ./genome.fa --gvcf --joint-calling
    
#     \b
#     # 📜 自定义Joint Calling输出文件名
#     biopytools parabricks -i ./samples -o ./output \\
#         -r ./ref.fa --gvcf --joint-calling \\
#         --combined-output "my_cohort.g.vcf"
    
#     \b
#     # ⚙️ 高性能配置（多线程+自定义程序路径）
#     biopytools parabricks -i ./clean_data -o ./parabricks_results \\
#         -r ./reference.fa -t 128 --parabricks-path /custom/path
    
#     \b
#     # 🎛️ 自定义质量控制参数
#     biopytools parabricks -i /data/clean -o /results \\
#         -r /genome/ref.fa --min-confidence 35 --min-base-quality 25
    
#     \b
#     # 💾 使用快速临时目录
#     biopytools parabricks -i ./fastq_files -o ./wgs_output \\
#         -r ./genome.fa --tmp-dir /nvme/tmp
    
#     \b
#     # 🏥 临床级别高精度GVCF分析
#     biopytools parabricks -i /clinical/samples -o /clinical/results \\
#         -r /clinical/ref.fa --gvcf --min-confidence 40 \\
#         --min-base-quality 25 --pcr-indel-model AGGRESSIVE
    
#     \b
#     # 👥 大型队列Joint Calling分析
#     biopytools parabricks -i /cohort/fastq -o /cohort/results \\
#         -r /ref/genome.fa -t 256 --gvcf --joint-calling \\
#         --combined-output "population_1000.g.vcf"
    
#     🧬 GVCF vs VCF 详解 | GVCF vs VCF Details:
    
#     \b
#     📄 VCF（Variant Call Format）:
#        • 仅记录变异位点
#        • 文件较小，处理快速
#        • 适合单样本分析
#        • 不利于后续样本添加
#        • 低频变异可能丢失
    
#     \b
#     🧬 GVCF（Genomic VCF）:
#        • 记录所有位点信息（变异+参考）
#        • 文件较大，更全面
#        • 适合群体/队列分析
#        • 支持后续样本增量添加
#        • 提高低频变异检测能力
#        • 推荐用于Joint Calling
    
#     \b
#     🔗 Joint Calling优势:
#        • 提高变异检测准确性
#        • 更好地识别低频变异
#        • 统一的变异质量评分
#        • 减少批次效应
#        • 适合大规模群体研究
    
#     \b
#     使用建议:
#        • 单样本分析 → VCF (--no-gvcf)
#        • 队列分析 → GVCF + Joint Calling
#        • 计划扩展样本 → GVCF
#        • 临床诊断 → VCF (快速结果)
#        • 科研项目 → GVCF (完整信息)
    
#     📁 输入文件要求 | Input File Requirements:
    
#     \b
#     📂 目录结构要求:
#        • 输入目录包含配对的FASTQ文件
#        • 文件命名遵循标准配对模式
#        • 支持压缩格式(.fq.gz, .fastq.gz)
#        • 建议预先进行质量控制
    
#     \b
#     📝 FASTQ文件命名规范:
#        默认模式:
#        • R1文件: sample_1.clean.fq.gz
#        • R2文件: sample_2.clean.fq.gz
       
#        自定义模式示例:
#        • --read1-pattern "*_R1_*.fastq.gz"
#        • --read2-pattern "*_R2_*.fastq.gz"
    
#     \b
#     🧬 参考基因组要求:
#        • 标准FASTA格式
#        • 建议使用BWA索引预处理
#        • 常用版本: GRCh38, hg19
#        • 确保与研究目标一致
    
#     📊 输出文件结构 | Output File Structure:
    
#     \b
#     标准VCF输出（--no-gvcf）:
#        📂 output_dir/
#        ├── bam/
#        │   ├── sample1.sorted.bam
#        │   └── sample1.sorted.bam.bai
#        ├── vcf/
#        │   ├── sample1.vcf.gz
#        │   └── sample1.vcf.gz.tbi
#        └── tmp/ (临时文件)
    
#     \b
#     GVCF输出（默认）:
#        📂 output_dir/
#        ├── bam/
#        │   ├── sample1.sorted.bam
#        │   └── sample1.sorted.bam.bai
#        ├── vcf/
#        │   ├── sample1.g.vcf.gz       # 个体GVCF
#        │   ├── sample1.g.vcf.gz.tbi
#        │   ├── sample2.g.vcf.gz
#        │   └── sample2.g.vcf.gz.tbi
#        └── tmp/
    
#     \b
#     GVCF + Joint Calling输出:
#        📂 output_dir/
#        ├── bam/ (同上)
#        ├── vcf/
#        │   ├── sample1.g.vcf.gz       # 个体GVCF
#        │   ├── sample2.g.vcf.gz
#        │   ├── ...
#        │   ├── combined.g.vcf.gz      # 合并的GVCF
#        │   └── combined.g.vcf.gz.tbi
#        └── tmp/
    
#     \b
#     📈 统计报告:
#        • parabricks_summary.txt: 批处理总结
#        • processing_statistics.txt: 处理统计
#        • parabricks_analysis.log: 详细日志
    
#     ⚡ 系统和硬件要求 | System Requirements:
    
#     \b
#     🎮 GPU硬件:
#        • NVIDIA GPU with Compute Capability 6.0+
#        • 推荐: Tesla V100, A100, RTX 3090/4090
#        • 显存: 最少8GB，推荐16GB+
#        • 多GPU可显著提升性能
    
#     \b
#     💻 系统配置:
#        • RAM: 至少64GB，大规模分析建议128GB+
#        • 存储: NVMe SSD推荐
#        • CPU: 多核支持并行IO
#        • 空间: 输入数据的3-5倍
    
#     \b
#     🔧 软件依赖:
#        • NVIDIA GPU Driver (最新版本)
#        • CUDA Toolkit (11.0+)
#        • NVIDIA Parabricks (3.0+)
#        • 容器环境(Singularity/Docker)
    
#     ⚙️ 参数详解 | Parameter Details:
    
#     \b
#     🧬 GVCF参数:
#        --gvcf (默认):
#        ✅ 输出GVCF格式
#        ✅ 适合群体分析
#        ✅ 支持Joint Calling
#        ✅ 保留完整基因组信息
       
#        --no-gvcf:
#        📄 输出标准VCF格式
#        ⚡ 文件更小，处理更快
#        🎯 适合单样本分析
    
#     \b
#     🔗 Joint Calling参数:
#        --joint-calling:
#        • 启用多样本联合检测
#        • 需要配合--gvcf使用
#        • 自动合并所有成功样本的GVCF
#        • 提高变异检测准确性
       
#        --combined-output:
#        • 指定合并GVCF的文件名
#        • 默认: combined.g.vcf
#        • 建议使用描述性命名
#        • 示例: cohort_100samples.g.vcf
    
#     \b
#     🎯 质量控制参数:
#        --min-confidence (默认30):
#        • 变异检测最小置信度
#        • 范围: 10-50
#        • 临床级别建议: 35-40
       
#        --min-base-quality (默认20):
#        • 碱基调用最小质量
#        • 范围: 10-30
#        • 高精度分析建议: 25+
       
#        --ploidy (默认2):
#        • 基因组倍性
#        • 人类通常为2(二倍体)
#        • 特殊研究可调整
       
#        --pcr-indel-model:
#        • CONSERVATIVE: 保守模式(默认)
#        • AGGRESSIVE: 激进模式
#        • NONE: 不校正
    
#     🚀 性能优化建议 | Performance Tips:
    
#     \b
#     🎮 GPU优化:
#        • 确保GPU内存充足
#        • 多GPU会自动负载均衡
#        • 监控GPU利用率
#        • 避免CPU-GPU瓶颈
    
#     \b
#     🗄️ 存储优化:
#        • 使用NVMe SSD作为临时目录
#        • 确保输入输出在高速存储
#        • 预估空间: 输入数据×3-5
#        • GVCF比VCF占用更多空间
    
#     \b
#     💾 内存优化:
#        • 线程数不要超过内存限制
#        • 大基因组分析建议64GB+
#        • Joint Calling需要更多内存
#        • 监控避免交换
    
#     🛠️ 故障排除 | Troubleshooting:
    
#     \b
#     常见问题:
#        1️⃣ GPU内存不足
#           → 减少并行样本数
#           → 使用更大显存GPU
       
#        2️⃣ Joint Calling失败
#           → 检查所有GVCF文件完整
#           → 确保有足够内存
#           → 验证参考基因组一致
       
#        3️⃣ GVCF文件过大
#           → 这是正常的，GVCF包含更多信息
#           → 考虑压缩存储
#           → 使用--no-gvcf输出VCF
       
#        4️⃣ 配对文件缺失
#           → 检查R1/R2文件配对
#           → 验证文件命名模式
    
#     🏆 最佳实践 | Best Practices:
    
#     \b
#     1️⃣ 数据准备:
#        • 使用高质量clean FASTQ
#        • 预先去除adapters
#        • 验证文件完整性
#        • 确保配对一致
    
#     \b
#     2️⃣ 分析策略选择:
#        单样本诊断:
#        → --no-gvcf (快速VCF)
       
#        队列研究:
#        → --gvcf --joint-calling
       
#        可能扩展项目:
#        → --gvcf (保留扩展性)
    
#     \b
#     3️⃣ Joint Calling建议:
#        • 样本数>10时使用
#        • 确保样本质量一致
#        • 使用相同参考基因组
#        • 预留足够内存和存储
    
#     \b
#     4️⃣ 结果验证:
#        • 检查BAM比对统计
#        • 验证VCF/GVCF质量分布
#        • 比较样本间统计指标
#        • 使用已知位点验证
    
#     📚 相关资源 | Resources:
    
#     \b
#     • NVIDIA Parabricks官网
#     • GATK Best Practices
#     • GVCF格式规范
#     • Joint Calling教程
#     • GPU加速基因组学
    
#     ⚠️ 注意事项 | Important Notes:
    
#     \b
#     • GVCF文件比VCF大5-10倍，预留足够空间
#     • Joint Calling需要所有样本使用相同参考基因组
#     • 临床应用建议使用VCF获得快速结果
#     • 科研项目推荐GVCF保留完整信息
#     • 大型队列(>100样本)Joint Calling需要大内存
#     • 定期更新Parabricks获得最新优化
#     """
    
#     # 🚀 懒加载
#     parabricks_main = _lazy_import_parabricks_main()
    
#     # 构建参数列表
#     args = ['parabricks.py']
    
#     # 必需参数
#     args.extend(['-i', input_dir])
#     args.extend(['-o', output_dir])
#     args.extend(['-r', reference])
    
#     # 性能参数
#     if threads != 88:
#         args.extend(['-t', str(threads)])
#     if parabricks_path != '/share/apps/containers/parabricks.sif':
#         args.extend(['--parabricks-path', parabricks_path])
#     if tmp_dir:
#         args.extend(['--tmp-dir', tmp_dir])
    
#     # GVCF和Joint Calling参数
#     if not gvcf:  # 默认是True，如果为False则添加--no-gvcf
#         args.append('--no-gvcf')
#     if not joint_calling:  # ✅ 修改：如果为False才添加--no-joint-calling
#         args.append('--no-joint-calling')
#     if combined_output != 'combined.g.vcf':
#         args.extend(['--combined-output', combined_output])
        
#     # 质量控制参数
#     if min_confidence != 30:
#         args.extend(['--min-confidence', str(min_confidence)])
#     if min_base_quality != 20:
#         args.extend(['--min-base-quality', str(min_base_quality)])
#     if ploidy != 2:
#         args.extend(['--ploidy', str(ploidy)])
#     if pcr_indel_model != 'CONSERVATIVE':
#         args.extend(['--pcr-indel-model', pcr_indel_model])
    
#     # 文件模式参数
#     if read1_pattern != '*_1.clean.fq.gz':
#         args.extend(['--read1-pattern', read1_pattern])
#     if read2_pattern != '*_2.clean.fq.gz':
#         args.extend(['--read2-pattern', read2_pattern])
    
#     # 保存并恢复sys.argv
#     original_argv = sys.argv
#     sys.argv = args
    
#     try:
#         parabricks_main()
#     except SystemExit as e:
#         if e.code != 0:
#             sys.exit(e.code)
#     except KeyboardInterrupt:
#         click.echo("\n🛑 Parabricks WGS分析被用户中断 | Analysis interrupted by user", err=True)
#         sys.exit(1)
#     except Exception as e:
#         click.echo(f"💥 Parabricks WGS分析失败 | Analysis failed: {e}", err=True)
#         sys.exit(1)
#     finally:
#         sys.argv = original_argv

"""
🚀 Parabricks WGS分析命令 | Parabricks WGS Analysis Command
高级优化版本：支持灵活的工作流程控制和分步骤执行
"""

import click
import sys
import os


def _lazy_import_parabricks_main():
    """懒加载Parabricks main函数 | Lazy load Parabricks main function"""
    try:
        from ...parabricks.main import main as parabricks_main
        return parabricks_main
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


@click.command(
    short_help='🚀 GPU加速WGS分析工具：支持灵活的工作流程控制',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='📁 输入目录路径 | Input directory path (containing FASTQ files)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='📂 输出目录路径 | Output directory path')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='🧬 参考基因组文件路径 | Reference genome file path')

# ===== 工作流程参数 | Workflow Parameters =====
@click.option('--workflow', '-w',
              type=click.Choice(['fq2bam', 'haplotypecaller', 'genotypegvcf', 'all']),
              default='all',
              help='🔀 工作流程 | Workflow (default: all)')

# ===== 性能参数 | Performance Parameters =====
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--parabricks-path',
              default='/share/apps/containers/parabricks.sif',
              type=str,
              help='🔧 Parabricks程序路径 | Parabricks program path')
@click.option('--tmp-dir',
              type=click.Path(),
              help='💾 临时目录路径 | Temporary directory path')

# ===== GVCF和Joint Calling参数 | GVCF and Joint Calling Parameters =====
@click.option('--gvcf/--no-gvcf',
              default=True,
              help='🧬 输出GVCF格式 | Output GVCF format (default: GVCF)')
@click.option('--joint-calling/--no-joint-calling',
              default=True,
              help='🔗 启用Joint Calling | Enable Joint Calling (default: enabled)')
@click.option('--combined-output',
              default='combined.g.vcf',
              help='📜 Joint Calling输出文件名 | Joint Calling output filename (default: combined.g.vcf)')

# ===== 质量控制参数 | Quality Control Parameters =====
@click.option('--min-confidence',
              default=30,
              type=int,
              help='🎯 最小置信度阈值 | Minimum confidence threshold (default: 30)')
@click.option('--min-base-quality',
              default=20,
              type=int,
              help='⭐ 最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
@click.option('--ploidy',
              default=2,
              type=int,
              help='🧬 倍性 | Ploidy (default: 2)')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              type=str,
              help='🔬 PCR indel模型 | PCR indel model (default: CONSERVATIVE)')

# ===== 文件模式参数 | File Pattern Parameters =====
@click.option('--read1-pattern',
              default='*_1.clean.fq.gz',
              type=str,
              help='🔍 R1文件匹配模式 | R1 file pattern (default: *_1.clean.fq.gz)')
@click.option('--read2-pattern',
              default='*_2.clean.fq.gz',
              type=str,
              help='🔍 R2文件匹配模式 | R2 file pattern (default: *_2.clean.fq.gz)')
def parabricks(input_dir, output_dir, reference, workflow, threads, parabricks_path, tmp_dir,
               gvcf, joint_calling, combined_output,
               min_confidence, min_base_quality, ploidy, pcr_indel_model,
               read1_pattern, read2_pattern):
    """
    🚀 Parabricks WGS批处理分析工具 | Parabricks WGS Batch Analysis Tool
    
    基于NVIDIA Parabricks的GPU加速全基因组测序分析工具，支持灵活的
    工作流程控制，可以分步骤执行或运行完整流程。
    
    ✨ 核心功能 | Core Features:
    
    \b
    🔀 灵活的工作流程控制:
       • fq2bam: FASTQ到BAM的比对流程
       • haplotypecaller: BAM到VCF/GVCF的变异检测
       • genotypegvcf: GVCF合并和Joint Calling
       • all: 完整的端到端分析流程
    
    \b
    🚀 GPU加速性能:
       • 比传统CPU工具快10-50倍
       • 支持多GPU并行处理
       • 自动负载均衡优化
       • 高效的内存管理
    
    \b
    🧬 完整的变异检测:
       • GVCF格式输出
       • Joint Calling联合检测
       • 提高低频变异准确性
       • 支持大规模队列分析
    
    🔀 工作流程详解 | Workflow Details:
    
    \b
    📌 fq2bam - 比对流程:
       输入: FASTQ文件对
       输出: BAM文件
       步骤:
       1. GPU加速序列比对(BWA-MEM)
       2. 读段排序(Sort)
       3. 标记重复(MarkDuplicates)
       4. 碱基质量重校准(BQSR)
       
       适用场景:
       • 仅需要比对结果
       • 分步骤处理大批量数据
       • 比对后需要自定义处理
    
    \b
    📌 haplotypecaller - 变异检测:
       输入: BAM文件
       输出: VCF或GVCF文件
       步骤:
       1. 读取已有的BAM文件
       2. GPU加速变异检测
       3. 生成VCF/GVCF文件
       
       适用场景:
       • 已有BAM文件需要重新call variants
       • 使用不同参数重新检测
       • 分批次处理变异检测
    
    \b
    📌 genotypegvcf - Joint Calling:
       输入: 多个GVCF文件
       输出: 合并的GVCF文件
       步骤:
       1. 自动检测所有GVCF文件
       2. 合并多样本GVCF
       3. 生成联合检测结果
       
       适用场景:
       • 队列分析的最后合并步骤
       • 增量添加新样本后重新合并
       • 仅需要Joint Calling结果
    
    \b
    📌 all - 完整流程(默认):
       输入: FASTQ文件对
       输出: BAM + GVCF + Joint Calling
       步骤:
       1. fq2bam全部步骤
       2. haplotypecaller全部步骤
       3. genotypegvcf合并(可选)
       
       适用场景:
       • 标准的端到端分析
       • 一次性完成所有步骤
       • 推荐用于常规分析
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 完整流程分析（最常用）
    biopytools parabricks -i /data/fastq -o /results -r /ref/genome.fa
    
    \b
    # 🔀 分步骤执行策略
    
    # 步骤1: 仅执行比对
    biopytools parabricks -i /data/fastq -o /results -r /ref/genome.fa \\
        --workflow fq2bam
    
    # 步骤2: 基于BAM执行变异检测
    biopytools parabricks -i /data/fastq -o /results -r /ref/genome.fa \\
        --workflow haplotypecaller
    
    # 步骤3: 合并所有GVCF
    biopytools parabricks -i /data/fastq -o /results -r /ref/genome.fa \\
        --workflow genotypegvcf
    
    \b
    # 📊 不同应用场景
    
    # 场景1: 大批量数据分批处理
    # 先完成所有样本的比对
    biopytools parabricks -i /batch1 -o /results -r /ref.fa -w fq2bam
    biopytools parabricks -i /batch2 -o /results -r /ref.fa -w fq2bam
    # 再统一进行变异检测
    biopytools parabricks -i /batch1 -o /results -r /ref.fa -w haplotypecaller
    biopytools parabricks -i /batch2 -o /results -r /ref.fa -w haplotypecaller
    # 最后Joint Calling
    biopytools parabricks -i /batch1 -o /results -r /ref.fa -w genotypegvcf
    
    \b
    # 场景2: 重新检测变异（已有BAM）
    biopytools parabricks -i /data -o /results -r /ref.fa \\
        -w haplotypecaller --min-confidence 40
    
    \b
    # 场景3: 仅执行Joint Calling
    biopytools parabricks -i /data -o /results -r /ref.fa \\
        -w genotypegvcf --combined-output "cohort_500.g.vcf"
    
    \b
    # 场景4: 自定义参数的完整流程
    biopytools parabricks -i /data -o /results -r /ref.fa \\
        -w all -t 128 --gvcf --joint-calling \\
        --min-confidence 35 --ploidy 2
    
    \b
    # 场景5: 禁用Joint Calling的GVCF输出
    biopytools parabricks -i /data -o /results -r /ref.fa \\
        --gvcf --no-joint-calling
    
    \b
    # 场景6: 输出VCF而不是GVCF
    biopytools parabricks -i /data -o /results -r /ref.fa \\
        --no-gvcf -w haplotypecaller
    
    🎯 工作流程选择指南 | Workflow Selection Guide:
    
    \b
    选择fq2bam的情况:
    ✅ 需要在变异检测前进行额外的BAM处理
    ✅ 分批次处理大量样本以管理资源
    ✅ 需要先完成所有样本的比对
    ✅ 比对和变异检测分开进行质控
    
    \b
    选择haplotypecaller的情况:
    ✅ 已有BAM文件需要重新call variants
    ✅ 使用不同参数重新检测变异
    ✅ 比对已在其他流程中完成
    ✅ 需要测试不同的过滤参数
    
    \b
    选择genotypegvcf的情况:
    ✅ 已有多个GVCF文件需要合并
    ✅ 队列分析的最后合并步骤
    ✅ 增量添加新样本后重新Joint Calling
    ✅ 仅需要更新Joint Calling结果
    
    \b
    选择all的情况:
    ✅ 标准的端到端分析
    ✅ 一次性完成所有步骤
    ✅ 小到中等规模的项目
    ✅ 推荐用于大多数常规分析
    
    📁 输入文件要求 | Input Requirements:
    
    \b
    fq2bam workflow需要:
       • 配对的FASTQ文件
       • 文件命名符合模式
       • 支持压缩格式(.fq.gz)
    
    \b
    haplotypecaller workflow需要:
       • 已存在的BAM文件(在output_dir/bam/)
       • BAM文件已完成排序和索引
       • BAM文件命名: {sample}.sorted.bam
    
    \b
    genotypegvcf workflow需要:
       • 已存在的GVCF文件(在output_dir/vcf/)
       • GVCF文件命名: {sample}.g.vcf.gz
       • 自动检测所有符合条件的GVCF
    
    \b
    all workflow需要:
       • 配对的FASTQ文件
       • 其他步骤自动生成
    
    📊 输出文件结构 | Output Structure:
    
    \b
    fq2bam workflow输出:
       📂 output_dir/
       ├── bam/
       │   ├── sample1.sorted.bam
       │   ├── sample1.sorted.bam.bai
       │   ├── sample2.sorted.bam
       │   └── sample2.sorted.bam.bai
       └── tmp/
    
    \b
    haplotypecaller workflow输出:
       📂 output_dir/
       ├── bam/ (已存在)
       ├── vcf/
       │   ├── sample1.g.vcf.gz (如果--gvcf)
       │   ├── sample1.g.vcf.gz.tbi
       │   ├── sample2.g.vcf.gz
       │   └── sample2.g.vcf.gz.tbi
       └── tmp/
    
    \b
    genotypegvcf workflow输出:
       📂 output_dir/
       ├── vcf/
       │   ├── sample1.g.vcf.gz (已存在)
       │   ├── sample2.g.vcf.gz (已存在)
       │   ├── combined.g.vcf.gz (新增)
       │   └── combined.g.vcf.gz.tbi
       └── tmp/
    
    \b
    all workflow输出:
       📂 output_dir/
       ├── bam/
       │   └── *.sorted.bam
       ├── vcf/
       │   ├── *.g.vcf.gz (个体GVCF)
       │   └── combined.g.vcf.gz (如果--joint-calling)
       └── tmp/
    
    ⚡ 性能优化建议 | Performance Tips:
    
    \b
    🎮 GPU资源管理:
       • fq2bam最消耗GPU资源
       • haplotypecaller次之
       • genotypegvcf主要消耗CPU和内存
       • 合理规划GPU使用时间
    
    \b
    💾 存储空间规划:
       • fq2bam需要: FASTQ大小×2-3倍
       • haplotypecaller需要: BAM大小×0.5-1倍
       • genotypegvcf需要: 所有GVCF大小×1.5倍
       • 总空间: 输入数据×3-5倍
    
    \b
    🔄 批处理策略:
       方案1: 串行完整流程
       → 适合小规模项目(<50样本)
       → 使用--workflow all
       
       方案2: 分步并行处理
       → 适合大规模项目(>100样本)
       → 先完成所有fq2bam
       → 再批量haplotypecaller
       → 最后统一genotypegvcf
       
       方案3: 混合策略
       → fq2bam和haplotypecaller分批
       → 定期执行genotypegvcf更新
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1️⃣ haplotypecaller找不到BAM文件
          → 确认已运行fq2bam workflow
          → 检查bam目录下是否有*.sorted.bam
          → 验证文件命名正确
       
       2️⃣ genotypegvcf找不到GVCF文件
          → 确认已运行haplotypecaller workflow
          → 检查vcf目录下是否有*.g.vcf.gz
          → 确认使用了--gvcf参数
       
       3️⃣ Joint Calling失败
          → 检查所有GVCF文件完整性
          → 确保有足够内存
          → 验证参考基因组一致
       
       4️⃣ 工作流程选择错误
          → fq2bam需要FASTQ输入
          → haplotypecaller需要已有BAM
          → genotypegvcf需要已有GVCF
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 标准项目流程:
       # 推荐使用all workflow
       biopytools parabricks -i /data -o /results -r /ref.fa -w all
    
    \b
    2️⃣ 大规模项目流程:
       # 第一阶段: 比对所有样本
       for batch in batch*; do
           biopytools parabricks -i $batch -o /results -r /ref.fa -w fq2bam
       done
       
       # 第二阶段: 变异检测
       biopytools parabricks -i /all -o /results -r /ref.fa -w haplotypecaller
       
       # 第三阶段: Joint Calling
       biopytools parabricks -i /all -o /results -r /ref.fa -w genotypegvcf
    
    \b
    3️⃣ 增量分析流程:
       # 初始队列
       biopytools parabricks -i /batch1 -o /results -r /ref.fa -w all
       
       # 新增样本
       biopytools parabricks -i /new_samples -o /results -r /ref.fa -w fq2bam
       biopytools parabricks -i /new_samples -o /results -r /ref.fa -w haplotypecaller
       
       # 重新Joint Calling（包含所有样本）
       biopytools parabricks -i /all -o /results -r /ref.fa -w genotypegvcf
    
    \b
    4️⃣ 参数优化流程:
       # 测试不同参数
       biopytools parabricks -i /test -o /test1 -r /ref.fa \\
           -w haplotypecaller --min-confidence 30
       biopytools parabricks -i /test -o /test2 -r /ref.fa \\
           -w haplotypecaller --min-confidence 40
       # 比较结果后选择最佳参数
    
    📚 相关资源 | Resources:
    
    \b
    • NVIDIA Parabricks官网
    • GATK Best Practices
    • GPU加速基因组学
    • WGS分析流程指南
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 确保每个workflow的输入文件已准备好
    • genotypegvcf会自动检测所有GVCF文件
    • 分步执行时注意文件路径一致性
    • Joint Calling默认启用，可用--no-joint-calling禁用
    • 不同workflow可能需要不同的系统资源
    • 建议先用小数据集测试工作流程
    • 大规模分析时合理规划GPU使用时间
    """
    
    # 🚀 懒加载
    parabricks_main = _lazy_import_parabricks_main()
    
    # 构建参数列表
    args = ['parabricks.py']
    
    # 必需参数
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])
    
    # 工作流程参数
    if workflow != 'all':
        args.extend(['-w', workflow])
    
    # 性能参数
    if threads != 88:
        args.extend(['-t', str(threads)])
    if parabricks_path != '/share/apps/containers/parabricks.sif':
        args.extend(['--parabricks-path', parabricks_path])
    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])
    
    # GVCF和Joint Calling参数
    if not gvcf:
        args.append('--no-gvcf')
    if not joint_calling:
        args.append('--no-joint-calling')
    if combined_output != 'combined.g.vcf':
        args.extend(['--combined-output', combined_output])
    
    # 质量控制参数
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])
    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])
    
    # 文件模式参数
    if read1_pattern != '*_1.clean.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
    if read2_pattern != '*_2.clean.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        parabricks_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 Parabricks分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 Parabricks分析失败 | Analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv