# """
# GTX WGS分析命令 | GTX WGS Analysis Command
# 优化版本：使用懒加载解决响应速度问题
# """

# import click
# import sys
# import os


# def _lazy_import_gtx_main():
#     """懒加载gtx main函数 | Lazy load gtx main function"""
#     try:
#         from ...gtx.main import main as gtx_main
#         return gtx_main
#     except ImportError as e:
#         click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
#         sys.exit(1)


# def _is_help_request():
#     """检查是否是帮助请求 | Check if this is a help request"""
#     help_flags = {'-h', '--help'}
#     return any(arg in help_flags for arg in sys.argv)


# def _validate_input_dir(dir_path):
#     """验证输入目录是否存在（仅在非帮助模式下）| Validate input directory existence (only in non-help mode)"""
#     if not _is_help_request():
#         if not os.path.exists(dir_path):
#             raise click.BadParameter(f"输入目录不存在 | Input directory does not exist: {dir_path}")
#         if not os.path.isdir(dir_path):
#             raise click.BadParameter(f"输入路径不是目录 | Input path is not a directory: {dir_path}")
#     return dir_path


# def _validate_reference_file(file_path):
#     """验证参考基因组文件是否存在（仅在非帮助模式下）| Validate reference file existence (only in non-help mode)"""
#     if not _is_help_request() and not os.path.exists(file_path):
#         raise click.BadParameter(f"参考基因组文件不存在 | Reference genome file does not exist: {file_path}")
#     return file_path


# @click.command(short_help='运行GTX WGS流程',
#                context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# @click.option('--input-dir', '-i',
#               required=True,
#               callback=lambda ctx, param, value: _validate_input_dir(value) if value else None,
#               help='📂 输入目录路径 (包含clean FASTQ文件) | Input directory path (containing clean FASTQ files)')
# @click.option('--output-dir', '-o',
#               required=True,
#               type=click.Path(),
#               help='📤 输出目录路径 | Output directory path')
# @click.option('--reference', '-r',
#               required=True,
#               callback=lambda ctx, param, value: _validate_reference_file(value) if value else None,
#               help='🧬 参考基因组文件路径 | Reference genome file path')
# @click.option('--threads', '-t',
#               default=88,
#               type=int,
#               help='🧵 线程数 | Number of threads (default: 88)')
# @click.option('--gtx-path',
#               default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx',
#               type=click.Path(),
#               help='💻 GTX程序路径 | GTX program path (default: /share/apps/gtx/GTX.CAT_2.2.1/bin/gtx)')
# @click.option('--tmp-dir',
#               type=click.Path(),
#               help='🗑️ 临时目录路径 (默认使用输出目录下的tmp) | Temporary directory path (default: tmp under output directory)')
# @click.option('--min-confidence',
#               default=30,
#               type=int,
#               help='🎯 最小置信度阈值 | Minimum confidence threshold (default: 30)')
# @click.option('--min-base-quality',
#               default=20,
#               type=int,
#               help='✨ 最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
# @click.option('--ploidy',
#               default=2,
#               type=int,
#               help='🧬 倍性 | Ploidy (default: 2)')
# @click.option('--pcr-indel-model',
#               default='CONSERVATIVE',
#               type=str,
#               help='🔬 PCR indel模型 | PCR indel model (default: CONSERVATIVE)')
# @click.option('--read1-pattern',
#               default='*_1.fq.gz',
#               type=str,
#               help='📄 R1文件匹配模式 | R1 file pattern (default: *_1.fq.gz)')
# @click.option('--read2-pattern',
#               default='*_2.fq.gz',
#               type=str,
#               help='📄 R2文件匹配模式 | R2 file pattern (default: *_2.fq.gz)')
# def gtx(input_dir, output_dir, reference, threads, gtx_path, tmp_dir,
#         min_confidence, min_base_quality, ploidy, pcr_indel_model,
#         read1_pattern, read2_pattern):
#     """
#     🧬 GTX WGS批处理分析工具
    
#     使用GTX软件对FASTQ文件进行全基因组测序分析，
#     包括比对、变异检测和质量控制等步骤。
    
#     示例 | Examples:
    
#     \b
#     # 🚀 基本分析
#     biopytools gtx -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
    
#     \b
#     # ⚡ 高线程数分析
#     biopytools gtx -i ./clean_data -o ./gtx_results -r ./reference.fa -t 64
    
#     \b
#     # 💻 自定义GTX路径
#     biopytools gtx -i /data/clean -o /results -r /genome/ref.fa --gtx-path /custom/gtx/path
    
#     \b
#     # 🔬 自定义质量控制参数
#     biopytools gtx -i ./data -o ./results -r ./ref.fa \\
#         --min-confidence 25 --min-base-quality 15 --ploidy 2
    
#     \b
#     # 📁 自定义文件匹配模式
#     biopytools gtx -i ./data -o ./results -r ./ref.fa \\
#         --read1-pattern "*_R1.fq.gz" --read2-pattern "*_R2.fq.gz"
    
#     \b
#     # 🗑️ 指定临时目录
#     biopytools gtx -i ./data -o ./results -r ./ref.fa \\
#         --tmp-dir /tmp/gtx_analysis --threads 32
    
#     \b
#     # 🔧 完整参数示例
#     biopytools gtx -i /data/fastq -o /results/gtx \\
#         -r /genome/reference.fa -t 128 \\
#         --gtx-path /opt/gtx/bin/gtx \\
#         --min-confidence 40 --min-base-quality 25 \\
#         --pcr-indel-model AGGRESSIVE \\
#         --read1-pattern "*_1.clean.fq.gz" \\
#         --read2-pattern "*_2.clean.fq.gz"
#     """
    
#     # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
#     gtx_main = _lazy_import_gtx_main()
    
#     # 构建参数列表传递给原始main函数
#     args = ['gtx.py']
    
#     # 必需参数
#     args.extend(['-i', input_dir])
#     args.extend(['-o', output_dir])
#     args.extend(['-r', reference])
    
#     # 可选参数（只在非默认值时添加）
#     if threads != 88:
#         args.extend(['-t', str(threads)])
    
#     if gtx_path != '/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx':
#         args.extend(['--gtx-path', gtx_path])
    
#     if tmp_dir:
#         args.extend(['--tmp-dir', tmp_dir])
    
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
#     if read1_pattern != '*_1.fq.gz':
#         args.extend(['--read1-pattern', read1_pattern])
    
#     if read2_pattern != '*_2.fq.gz':
#         args.extend(['--read2-pattern', read2_pattern])
    
#     # 保存并恢复sys.argv
#     original_argv = sys.argv
#     sys.argv = args
    
#     try:
#         # 调用原始的main函数
#         gtx_main()
#     except SystemExit as e:
#         # 处理程序正常退出
#         if e.code != 0:
#             sys.exit(e.code)
#     except KeyboardInterrupt:
#         click.echo("\n⚠️ 用户中断操作 | User interrupted", err=True)
#         sys.exit(1)
#     except Exception as e:
#         click.echo(f"💥 运行错误 | Runtime error: {e}", err=True)
#         sys.exit(1)
#     finally:
#         sys.argv = original_argva

"""
GTX WGS分析命令 | GTX WGS Analysis Command
增强版本：支持Joint Calling功能，使用懒加载优化响应速度
"""

import click
import sys
import os


def _lazy_import_gtx_main():
    """懒加载gtx main函数 | Lazy load gtx main function"""
    try:
        from ...gtx.main import main as gtx_main
        return gtx_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_dir(dir_path):
    """验证输入目录是否存在（仅在非帮助模式下）| Validate input directory existence (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(dir_path):
            raise click.BadParameter(f"输入目录不存在 | Input directory does not exist: {dir_path}")
        if not os.path.isdir(dir_path):
            raise click.BadParameter(f"输入路径不是目录 | Input path is not a directory: {dir_path}")
    return dir_path


def _validate_reference_file(file_path):
    """验证参考基因组文件是否存在（仅在非帮助模式下）| Validate reference file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"参考基因组文件不存在 | Reference genome file does not exist: {file_path}")
    return file_path


def _validate_joint_output(joint_output):
    """验证Joint输出文件名格式 | Validate joint output filename format"""
    if joint_output and not joint_output.endswith('.vcf.gz'):
        raise click.BadParameter(f"Joint输出文件必须以.vcf.gz结尾 | Joint output file must end with .vcf.gz: {joint_output}")
    return joint_output


@click.command(short_help='GTX全基因组测序批处理分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_dir(value) if value else None,
              help='📂 输入目录路径 (包含clean FASTQ文件) | Input directory path (containing clean FASTQ files)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='📤 输出目录路径 | Output directory path')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_reference_file(value) if value else None,
              help='🧬 参考基因组文件路径 | Reference genome file path')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--gtx-path',
              default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx',
              type=click.Path(),
              help='💻 GTX程序路径 | GTX program path (default: /share/apps/gtx/GTX.CAT_2.2.1/bin/gtx)')
@click.option('--tmp-dir',
              type=click.Path(),
              help='🗑️  临时目录路径 (默认使用输出目录下的tmp) | Temporary directory path (default: tmp under output directory)')
@click.option('--enable-joint',
              is_flag=True,
              help='🤝 启用Joint Calling多样品联合变异检测 | Enable Joint Calling for multi-sample variant detection')
@click.option('--joint-output',
              default='merged_gtx.vcf.gz',
              callback=lambda ctx, param, value: _validate_joint_output(value) if value else None,
              help='📜 Joint calling输出VCF文件名 | Joint calling output VCF filename (default: merged_gtx.vcf.gz)')
@click.option('--joint-threads',
              default=88,
              type=int,
              help='🧵 Joint calling专用线程数 | Number of threads for joint calling (default: 88)')
@click.option('--min-confidence',
              default=30,
              type=int,
              help='🎯 最小置信度阈值 | Minimum confidence threshold (default: 30)')
@click.option('--min-base-quality',
              default=20,
              type=int,
              help='✨ 最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
@click.option('--ploidy',
              default=2,
              type=int,
              help='🧬 倍性 | Ploidy (default: 2)')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              type=click.Choice(['CONSERVATIVE', 'AGGRESSIVE']),
              help='🔬 PCR indel模型 | PCR indel model (default: CONSERVATIVE)')
@click.option('--read1-pattern',
              default='*_1.clean.fq.gz',
              type=str,
              help='📄 R1文件匹配模式 | R1 file pattern (default: *_1.clean.fq.gz)')
@click.option('--read2-pattern',
              default='*_2.clean.fq.gz',
              type=str,
              help='📄 R2文件匹配模式 | R2 file pattern (default: *_2.clean.fq.gz)')
def gtx(input_dir, output_dir, reference, threads, gtx_path, tmp_dir,
            enable_joint, joint_output, joint_threads, min_confidence, min_base_quality,
            ploidy, pcr_indel_model, read1_pattern, read2_pattern):
    """
    🧬 GTX全基因组测序(WGS)批处理分析工具 (增强版 - 支持Joint Calling)
    
    一个完整的WGS分析流程工具，基于GTX（Genome Tools eXtended）软件，支持批量样品
    处理、质量控制、变异检测和多样品联合调用(Joint Calling)功能，适用于大规模全基因
    组测序项目的数据处理和群体遗传学分析。
    
    🔄 分析流程 | Analysis Pipeline:
    
    第一阶段 - 单样品分析 | Phase 1 - Individual Sample Analysis:
    1. 📂 检测和验证输入FASTQ文件对
    2. 🧬 高效序列比对到参考基因组
    3. 📜 精确单样品变异检测
    4. 💾 生成高质量BAM和VCF文件
    5. 🔬 质量控制和统计报告
    
    第二阶段 - 联合分析 (可选) | Phase 2 - Joint Analysis (Optional):
    1. 🤝 多样品联合变异调用
    2. 📋 智能生成样品映射文件
    3. 📜 合并所有样品的变异信息
    4. ✨ 提高变异检测准确性和基因型质量
    5. 🧬 生成适合群体分析的联合VCF文件
    
    🎯 核心优势 | Core Advantages:
    - ⚡ 高性能并行处理：支持大规模多线程计算
    - 🤝 智能联合分析：Joint Calling提升变异检测精度
    - 🔬 严格质量控制：多层次质量过滤确保数据可靠性
    - 📊 详细统计报告：全面的分析结果和性能指标
    - 🏗️ 模块化设计：灵活配置适应不同项目需求
    - 🧹 智能资源管理：自动清理临时文件节省存储空间
    
    📱 应用场景 | Use Cases:
    - 🧬 全基因组关联分析(GWAS)预处理
    - 🧪 群体遗传学和进化研究
    - 🔬 疾病基因组学和精准医学
    - 🌱 农作物和家畜育种改良
    - 📊 系统发育和比较基因组学
    - 🧬 个人基因组学和族谱分析
    
    示例用法 | Example Usage:
    
    \b
    # 🚀 基本单样品WGS分析
    biopytools gtx -i /data/clean_fastq -o /results -r /genome/reference.fa
    
    \b
    # 🤝 启用Joint Calling的群体分析
    biopytools gtx -i /data/population_data -o /results \\
                       -r /genome/reference.fa --enable-joint \\
                       --joint-output population_variants.vcf.gz
    
    \b
    # ⚡ 高性能计算配置
    biopytools gtx -i ./clean_data -o ./gtx_results \\
                       -r ./reference.fa -t 128 \\
                       --joint-threads 64 --enable-joint
    
    \b
    # 🔬 自定义质量控制参数
    biopytools gtx -i ./data -o ./results -r ./ref.fa \\
                       --min-confidence 35 --min-base-quality 25 \\
                       --ploidy 2 --pcr-indel-model AGGRESSIVE \\
                       --enable-joint
    
    \b
    # 📁 自定义文件匹配模式
    biopytools gtx -i ./sequencing_data -o ./analysis_results \\
                       -r ./genome.fa \\
                       --read1-pattern "*_R1_*.fastq.gz" \\
                       --read2-pattern "*_R2_*.fastq.gz" \\
                       --enable-joint --joint-output cohort_merged.vcf.gz
    
    \b
    # 💻 完全自定义配置
    biopytools gtx -i /project/WGS_data \\
                       -o /project/gtx_analysis \\
                       -r /reference/GRCh38.fa \\
                       -t 88 --joint-threads 44 \\
                       --gtx-path /custom/gtx/bin/gtx \\
                       --tmp-dir /fast_storage/tmp \\
                       --enable-joint \\
                       --joint-output final_cohort.vcf.gz \\
                       --min-confidence 40 --min-base-quality 30 \\
                       --read1-pattern "*_1.clean.fq.gz" \\
                       --read2-pattern "*_2.clean.fq.gz"
    
    \b
    # 🌱 植物基因组分析 (四倍体)
    biopytools gtx -i /plant_data -o /plant_results \\
                       -r /plant_genome.fa --ploidy 4 \\
                       --enable-joint --joint-output plant_population.vcf.gz \\
                       -t 64 --min-confidence 25
    
    📂 输出目录结构 | Output Directory Structure:
    
    output_directory/
    ├── bam/              # 💾 BAM比对文件目录
    │   ├── sample1.sorted.bam
    │   ├── sample2.sorted.bam
    │   └── ...
    ├── vcf/              # 📜 单样品VCF文件目录
    │   ├── sample1.vcf.gz
    │   ├── sample2.vcf.gz
    │   └── ...
    ├── joint/            # 🤝 Joint calling结果目录 (如果启用)
    │   ├── merged_gtx.vcf.gz    # 联合VCF文件
    │   └── sample_map.txt       # 样品映射文件
    ├── tmp/              # 🗑️ 临时文件目录
    └── gtx_analysis_summary.txt # 📄 详细分析报告
    
    🤝 Joint Calling详解 | Joint Calling Details:
    
    什么是Joint Calling?
    - 📊 同时分析多个样品的变异检测方法
    - 🎯 基于群体信息优化每个样品的基因型调用
    - 🔍 显著提高低频变异的检测准确性
    - ✨ 减少单样品分析中的假阳性结果
    
    Joint Calling的优势:
    - 🎯 提高变异检测精度：群体信息辅助基因型判断
    - 📊 改善质量评分：基于多样品数据的置信度计算
    - 🔍 检测稀有变异：更好识别群体中的低频变异
    - 🧬 统一基因型格式：便于后续群体分析和GWAS
    - 📈 提升统计功效：为关联分析提供高质量数据
    
    何时使用Joint Calling?
    - ✅ 群体遗传学研究项目
    - ✅ GWAS关联分析预处理
    - ✅ 家系分析和遗传咨询
    - ✅ 育种项目的基因型分析
    - ✅ 进化和系统发育研究
    - ❌ 单一样品的临床诊断
    - ❌ 体细胞突变检测
    
    ⚙️ 性能优化建议 | Performance Optimization:
    
    硬件配置:
    - 💾 内存：建议≥32GB，大基因组项目需≥64GB
    - 💽 存储：使用SSD提升I/O性能，预留足够空间
    - ⚡ CPU：多核处理器，建议≥16核心
    - 🌐 网络：高速网络存储可提升大文件访问速度
    
    参数调优:
    - 🧵 --threads: 根据CPU核心数调整，避免过度订阅
    - 🤝 --joint-threads: 可适当减少以节省内存
    - 🗑️ --tmp-dir: 使用高速存储设备作为临时目录
    - 🎯 质控参数: 根据数据质量和分析目标调整阈值
    
    批处理策略:
    - 📦 分批处理：大项目可按染色体或样品数分批
    - ⏰ 错峰运行：避免系统负载高峰期运行
    - 📊 监控资源：实时监控内存和磁盘使用情况
    - 🔄 断点续传：利用输出文件检查功能避免重复计算
    
    🛠️ 故障排除 | Troubleshooting:
    
    常见问题及解决方案:
    
    1. "GTX program not found"
       - 🔍 检查GTX安装路径：ls -la /share/apps/gtx/
       - 🔧 使用--gtx-path指定正确路径
       - 📦 确认GTX软件完整安装
    
    2. "No FASTQ file pairs found"
       - 📂 验证输入目录路径和权限
       - 🔍 检查文件命名是否匹配模式
       - 📝 使用自定义--read1-pattern和--read2-pattern
    
    3. "Memory allocation failed"
       - 💾 增加系统可用内存或使用内存更大的节点
       - 🧵 减少线程数降低内存占用
       - 🗑️ 清理磁盘空间，确保有足够虚拟内存
    
    4. "Joint calling failed"
       - 🔢 确保至少有2个成功处理的样品
       - ✅ 检查单样品VCF文件完整性
       - 📋 验证样品映射文件格式
    
    5. "Reference genome index error"
       - 🧬 确认参考基因组文件完整性
       - 🔍 检查文件权限和路径
       - 🔄 重新构建基因组索引文件
    
    📊 质量控制参数说明 | Quality Control Parameters:
    
    --min-confidence (置信度阈值):
    - 🎯 控制变异调用的严格程度
    - 📈 较高值(35-50)：更严格，假阳性少但可能漏检
    - 📉 较低值(20-30)：更宽松，检出率高但假阳性多
    - 💡 推荐：群体分析用30-40，临床分析用40-50
    
    --min-base-quality (碱基质量阈值):
    - ✨ 过滤低质量碱基对变异调用的影响
    - 📊 Phred质量分数，20=99%准确性，30=99.9%准确性
    - 🔬 推荐值：一般数据20-25，高质量数据25-30
    
    --ploidy (倍性设置):
    - 🧬 指定目标物种的染色体倍数
    - 👥 人类和大多数动物：2 (二倍体)
    - 🌱 许多植物：4, 6, 8 (多倍体)
    - ⚠️ 错误设置会影响基因型调用准确性
    
    --pcr-indel-model (PCR模型):
    - 🔬 CONSERVATIVE: 保守模式，减少PCR相关假阳性
    - ⚡ AGGRESSIVE: 激进模式，提高indel检出率
    - 💡 推荐：WGS用CONSERVATIVE，靶向测序用AGGRESSIVE
    
    🎓 最佳实践 | Best Practices:
    
    项目规划:
    - 📋 制定详细的分析计划和质控标准
    - 🗂️ 建立规范的文件命名和目录结构
    - 📊 预估计算资源需求和时间安排
    - 💾 制定数据备份和存档策略
    
    数据准备:
    - 🔬 使用高质量的clean FASTQ数据
    - ✅ 预先验证数据完整性和格式
    - 🧬 确保参考基因组版本一致性
    - 📝 准备详细的样品信息表
    
    分析执行:
    - 🔄 先用小样本测试流程和参数
    - 📊 监控分析过程和资源使用
    - 📝 详细记录分析参数和环境
    - ✅ 定期检查中间结果质量
    
    结果验证:
    - 🔍 检查关键统计指标是否合理
    - 🧬 验证已知变异位点的检出情况
    - 📊 比较不同样品间的结果一致性
    - 📈 评估Joint calling前后的改进效果
    
    后续分析:
    - 🧬 使用联合VCF进行群体遗传学分析
    - 📊 结合表型数据进行关联分析
    - 🔍 进行变异注释和功能预测
    - 📈 可视化结果并生成报告
    
    引用信息 | Citation Information:
    如果在学术研究中使用此工具，请引用GTX软件和相关方法学文献。
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    gtx_main = _lazy_import_gtx_main()
    
    # 构建参数列表传递给原始main函数
    args = ['gtx_wgs.py']
    
    # 必需参数
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])
    
    # 可选参数（只在非默认值时添加）
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if gtx_path != '/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx':
        args.extend(['--gtx-path', gtx_path])
    
    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])
    
    # Joint calling参数
    if enable_joint:
        args.append('--enable-joint')
    
    if joint_output != 'merged_gtx.vcf.gz':
        args.extend(['--joint-output', joint_output])
    
    if joint_threads != 88:
        args.extend(['--joint-threads', str(joint_threads)])
    
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
    if read1_pattern != '*_1.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
    
    if read2_pattern != '*_2.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        gtx_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n⚠️ GTX WGS分析被用户中断 | GTX WGS analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 GTX WGS分析失败 | GTX WGS analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv