"""
🚀 SRA转FASTQ转换命令 | SRA to FASTQ Conversion Command
parallel-fastq-dump加速版：自动并行处理，速度提升3-5倍
"""

import click
import sys
import os


def _lazy_import_sra_main():
    """懒加载SRA转换main函数 | Lazy load SRA conversion main function"""
    try:
        from ...sra_converter.main import main as sra_main
        return sra_main
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
    short_help='🚀 SRA转FASTQ工具：使用parallel-fastq-dump自动并行加速',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='📂 输入SRA文件或文件夹路径 | Input SRA file or folder path')

# ===== 输出控制 | Output Control =====
@click.option('--output', '-o',
              default='./fastq_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./fastq_output)')

# ===== 性能参数 | Performance Parameters =====
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--tmpdir',
              type=click.Path(),
              help='💾 临时目录 | Temporary directory (use fast storage like SSD for acceleration)')

# ===== 转换参数 | Conversion Parameters =====
@click.option('--compress/--no-compress',
              default=True,
              help='🗜️  压缩输出 | Compress output to .gz format (default: compress)')
@click.option('--split/--no-split',
              'split_files',
              default=True,
              help='✂️  拆分双端测序文件 | Split paired-end reads (default: split)')

# ===== 过滤参数 | Filtering Parameters =====
@click.option('--min-len',
              type=int,
              default=0,
              help='📏 最小读长过滤 | Minimum read length filter (default: 0, no filter)')
@click.option('--clip',
              is_flag=True,
              help='✂️  剪切adapters | Clip adapters')
def sra2fastq(input, output, threads, tmpdir, compress, split_files, min_len, clip):
    """
    🚀 SRA转FASTQ高速转换工具 | SRA to FASTQ High-Speed Conversion Tool
    
    使用parallel-fastq-dump自动并行处理，显著提升转换速度！
    自动检测并充分利用CPU多核性能，让SRA转换更快更简单。
    
    ✨ 核心优势 | Key Advantages:
    
    \b
    ⚡ 自动并行加速:
       • parallel-fastq-dump: 自动多线程并行
       • 智能任务分配: 充分利用CPU资源
       • 零配置加速: 无需手动优化参数
       • 速度提升: 比传统fastq-dump快3-5倍
    
    \b
    🎯 简单易用:
       • 一键转换: 自动选择最优工具
       • 批量处理: 支持文件夹批量转换
       • 智能压缩: 自动使用gzip压缩
       • 实时进度: 清晰的转换状态显示
    
    \b
    📦 完整功能:
       • 双端数据自动拆分
       • 最小读长过滤
       • Adapter自动剪切
       • 详细统计报告
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    SRA文件格式:
       • 标准SRA格式 (*.sra)
       • NCBI/ENA/DDBJ数据库文件
       • 单个文件或批量文件夹
       • 命名规范: SRR/ERR/DRR前缀
    
    \b
    文件来源:
       • NCBI SRA数据库下载
       • 使用prefetch预下载
       • 本地存档的SRA文件
       • 公共数据集
    
    📊 输出文件说明 | Output Files:
    
    \b
    单端测序输出:
       • SRR12345678.fastq.gz       # 压缩格式（推荐）
       • SRR12345678.fastq          # 非压缩格式
    
    \b
    双端测序输出(--split):
       • SRR12345678_1.fastq.gz     # Read 1（前端）
       • SRR12345678_2.fastq.gz     # Read 2（后端）
    
    \b
    双端测序输出(--no-split):
       • SRR12345678.fastq.gz       # 交错格式（interleaved）
    
    \b
    统计报告:
       • conversion_summary.txt     # 转换统计摘要
       • conversion.log             # 详细运行日志
       • failed_files.txt           # 失败文件列表
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本转换（自动加速）
    biopytools sra2fastq -i SRR12345678.sra -o fastq_output/
    
    \b
    # 📁 批量转换文件夹
    biopytools sra2fastq -i sra_files/ -o fastq_results/
    
    \b
    # 🚀 多线程高速转换
    biopytools sra2fastq -i input.sra -o output/ -t 64
    
    \b
    # 💾 使用SSD临时目录加速
    biopytools sra2fastq -i data.sra -o fastq/ \\
        --tmpdir /mnt/ssd/tmp -t 32
    
    \b
    # 📝 不压缩输出（便于查看）
    biopytools sra2fastq -i test.sra -o output/ --no-compress
    
    \b
    # 🔗 保持双端交错格式
    biopytools sra2fastq -i paired.sra -o output/ --no-split
    
    \b
    # 📏 过滤短reads
    biopytools sra2fastq -i input.sra -o filtered/ --min-len 50
    
    \b
    # ✂️  剪切adapters并过滤
    biopytools sra2fastq -i seqs.sra -o clean/ --clip --min-len 36
    
    \b
    # 🎨 完整配置示例
    biopytools sra2fastq -i sra_dir/ -o results/ \\
        -t 88 --tmpdir /fast/tmp --compress --split \\
        --min-len 50 --clip
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 🧬 公共数据库数据批量下载处理
    • 📊 转录组测序数据转换
    • 🔬 基因组重测序数据准备
    • 💻 RNA-seq/DNA-seq项目数据处理
    • 🗄️ 长期存储SRA数据激活使用
    • ⚡ 时间敏感的快速转换需求
    
    🚀 工具对比 | Tool Comparison:
    
    \b
    转换工具性能对比（1GB SRA文件）:
       fastq-dump:           ~20分钟  ⭐
       parallel-fastq-dump:  ~5分钟   ⭐⭐⭐⭐⭐
       fasterq-dump:         ~3分钟   ⭐⭐⭐⭐⭐
    
    \b
    工具选择建议:
       • parallel-fastq-dump: 推荐！易用且快速
       • fasterq-dump: 最快，需要更多配置
       • fastq-dump: 传统工具，较慢
    
    \b
    本工具特点:
       ✅ 自动选择parallel-fastq-dump
       ✅ 无需复杂配置，开箱即用
       ✅ 智能并行处理
       ✅ 兼容性好，稳定可靠
    
    ⚙️ 参数说明 | Parameter Details:
    
    \b
    线程数 (--threads):
       • 默认: 88（充分利用服务器）
       • 推荐设置:
         - 个人电脑: 4-8
         - 工作站: 16-32
         - 服务器: 32-88
       • 注意: parallel-fastq-dump会自动优化
    
    \b
    临时目录 (--tmpdir):
       • 默认: 当前目录
       • 推荐: 使用高速存储
         - SSD: 显著提升速度
         - /tmp: Linux系统临时目录
         - RAM disk: 极致性能
       • 空间需求: SRA文件的2-3倍
    
    \b
    压缩选项:
       • --compress (默认):
         ✅ 节省70-80%存储空间
         ✅ 大多数工具直接支持.gz
         ✅ 推荐长期存储使用
       
       • --no-compress:
         📝 便于直接查看和编辑
         💾 需要更多磁盘空间
         ⚡ 临时文件或小数据集
    
    \b
    拆分选项:
       • --split (默认):
         ✂️  双端数据分为_1和_2文件
         ✅ 符合大多数分析工具要求
         ✅ 便于质控和比对
       
       • --no-split:
         🔗 输出交错格式(interleaved)
         📦 某些工具特定需求
         💾 减少文件数量
    
    💾 存储空间需求 | Storage Requirements:
    
    \b
    空间计算:
       • SRA原始文件: X GB
       • 临时文件: 2-3X GB
       • 未压缩FASTQ: 3-4X GB
       • 压缩FASTQ.gz: 0.8-1.2X GB
    
    \b
    示例（5GB SRA文件）:
       1. SRA文件: 5GB
       2. 临时空间: 10-15GB
       3. 压缩输出: 4-6GB
       4. 峰值需求: ~20GB
    
    \b
    空间优化建议:
       ✅ 使用--compress压缩输出
       ✅ 设置--tmpdir到独立分区
       ✅ 转换后及时清理临时文件
       ✅ 批量转换时分批处理
    
    ⚡ 性能优化 | Performance Optimization:
    
    \b
    速度优化技巧:
       1️⃣ 使用SSD作为--tmpdir
          提速效果: 2-3倍
       
       2️⃣ 合理设置线程数
          推荐: CPU核心数×1.5
       
       3️⃣ 足够的内存
          建议: 至少8GB可用内存
       
       4️⃣ 避免网络存储
          本地存储比NFS快得多
    
    \b
    实测性能数据（32线程）:
       • 1GB SRA → ~3-5分钟
       • 5GB SRA → ~15-20分钟
       • 10GB SRA → ~30-40分钟
       • 批量10个文件 → 并行处理
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题及解决:
       1️⃣ "parallel-fastq-dump: command not found"
          解决:
          → pip install parallel-fastq-dump
          → conda install -c bioconda parallel-fastq-dump
       
       2️⃣ "No space left on device"
          解决:
          → 清理磁盘空间
          → 使用--tmpdir指定大空间目录
          → 启用--compress减少输出大小
       
       3️⃣ 转换速度慢
          解决:
          → 增加--threads参数
          → 使用--tmpdir到SSD
          → 检查CPU和I/O负载
       
       4️⃣ 输出文件损坏
          解决:
          → 检查磁盘空间充足
          → 验证SRA文件完整性
          → 重新下载SRA文件
       
       5️⃣ 压缩失败
          解决:
          → 检查gzip是否安装
          → 使用--no-compress跳过压缩
          → 手动压缩: gzip *.fastq
    
    \b
    性能检查清单:
       ✅ parallel-fastq-dump已安装
       ✅ 足够的磁盘空间（5-8倍SRA大小）
       ✅ --tmpdir指向快速存储
       ✅ 线程数设置合理（≤CPU核心×2）
       ✅ 内存充足（建议8GB+）
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据准备
       • 使用prefetch批量下载SRA
       • 验证下载完整性（vdb-validate）
       • 规划足够的存储空间
       • 组织好文件目录结构
    
    \b
    2️⃣ 转换策略
       • 小文件(<5GB): 直接转换
       • 大文件(>10GB): 分配更多线程
       • 批量处理: 使用shell脚本循环
       • 关键数据: 转换后验证完整性
    
    \b
    3️⃣ 性能配置
       # 高性能配置（服务器）
       biopytools sra2fastq -i input.sra -o output/ \\
           -t 64 --tmpdir /ssd/tmp
       
       # 标准配置（工作站）
       biopytools sra2fastq -i input.sra -o output/ \\
           -t 16 --tmpdir /tmp
       
       # 低资源配置（个人电脑）
       biopytools sra2fastq -i input.sra -o output/ \\
           -t 4
    
    \b
    4️⃣ 批量处理脚本
       # 批量转换目录下所有SRA文件
       for sra in sra_dir/*.sra; do
           biopytools sra2fastq -i "$sra" -o fastq_out/ \\
               -t 32 --tmpdir /fast/tmp
       done
       
       # 从列表文件批量下载并转换
       cat accession_list.txt | while read acc; do
           prefetch $acc
           biopytools sra2fastq -i ${acc}.sra -o results/
       done
    
    \b
    5️⃣ 质量控制
       • 转换后运行FastQC检查质量
       • 验证reads数量是否正确
       • 检查文件完整性
       • 对比原始SRA统计信息
    
    📚 相关资源 | Related Resources:
    
    \b
    工具和文档:
    • parallel-fastq-dump: https://github.com/rvalieris/parallel-fastq-dump
    • SRA Toolkit: https://github.com/ncbi/sra-tools
    • NCBI SRA: https://www.ncbi.nlm.nih.gov/sra
    • ENA Browser: https://www.ebi.ac.uk/ena
    
    \b
    下游分析工具:
    • FastQC: 质量控制和评估
    • fastp: 高速数据清理
    • Trimmomatic: Adapter修剪
    • MultiQC: 批量质控报告
    
    \b
    学习资源:
    • SRA数据下载教程
    • FASTQ格式说明
    • 测序数据质控流程
    • 比对分析入门指南
    
    ⚠️ 重要提示 | Important Notes:
    
    \b
    • 确保安装了parallel-fastq-dump (pip/conda)
    • 磁盘空间需要SRA文件的5-8倍
    • 临时文件转换后会自动清理
    • SSD作为临时目录可显著提速
    • 大文件建议在后台运行(nohup/screen)
    • 批量转换前先测试单个文件
    • 转换完成后验证数据完整性
    • 公共数据使用时注意引用规范
    
    🔐 数据安全 | Data Security:
    
    \b
    • 敏感数据注意访问权限设置
    • 不在共享服务器存储个人数据
    • 遵守数据使用协议和条款
    • 及时删除不需要的临时文件
    • 重要数据做好备份
    """
    
    # 🚀 懒加载
    sra_main = _lazy_import_sra_main()
    
    # 构建参数列表
    args = ['sra2fastq.py']
    args.extend(['-i', input])
    
    # 输出控制
    if output != './fastq_output':
        args.extend(['-o', output])
    
    # 性能参数
    if threads != 88:
        args.extend(['-t', str(threads)])
    if tmpdir:
        args.extend(['--tmpdir', tmpdir])
    
    # 转换参数
    if not compress:
        args.append('--no-compress')
    if not split_files:
        args.append('--no-split')
    
    # 过滤参数
    if min_len != 0:
        args.extend(['--min-len', str(min_len)])
    if clip:
        args.append('--clip')
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        sra_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 SRA转换被用户中断 | Conversion interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 SRA转换执行失败 | Conversion execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv