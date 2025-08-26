"""
Parabricks WGS分析命令 | Parabricks WGS Analysis Command
"""

import click
import sys
from ...parabricks.main import main as parabricks_main


@click.command(short_help='基于GPU加速的全基因组测序(WGS)批处理分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='输入目录路径 (包含clean FASTQ文件) | Input directory path (containing clean FASTQ files)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录路径 | Output directory path')
@click.option('--reference', '-r',
              required=True,
              type=click.Path(exists=True),
              help='参考基因组文件路径 | Reference genome file path')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='线程数 | Number of threads (default: 88)')
@click.option('--parabricks-path',
              default='/share/apps/parabricks/parabricks.CAT_2.2.1/bin/parabricks',
              type=str,
              help='parabricks程序路径 | parabricks program path (default: /share/apps/parabricks/parabricks.CAT_2.2.1/bin/parabricks)')
@click.option('--tmp-dir',
              type=click.Path(),
              help='临时目录路径 (默认使用输出目录下的tmp) | Temporary directory path (default: tmp under output directory)')
@click.option('--min-confidence',
              default=30,
              type=int,
              help='最小置信度阈值 | Minimum confidence threshold (default: 30)')
@click.option('--min-base-quality',
              default=20,
              type=int,
              help='最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
@click.option('--ploidy',
              default=2,
              type=int,
              help='倍性 | Ploidy (default: 2)')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              type=str,
              help='PCR indel模型 | PCR indel model (default: CONSERVATIVE)')
@click.option('--read1-pattern',
              default='*_1.clean.fq.gz',
              type=str,
              help='R1文件匹配模式 | R1 file pattern (default: *_1.clean.fq.gz)')
@click.option('--read2-pattern',
              default='*_2.clean.fq.gz',
              type=str,
              help='R2文件匹配模式 | R2 file pattern (default: *_2.clean.fq.gz)')
def parabricks(input_dir, output_dir, reference, threads, parabricks_path, tmp_dir,
                   min_confidence, min_base_quality, ploidy, pcr_indel_model,
                   read1_pattern, read2_pattern):
    """
    Parabricks WGS批处理分析脚本 (模块化版本)
    
    基于NVIDIA Parabricks的GPU加速全基因组测序(WGS)批处理分析工具，
    提供从FASTQ到VCF的完整分析流程，适用于大规模基因组数据的高效处理。
    
    功能特点 | Features:
    - 🚀 GPU加速的高性能计算
    - 📦 从FASTQ到VCF的完整流程
    - 🔄 智能批处理文件管理
    - 🎯 灵活的质量控制参数
    - 📊 详细的处理统计报告
    - 💾 高效的存储空间管理
    - 🛡️ 稳定的错误处理机制
    
    分析流程 | Analysis Pipeline:
    1. FASTQ文件自动发现和配对
    2. 读段质量评估和预处理
    3. GPU加速序列比对 (BWA-MEM)
    4. 读段排序和去重处理
    5. 碱基质量分数重校准
    6. GPU加速变异检测 (HaplotypeCaller)
    7. VCF文件生成和质量过滤
    8. 批处理统计报告生成
    
    应用场景 | Use Cases:
    - 🧬 大规模群体基因组学研究
    - 🏥 临床基因组学诊断流程
    - 🔬 基因组关联研究 (GWAS)
    - 🌱 农业基因组育种项目
    - 🧪 药物基因组学研究
    - 📈 进化基因组学分析
    - 💊 精准医疗数据处理
    
    示例 | Examples:
    
    \b
    # 基本WGS批处理分析
    biopytools parabricks -i /path/to/clean/data \\
        -o /path/to/output -r /path/to/reference.fa
    
    \b
    # 指定线程数和自定义parabricks路径
    biopytools parabricks -i ./clean_data -o ./parabricks_results \\
        -r ./reference.fa -t 64 --parabricks-path /custom/parabricks/path
    
    \b
    # 自定义质量控制参数
    biopytools parabricks -i /data/clean -o /results \\
        -r /genome/ref.fa --min-confidence 25 --min-base-quality 15
    
    \b
    # 指定临时目录和特殊文件模式
    biopytools parabricks -i ./fastq_files -o ./wgs_output \\
        -r ./genome.fa --tmp-dir /fast_ssd/tmp \\
        --read1-pattern "*_R1_*.fastq.gz" --read2-pattern "*_R2_*.fastq.gz"
    
    \b
    # 高通量批处理 (最大线程数)
    biopytools parabricks -i /large_cohort/fastq -o /large_cohort/results \\
        -r /reference/genome.fa -t 128 --min-confidence 35
    
    \b
    # 临床级别高精度分析
    biopytools parabricks -i /clinical/samples -o /clinical/results \\
        -r /clinical/reference.fa --min-confidence 40 \\
        --min-base-quality 25 --pcr-indel-model AGGRESSIVE
    
    \b
    # 单倍体样本分析
    biopytools parabricks -i /haploid/samples -o /haploid/results \\
        -r /reference.fa --ploidy 1 --min-confidence 30
    
    输入文件要求 | Input File Requirements:
    
    目录结构要求:
    - 输入目录应包含配对的FASTQ文件
    - 文件命名遵循标准配对模式
    - 支持压缩格式 (.fq.gz, .fastq.gz)
    - 建议预先进行质量控制和清理
    
    FASTQ文件命名规范:
    默认模式:
    - R1文件: sample_name_1.clean.fq.gz
    - R2文件: sample_name_2.clean.fq.gz
    
    自定义模式示例:
    - --read1-pattern "*_R1_*.fastq.gz"
    - --read2-pattern "*_R2_*.fastq.gz"
    - --read1-pattern "*.1.fq.gz"
    - --read2-pattern "*.2.fq.gz"
    
    参考基因组要求:
    - 标准FASTA格式
    - 建议使用BWA索引预处理
    - 常用版本: GRCh38, hg19, GRCh37
    - 确保与研究目标一致
    
    系统和硬件要求 | System & Hardware Requirements:
    
    GPU硬件要求:
    - NVIDIA GPU with CUDA Compute Capability 6.0+
    - 建议: Tesla V100, A100, RTX 系列
    - 显存: 最少8GB，建议16GB+
    - 多GPU支持可显著提升性能
    
    系统配置建议:
    - RAM: 至少64GB，大规模分析建议128GB+
    - 存储: NVMe SSD推荐，确保足够空间
    - CPU: 多核处理器支持并行IO操作
    - 网络: 高速存储访问对于大文件处理至关重要
    
    软件依赖:
    - NVIDIA GPU Driver (最新版本)
    - CUDA Toolkit (11.0+)
    - NVIDIA Parabricks (3.0+)
    - BWA索引工具 (可选，用于参考基因组预处理)
    
    参数详解 | Parameter Details:
    
    核心参数 | Core Parameters:
    --threads (-t):
    - 控制CPU线程数量
    - 影响内存使用和处理速度
    - 建议设置为CPU核心数的1-2倍
    - 大内存系统可设置更高值
    
    --tmp-dir:
    - 指定临时文件存储位置
    - 建议使用高速SSD存储
    - 确保有足够空间存储中间文件
    - 分析完成后会自动清理
    
    质量控制参数 | Quality Control Parameters:
    
    --min-confidence:
    - 变异检测的最小置信度
    - 范围: 10-50, 默认30
    - 更高值提供更高精度但可能丢失真实变异
    - 临床应用建议使用35-40
    
    --min-base-quality:
    - 碱基调用的最小质量分数
    - 范围: 10-30, 默认20
    - 影响变异检测的敏感性
    - 高通量研究可适当降低到15
    
    --ploidy:
    - 指定基因组倍性
    - 人类通常为2 (二倍体)
    - 某些研究可能需要1 (单倍体) 或其他值
    - 影响基因型调用策略
    
    --pcr-indel-model:
    - PCR扩增引起的indel处理模式
    - CONSERVATIVE: 保守模式，适用于高质量数据
    - AGGRESSIVE: 激进模式，适用于PCR偏差较大的数据
    - NONE: 不进行PCR indel校正
    
    文件模式配置 | File Pattern Configuration:
    
    --read1-pattern / --read2-pattern:
    - 使用shell通配符语法
    - 确保能正确匹配配对的FASTQ文件
    - 常用模式:
      * "*_1.clean.fq.gz" / "*_2.clean.fq.gz"
      * "*_R1_*.fastq.gz" / "*_R2_*.fastq.gz"
      * "*.1.fq.gz" / "*.2.fq.gz"
    
    输出文件结构 | Output File Structure:
    
    主要输出目录:
    - bam/: 存储排序后的BAM文件
      * {sample}.sorted.bam: 排序的比对文件
      * {sample}.sorted.bam.bai: BAM索引文件
    
    - vcf/: 存储变异检测结果
      * {sample}.vcf.gz: 压缩的VCF文件
      * {sample}.vcf.gz.tbi: VCF索引文件
    
    - tmp/: 临时文件目录
      * 中间处理文件
      * 分析完成后可以删除
    
    分析报告文件:
    - parabricks_summary.txt: 批处理总结报告
    - processing_statistics.txt: 处理统计信息
    - parabricks_analysis.log: 详细的分析日志
    - failed_samples.txt: 失败样本列表 (如果有)
    
    性能优化建议 | Performance Optimization:
    
    GPU优化:
    - 确保GPU内存充足避免频繁交换
    - 多GPU环境下parabricks会自动负载均衡
    - 监控GPU利用率确保资源充分使用
    - 避免CPU和GPU计算瓶颈不匹配
    
    存储优化:
    - 使用高速NVMe SSD作为临时目录
    - 确保输入输出路径在高性能存储上
    - 预估存储空间需求: 输入数据的3-5倍
    - 考虑使用网络存储时的带宽限制
    
    内存优化:
    - 线程数不要超过可用内存/样本大小的比值
    - 大基因组分析建议至少64GB内存
    - 监控内存使用避免系统交换
    - 必要时减少并行处理的样本数量
    
    质量控制策略 | Quality Control Strategy:
    
    数据预处理建议:
    - 使用高质量的clean FASTQ作为输入
    - 预先去除adapters和低质量读段
    - 验证FASTQ文件完整性和格式
    - 确保配对文件的读段数量一致
    
    参数调优指南:
    - 研究级分析: 默认参数通常足够
    - 临床级分析: 提高置信度和质量阈值
    - 高通量筛选: 可适当降低质量要求提高速度
    - 特殊样本: 根据样本特性调整倍性等参数
    
    错误诊断和解决 | Error Diagnosis & Solutions:
    
    常见错误类型:
    1. GPU内存不足:
       - 减少并行处理的样本数
       - 使用更大显存的GPU
       - 调整parabricks内存设置
    
    2. 存储空间不足:
       - 清理临时文件释放空间
       - 使用更大容量的存储设备
       - 分批处理减少同时占用的空间
    
    3. 文件格式错误:
       - 验证FASTQ文件格式和完整性
       - 检查文件命名是否符合模式
       - 确认参考基因组文件正确性
    
    4. 配对文件缺失:
       - 检查R1和R2文件配对完整性
       - 验证文件命名模式匹配
       - 确认文件权限和可访问性
    
    批处理管理 | Batch Processing Management:
    
    样本管理策略:
    - 按批次组织样本避免系统负载过大
    - 优先处理重要或时间敏感的样本
    - 保留失败样本列表便于重新处理
    - 定期备份重要的中间结果
    
    监控和日志:
    - 实时监控GPU和CPU使用率
    - 跟踪每个样本的处理进度
    - 记录详细的错误信息便于调试
    - 统计处理时间优化工作流程
    
    结果验证建议:
    - 检查输出BAM文件的比对统计
    - 验证VCF文件的变异质量分布
    - 比较不同样本间的统计指标
    - 使用已知变异位点验证检测准确性
    
    高级配置技巧 | Advanced Configuration Tips:
    
    大规模数据处理:
    - 使用分布式计算集群
    - 配置网络文件系统提高访问效率
    - 实施数据管道自动化处理
    - 建立质量监控和报警系统
    
    特殊样本处理:
    - 肿瘤样本: 调整变异检测敏感性
    - 古DNA样本: 降低质量要求增加覆盖度
    - 近亲繁殖样本: 注意杂合子缺失问题
    - 混合样本: 使用特殊的基因型调用策略
    
    流程集成建议:
    - 与上游质控流程无缝对接
    - 为下游分析准备标准化输出
    - 建立自动化的质量检查点
    - 集成现有的数据管理系统
    
    成本效益优化 | Cost-Benefit Optimization:
    
    计算资源规划:
    - 根据样本数量规划GPU使用时间
    - 平衡处理速度与资源成本
    - 考虑云计算vs本地计算的经济性
    - 优化批处理大小减少空闲时间
    
    数据管理策略:
    - 及时清理不需要的中间文件
    - 压缩长期存储的文件减少空间成本
    - 建立数据生命周期管理策略
    - 考虑冷热数据的分层存储
    
    技术规范 | Technical Specifications:
    
    支持的文件格式:
    - 输入: FASTQ (.fq, .fastq), 压缩格式 (.gz, .bz2)
    - 输出: BAM (sorted, indexed), VCF (compressed, indexed)
    - 参考: FASTA (.fa, .fasta), 支持多种索引格式
    
    处理能力范围:
    - 样本数量: 1-1000+ (取决于系统配置)
    - 基因组大小: 支持从细菌到哺乳动物基因组
    - 数据通量: TB级数据批处理能力
    - 并发处理: 根据GPU数量自动扩展
    
    兼容性说明:
    - 兼容GATK最佳实践流程
    - 输出格式符合标准规范
    - 支持多种下游分析工具
    - 与现有生信流程无缝集成
    """
    
    # 构建参数列表传递给原始main函数
    args = ['parabricks.py']
    
    # 必需参数
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])
    
    # 可选参数（只在非默认值时添加）
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if parabricks_path != '/share/apps/parabricks/parabricks.CAT_2.2.1/bin/parabricks':
        args.extend(['--parabricks-path', parabricks_path])
    
    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])
    
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])
    
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    
    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])
    
    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])
    
    if read1_pattern != '*_1.clean.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
    
    if read2_pattern != '*_2.clean.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        parabricks_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\nParabricks WGS分析被用户中断 | Parabricks WGS analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Parabricks WGS分析失败 | Parabricks WGS analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv