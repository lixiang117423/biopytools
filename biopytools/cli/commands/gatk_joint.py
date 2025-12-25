"""
🧬 GATK Joint Genotyping命令 | GATK Joint Genotyping Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_joint_main():
    """懒加载Joint Genotyping main函数 | Lazy load joint genotyping main function"""
    try:
        from ...gatk_joint.main import main as joint_main
        return joint_main
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
    short_help='🧬 GATK Joint Genotyping：多样本联合变异检测和过滤',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='📂 输入目录 | Input directory (containing VCF/GVCF files)')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='🧬 参考基因组文件 | Reference genome file (.fasta/.fa)')

# ===== 输出控制 | Output Control =====
@click.option('--output', '-o',
              default='./joint_genotyping_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./joint_genotyping_output)')

# ===== 计算资源 | Computing Resources =====
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='⚙️  线程数 | Number of threads (default: 88)')
@click.option('--memory', '-m',
              default='100g',
              help='💾 Java内存设置 | Java memory setting (default: 100g)')

# ===== 区间设置 | Interval Settings =====
@click.option('--intervals', '-L',
              help='📍 分析区间 | Analysis intervals (chromosome or interval file)')

# ===== SNP过滤参数 | SNP Filtering Parameters =====
@click.option('--snp-qd',
              type=float,
              default=2.0,
              help='🔹 SNP QD阈值 | SNP QD threshold (default: 2.0)')
@click.option('--snp-fs',
              type=float,
              default=60.0,
              help='⚖️  SNP FS阈值 | SNP FS threshold (default: 60.0)')
@click.option('--snp-mq',
              type=float,
              default=40.0,
              help='🎯 SNP MQ阈值 | SNP MQ threshold (default: 40.0)')
@click.option('--snp-mqrs',
              type=float,
              default=-12.5,
              help='🔍 SNP MQRankSum阈值 | SNP MQRankSum threshold (default: -12.5)')
@click.option('--snp-rprs',
              type=float,
              default=-8.0,
              help='📍 SNP ReadPosRankSum阈值 | SNP ReadPosRankSum threshold (default: -8.0)')
@click.option('--snp-sor',
              type=float,
              default=3.0,
              help='📉 SNP SOR阈值 | SNP SOR threshold (default: 3.0)')

# ===== INDEL过滤参数 | INDEL Filtering Parameters =====
@click.option('--indel-qd',
              type=float,
              default=2.0,
              help='🔸 INDEL QD阈值 | INDEL QD threshold (default: 2.0)')
@click.option('--indel-fs',
              type=float,
              default=200.0,
              help='⚖️  INDEL FS阈值 | INDEL FS threshold (default: 200.0)')
@click.option('--indel-rprs',
              type=float,
              default=-20.0,
              help='📍 INDEL ReadPosRankSum阈值 | INDEL ReadPosRankSum threshold (default: -20.0)')
@click.option('--indel-sor',
              type=float,
              default=10.0,
              help='📉 INDEL SOR阈值 | INDEL SOR threshold (default: 10.0)')

# ===== 工具路径 | Tool Paths =====
@click.option('--gatk-path',
              default='gatk',
              help='🔧 GATK软件路径 | GATK software path (default: gatk)')
@click.option('--bcftools-path',
              default='bcftools',
              help='🔧 BCFtools软件路径 | BCFtools software path (default: bcftools)')
def gatk_joint(input_dir, reference, output, threads, memory, intervals,
                     snp_qd, snp_fs, snp_mq, snp_mqrs, snp_rprs, snp_sor,
                     indel_qd, indel_fs, indel_rprs, indel_sor,
                     gatk_path, bcftools_path):
    """
    🧬 GATK Joint Genotyping工具 | GATK Joint Genotyping Tool
    
    基于GATK最佳实践的多样本联合变异检测工具，支持从多个GVCF文件进行
    联合分型（Joint Calling），提高变异检测的准确性和一致性。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🔗 Joint Genotyping优势:
       • 多样本联合分型
       • 提高低频变异检测准确性
       • 统一的变异质量评分
       • 减少批次效应
       • 支持增量添加样本
    
    \b
    📦 完整分析流程:
       • 自动检测GVCF/VCF文件
       • GenomicsDB高效存储
       • GATK GenotypeGVCFs联合分型
       • SNP/INDEL分离和过滤
       • 高质量变异集合生成
    
    \b
    🎯 精准过滤:
       • GATK最佳实践标准
       • SNP/INDEL独立过滤阈值
       • 多维度质量评估
       • 灵活的参数定制
    
    🔄 分析流程 | Analysis Pipeline:
    
    \b
    步骤1: 文件检测 | File Detection
       • 自动检测输入目录中的GVCF/VCF文件
       • 支持.g.vcf.gz和.vcf.gz格式
       • 创建样本映射文件
    
    \b
    步骤2: GenomicsDB导入 | GenomicsDB Import
       • 使用GATK GenomicsDBImport
       • 高效存储多样本数据
       • 支持区间并行处理
       • 优化内存使用
    
    \b
    步骤3: 联合分型 | Joint Genotyping
       • 使用GATK GenotypeGVCFs
       • 多样本联合变异检测
       • 统一质量评分
       • 生成联合VCF文件
    
    \b
    步骤4: SNP提取和过滤 | SNP Extraction & Filtering
       • SelectVariants提取SNP
       • VariantFiltration应用过滤条件
       • GATK最佳实践标准
       • 生成高质量SNP集合
    
    \b
    步骤5: INDEL提取和过滤 | INDEL Extraction & Filtering
       • SelectVariants提取INDEL
       • VariantFiltration应用过滤条件
       • INDEL特定过滤标准
       • 生成高质量INDEL集合
    
    \b
    步骤6: 变异合并 | Variant Merging
       • 合并过滤后的SNP和INDEL
       • 生成最终高质量变异集合
       • 创建索引文件
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    GVCF/VCF文件格式:
       • 标准GVCF格式: *.g.vcf.gz
       • 标准VCF格式: *.vcf.gz
       • 必须包含索引文件: *.tbi
       • 所有文件使用相同参考基因组
    
    \b
    文件来源:
       • GATK HaplotypeCaller输出的GVCF
       • Parabricks haplotypecaller输出
       • 其他符合GATK标准的GVCF
    
    \b
    目录结构示例:
       input_dir/
       ├── sample1.g.vcf.gz
       ├── sample1.g.vcf.gz.tbi
       ├── sample2.g.vcf.gz
       ├── sample2.g.vcf.gz.tbi
       └── sample3.g.vcf.gz
           └── sample3.g.vcf.gz.tbi
    
    \b
    参考基因组要求:
       • FASTA格式(.fa/.fasta)
       • 必须有索引文件(.fai)
       • 必须有字典文件(.dict)
       • 版本与GVCF一致
    
    📊 输出文件结构 | Output File Structure:
    
    \b
    主要输出文件:
       📂 output_dir/
       ├── genomicsdb/              # GenomicsDB数据库
       ├── joint_genotyped.vcf.gz   # 联合分型原始结果
       ├── raw.snp.vcf.gz           # 原始SNP
       ├── filtered.snp.vcf.gz      # 过滤后SNP
       ├── raw.indel.vcf.gz         # 原始INDEL
       ├── filtered.indel.vcf.gz    # 过滤后INDEL
       ├── filtered.merged.vcf.gz   # 合并的高质量变异
       ├── sample_map.txt           # 样本映射文件
       ├── summary_report.txt       # 统计报告
       └── joint_genotyping.log     # 运行日志
    
    \b
    索引文件:
       • 所有.vcf.gz文件都有对应的.tbi索引
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本用法（使用默认参数）
    biopytools gatk-joint -i gvcf_folder/ -r ref.fasta -o results/
    
    \b
    # ⚙️ 自定义计算资源
    biopytools gatk-joint -i vcf_dir/ -r genome.fa -o output/ \\
        -t 64 -m 200g
    
    \b
    # 📍 指定染色体区间
    biopytools gatk-joint -i data/ -r ref.fa -o out/ \\
        -L chr1 -t 32
    
    \b
    # 🔹 严格的SNP过滤
    biopytools gatk-joint -i gvcf/ -r ref.fa -o strict/ \\
        --snp-qd 3.0 --snp-fs 50.0 --snp-mq 50.0
    
    \b
    # 🔸 自定义INDEL过滤
    biopytools gatk-joint -i input/ -r ref.fa -o output/ \\
        --indel-qd 3.0 --indel-fs 180.0
    
    \b
    # 🎨 完全自定义过滤参数
    biopytools gatk-joint -i gvcf_dir/ -r genome.fa -o custom/ \\
        --snp-qd 2.5 --snp-fs 55.0 --snp-mq 45.0 \\
        --snp-mqrs -10.0 --snp-rprs -5.0 --snp-sor 2.5 \\
        --indel-qd 2.5 --indel-fs 180.0 \\
        --indel-rprs -15.0 --indel-sor 8.0
    
    \b
    # 🏥 临床级别高严格过滤
    biopytools gatk-joint -i clinical/ -r ref.fa -o clinical_out/ \\
        --snp-qd 4.0 --snp-mq 60.0 --indel-qd 4.0 -m 300g
    
    \b
    # 📊 大规模队列分析
    biopytools gatk-joint -i cohort_gvcf/ -r ref.fa -o cohort_results/ \\
        -t 128 -m 500g --intervals intervals.bed
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 🧬 群体基因组学研究
    • 👥 大规模队列测序项目
    • 🏥 临床基因组学诊断
    • 🔬 疾病关联研究(GWAS)
    • 🌱 育种项目变异鉴定
    • 📊 进化基因组学分析
    • 💊 药物基因组学研究
    
    ⚙️ 参数详解 | Parameter Details:
    
    \b
    💾 计算资源参数:
       --threads (-t):
       • 控制并行处理的线程数
       • 建议: CPU核心数×0.8-1.0
       • GenomicsDB导入主要使用
       • 过多线程可能导致内存压力
       
       --memory (-m):
       • GATK Java虚拟机内存限制
       • 格式: 数字+单位(g/G/m/M)
       • 建议: 样本数×5-10GB
       • 大队列(>100样本)需200GB+
    
    \b
    📍 区间设置:
       --intervals (-L):
       • 限制分析的基因组区间
       • 格式1: 染色体名(chr1, chr2)
       • 格式2: 区间文件(.bed, .list)
       • 用途: 
         - 分染色体并行处理
         - 仅分析外显子区域
         - 测试小区间
    
    \b
    🔹 SNP过滤参数详解:
       QD (Quality by Depth):
       • 质量分数/深度比值
       • 归一化的质量指标
       • 默认: ≥2.0
       • 严格: ≥3.0
       
       FS (FisherStrand):
       • 链偏好性检验
       • 检测测序偏差
       • 默认: ≤60.0
       • 严格: ≤50.0
       
       MQ (Mapping Quality):
       • 比对质量均值
       • 反映比对可靠性
       • 默认: ≥40.0
       • 严格: ≥50.0
       
       MQRankSum:
       • 比对质量秩和检验
       • 比较参考/变异reads
       • 默认: ≥-12.5
       • 严格: ≥-10.0
       
       ReadPosRankSum:
       • 读段位置秩和检验
       • 检测位置偏差
       • 默认: ≥-8.0
       • 严格: ≥-5.0
       
       SOR (StrandOddsRatio):
       • 链比值比
       • 更敏感的链偏检测
       • 默认: ≤3.0
       • 严格: ≤2.5
    
    \b
    🔸 INDEL过滤参数详解:
       INDEL过滤通常比SNP更宽松:
       
       QD: ≥2.0 (与SNP相同)
       FS: ≤200.0 (SNP的3倍多)
       ReadPosRankSum: ≥-20.0 (更宽松)
       SOR: ≤10.0 (SNP的3倍多)
       
       原因: INDEL更易受技术偏差影响
    
    🔗 Joint Genotyping vs 单样本分析 | Comparison:
    
    \b
    Joint Genotyping优势:
       ✅ 提高低频变异检测准确性
       ✅ 统一的变异质量评分标准
       ✅ 减少批次间的差异
       ✅ 更准确的基因型推断
       ✅ 支持增量添加新样本
       ✅ 提供群体水平的变异信息
    
    \b
    单样本分析特点:
       • 处理速度快
       • 适合临床快速诊断
       • 独立样本分析
       • 不需要其他样本
    
    \b
    使用建议:
       • 队列研究(>10样本) → Joint Genotyping
       • 群体基因组学 → Joint Genotyping
       • 单样本临床诊断 → 单样本分析
       • 计划扩展样本 → Joint Genotyping
    
    💾 存储空间需求 | Storage Requirements:
    
    \b
    空间计算:
       • GenomicsDB: GVCF总大小×1.5-2倍
       • 联合VCF: 样本数×单样本VCF大小×0.3
       • 过滤后VCF: 联合VCF×0.5-0.8
       • 临时文件: 上述总和×0.3
    
    \b
    示例(100样本, 每个GVCF 2GB):
       • GVCF总大小: 200GB
       • GenomicsDB: 300-400GB
       • 联合VCF: 60GB
       • 过滤后文件: 30GB
       • 峰值需求: 600-700GB
    
    ⚡ 性能优化建议 | Performance Tips:
    
    \b
    🚀 提升速度:
       • 使用SSD存储GenomicsDB
       • 增加内存减少磁盘I/O
       • 使用--intervals分染色体并行
       • 合理设置线程数
    
    \b
    💾 节省内存:
       • 减少并行线程数
       • 使用分批处理策略
       • 设置合理的Java内存上限
       • 定期清理临时文件
    
    \b
    📊 大规模优化:
       样本数>500时建议:
       • 分染色体并行处理
       • 每个染色体独立运行
       • 最后合并所有染色体
       • 使用高性能计算集群
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1️⃣ 内存不足错误
          → 增加--memory参数
          → 减少--threads数量
          → 使用--intervals分批处理
       
       2️⃣ GenomicsDB导入失败
          → 检查所有GVCF文件完整性
          → 验证索引文件(.tbi)存在
          → 确认参考基因组一致
       
       3️⃣ 联合分型速度慢
          → 使用SSD存储
          → 增加内存设置
          → 考虑分染色体处理
       
       4️⃣ 过滤后变异过少
          → 检查原始变异质量
          → 适当放宽过滤阈值
          → 评估测序深度
       
       5️⃣ 参考基因组错误
          → 确保有.fai索引文件
          → 确保有.dict字典文件
          → 版本与GVCF匹配
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 数据准备:
       • 确保所有GVCF使用相同参考基因组
       • 验证GVCF文件完整性
       • 准备好索引文件
       • 检查样本命名一致性
    
    \b
    2️⃣ 资源规划:
       # 小队列(<50样本)
       -t 32 -m 100g
       
       # 中等队列(50-200样本)
       -t 64 -m 300g
       
       # 大队列(>200样本)
       -t 128 -m 500g --intervals chr1
    
    \b
    3️⃣ 分批处理策略:
       # 分染色体处理
       for chr in chr{1..22} chrX chrY; do
           biopytools gatk-joint -i gvcf/ -r ref.fa \\
               -o results_${chr} -L ${chr} -t 64 -m 200g
       done
       
       # 合并所有染色体结果
       bcftools concat -o final.vcf.gz results_*/filtered.merged.vcf.gz
    
    \b
    4️⃣ 质量控制:
       • 检查联合VCF的样本数
       • 评估Ti/Tv比值
       • 验证已知变异位点
       • 检查过滤统计报告
    
    \b
    5️⃣ 增量分析:
       # 初始队列
       biopytools gatk-joint -i batch1/ -r ref.fa -o cohort1/
       
       # 添加新样本后重新分析
       # 将新GVCF复制到同一目录
       cp new_samples/*.g.vcf.gz batch1/
       
       # 重新运行Joint Genotyping
       biopytools gatk-joint -i batch1/ -r ref.fa -o cohort2/
    
    📚 相关资源 | Related Resources:
    
    \b
    • GATK官网: https://gatk.broadinstitute.org/
    • GATK最佳实践
    • GenomicsDB文档
    • Joint Genotyping教程
    • VCF格式规范
    
    \b
    推荐工具:
    • SnpEff/VEP: 变异注释
    • vcftools: VCF处理
    • IGV: 可视化验证
    • Plink: 群体遗传学分析
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 所有GVCF必须使用相同参考基因组
    • GenomicsDB需要大量磁盘空间
    • 大队列建议分染色体并行处理
    • 内存设置应考虑样本数量
    • 过滤标准应根据研究目的调整
    • 保留原始联合VCF以便重新过滤
    • 定期清理GenomicsDB临时文件
    • 建议在测试数据上先验证参数
    """
    
    # 🚀 懒加载
    joint_main = _lazy_import_joint_main()
    
    # 构建参数列表
    args = ['gatk_joint.py']
    args.extend(['-i', input_dir])
    args.extend(['-r', reference])
    
    # 输出控制
    if output != './joint_genotyping_output':
        args.extend(['-o', output])
    
    # 计算资源
    if threads != 88:
        args.extend(['-t', str(threads)])
    if memory != '100g':
        args.extend(['-m', memory])
    
    # 区间设置
    if intervals:
        args.extend(['-L', intervals])
    
    # SNP过滤参数
    if snp_qd != 2.0:
        args.extend(['--snp-qd', str(snp_qd)])
    if snp_fs != 60.0:
        args.extend(['--snp-fs', str(snp_fs)])
    if snp_mq != 40.0:
        args.extend(['--snp-mq', str(snp_mq)])
    if snp_mqrs != -12.5:
        args.extend(['--snp-mqrs', str(snp_mqrs)])
    if snp_rprs != -8.0:
        args.extend(['--snp-rprs', str(snp_rprs)])
    if snp_sor != 3.0:
        args.extend(['--snp-sor', str(snp_sor)])
    
    # INDEL过滤参数
    if indel_qd != 2.0:
        args.extend(['--indel-qd', str(indel_qd)])
    if indel_fs != 200.0:
        args.extend(['--indel-fs', str(indel_fs)])
    if indel_rprs != -20.0:
        args.extend(['--indel-rprs', str(indel_rprs)])
    if indel_sor != 10.0:
        args.extend(['--indel-sor', str(indel_sor)])
    
    # 工具路径
    if gatk_path != 'gatk':
        args.extend(['--gatk-path', gatk_path])
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        joint_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 Joint Genotyping被用户中断 | Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 Joint Genotyping执行失败 | Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv