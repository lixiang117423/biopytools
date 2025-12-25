"""
🧬 VCF SNP/INDEL过滤命令 | VCF SNP/INDEL Filtering Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_vcf_filter_main():
    """懒加载VCF过滤main函数 | Lazy load VCF filter main function"""
    try:
        from ...filter_snp_indel.main import main as vcf_filter_main
        return vcf_filter_main
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


@click.command(
    short_help='🧬 VCF变异过滤工具：分离和过滤SNP/INDEL变异',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# ===== 必需参数 | Required Parameters =====
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 输入VCF文件路径 | Input VCF file path (supports .vcf/.vcf.gz)')

# ===== 输出控制 | Output Control =====
@click.option('--output', '-o',
              default='./filtered_vcf',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./filtered_vcf)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')

# ===== SNP过滤参数 | SNP Filtering Parameters =====
@click.option('--snp-qual',
              type=float,
              default=30.0,
              help='🔹 SNP最小质量值 | SNP minimum QUAL (default: 30.0)')
@click.option('--snp-dp',
              type=int,
              default=10,
              help='📊 SNP最小测序深度 | SNP minimum DP (default: 10)')
@click.option('--snp-mq',
              type=float,
              default=40.0,
              help='🎯 SNP最小比对质量 | SNP minimum MQ (default: 40.0)')
@click.option('--snp-qd',
              type=float,
              default=2.0,
              help='📈 SNP最小质量/深度比 | SNP minimum QD (default: 2.0)')
@click.option('--snp-fs',
              type=float,
              default=60.0,
              help='⚖️ SNP最大FisherStrand值 | SNP maximum FS (default: 60.0)')
@click.option('--snp-sor',
              type=float,
              default=3.0,
              help='📉 SNP最大StrandOddsRatio | SNP maximum SOR (default: 3.0)')
@click.option('--snp-mqrs',
              type=float,
              default=-12.5,
              help='🔍 SNP最小MappingQualityRankSum | SNP minimum MQRankSum (default: -12.5)')
@click.option('--snp-rprs',
              type=float,
              default=-8.0,
              help='📍 SNP最小ReadPosRankSum | SNP minimum ReadPosRankSum (default: -8.0)')
# [新增] 这里添加 MAF 参数
@click.option('--snp-maf',
              type=float,
              default=0.05,
              help='📉 SNP最小次等位基因频率 | SNP minimum MAF (default: 0.05)')

# ===== INDEL过滤参数 | INDEL Filtering Parameters =====
@click.option('--indel-qual',
              type=float,
              default=30.0,
              help='🔸 INDEL最小质量值 | INDEL minimum QUAL (default: 30.0)')
@click.option('--indel-dp',
              type=int,
              default=10,
              help='📊 INDEL最小测序深度 | INDEL minimum DP (default: 10)')
@click.option('--indel-mq',
              type=float,
              default=40.0,
              help='🎯 INDEL最小比对质量 | INDEL minimum MQ (default: 40.0)')
@click.option('--indel-qd',
              type=float,
              default=2.0,
              help='📈 INDEL最小质量/深度比 | INDEL minimum QD (default: 2.0)')
@click.option('--indel-fs',
              type=float,
              default=200.0,
              help='⚖️ INDEL最大FisherStrand值 | INDEL maximum FS (default: 200.0)')
@click.option('--indel-sor',
              type=float,
              default=10.0,
              help='📉 INDEL最大StrandOddsRatio | INDEL maximum SOR (default: 10.0)')
@click.option('--indel-rprs',
              type=float,
              default=-20.0,
              help='📍 INDEL最小ReadPosRankSum | INDEL minimum ReadPosRankSum (default: -20.0)')

# ===== 工具路径 | Tool Path =====
@click.option('--bcftools-path',
              default='bcftools',
              help='🔧 BCFtools软件路径 | BCFtools software path (default: bcftools)')
def filter_snp_indel(input, output, threads, 
               snp_qual, snp_dp, snp_mq, snp_qd, snp_fs, snp_sor, snp_mqrs, snp_rprs, snp_maf,
               indel_qual, indel_dp, indel_mq, indel_qd, indel_fs, indel_sor, indel_rprs,
               bcftools_path):
    """
    🧬 VCF SNP/INDEL过滤工具 | VCF SNP/INDEL Filtering Tool
    
    基于GATK最佳实践的VCF变异过滤工具，自动分离SNP和INDEL并应用
    不同的质量控制标准，生成高质量的变异集合。
    
    ✨ 功能特点 | Key Features:
    
    \b
    🔬 智能分离:
       • 自动识别SNP和INDEL
       • 使用BCFtools高效分离
       • 保留完整注释信息
       • 支持压缩格式输入
    
    \b
    🎯 精准过滤:
       • GATK最佳实践标准
       • SNP/INDEL独立过滤阈值
       • 多维度质量评估
       • 灵活的参数定制
    
    \b
    📊 完整流程:
       • 变异类型分离
       • 质量过滤处理
       • 过滤后变异合并
       • 详细统计报告
    
    🔄 分析流程 | Analysis Pipeline:
    
    \b
    1️⃣ 变异分离 | Variant Separation:
       输入VCF → 分离SNP和INDEL
       • filtered.raw.snp.vcf.gz
       • filtered.raw.indel.vcf.gz
    
    \b
    2️⃣ SNP过滤 | SNP Filtering:
       应用SNP特定过滤标准
       • QUAL ≥ 30
       • DP ≥ 10
       • MQ ≥ 40
       • QD ≥ 2.0
       • FS ≤ 60.0
       • SOR ≤ 3.0
       • MQRankSum ≥ -12.5
       • ReadPosRankSum ≥ -8.0
       • MAF ≥ 0.05            <-- [建议] 在文档字符串中也可以加上这一行说明
    
    \b
    3️⃣ INDEL过滤 | INDEL Filtering:
       应用INDEL特定过滤标准
       • QUAL ≥ 30
       • DP ≥ 10
       • MQ ≥ 40
       • QD ≥ 2.0
       • FS ≤ 200.0
       • SOR ≤ 10.0
       • ReadPosRankSum ≥ -20.0
    
    \b
    4️⃣ 变异合并 | Variant Merging:
       合并过滤后的SNP和INDEL
       • filtered.filtered.merged.vcf.gz
    
    \b
    5️⃣ 统计报告 | Statistics Report:
       生成详细的过滤统计
       • 变异数量统计
       • 过滤率计算
       • 质量分布分析
    
    📁 输入文件要求 | Input File Requirements:
    
    \b
    VCF文件格式:
       • 标准VCF格式(v4.0+)
       • 支持压缩(.vcf.gz)和未压缩(.vcf)
       • 必须包含必要的注释字段:
         - QUAL: 变异质量分数
         - DP: 测序深度
         - MQ: 比对质量
         - QD: 质量/深度比
         - FS: FisherStrand
         - SOR: StrandOddsRatio
         - MQRankSum: 比对质量秩和检验
         - ReadPosRankSum: 读段位置秩和检验
    
    \b
    VCF来源:
       • GATK HaplotypeCaller输出
       • GATK GenotypeGVCFs输出
       • 其他标准VCF calling工具
       • 确保包含完整的INFO字段
    
    📊 输出文件说明 | Output Files:
    
    \b
    主要输出文件:
       📄 filtered.raw.snp.vcf.gz
          → 原始SNP变异（分离后）
       
       📄 filtered.filtered.snp.vcf.gz
          → 过滤后的高质量SNP
       
       📄 filtered.raw.indel.vcf.gz
          → 原始INDEL变异（分离后）
       
       📄 filtered.filtered.indel.vcf.gz
          → 过滤后的高质量INDEL
       
       📄 filtered.filtered.merged.vcf.gz
          → 合并的高质量变异集合
       
       📊 filtering_statistics.txt
          → 详细的过滤统计报告
       
       📋 vcf_filtering.log
          → 完整的运行日志
    
    \b
    索引文件:
       • *.vcf.gz.tbi - 所有压缩VCF的索引
    
    💡 使用示例 | Usage Examples:
    
    \b
    # 🎯 基本用法（使用默认参数）
    biopytools filter-snp-indel -i variants.vcf -o filtered_output/
    
    \b
    # 🔧 自定义线程数
    biopytools filter-snp-indel -i variants.vcf.gz -o output/ -t 64
    
    \b
    # 🔹 严格的SNP过滤
    biopytools filter-snp-indel -i input.vcf -o strict_snp/ \\
        --snp-qual 40 --snp-dp 15 --snp-mq 50
    
    \b
    # 🔸 宽松的INDEL过滤
    biopytools filter-snp-indel -i input.vcf -o lenient_indel/ \\
        --indel-qual 25 --indel-fs 250
    
    \b
    # 🎨 完全自定义所有参数
    biopytools filter-snp-indel -i variants.vcf -o custom_filter/ \\
        --snp-qual 35 --snp-dp 12 --snp-mq 45 --snp-qd 2.5 \\
        --snp-fs 50 --snp-sor 2.5 --snp-mqrs -10 --snp-rprs -5 \\
        --indel-qual 35 --indel-dp 12 --indel-qd 2.5 \\
        --indel-fs 180 --indel-sor 8
    
    \b
    # 🏥 临床级别高严格过滤
    biopytools filter-snp-indel -i clinical.vcf -o clinical_filtered/ \\
        --snp-qual 50 --snp-dp 20 --snp-mq 60 --snp-qd 3.0 \\
        --indel-qual 50 --indel-dp 20 --indel-qd 3.0
    
    \b
    # 🔬 研究级别标准过滤
    biopytools filter-snp-indel -i research.vcf -o research_filtered/ \\
        --snp-qual 30 --snp-dp 10 --indel-qual 30 --indel-dp 10
    
    🎯 应用场景 | Use Cases:
    
    \b
    • 🧬 全基因组测序(WGS)变异过滤
    • 🔬 外显子组测序(WES)变异质控
    • 🏥 临床基因检测数据处理
    • 📊 群体遗传学变异筛选
    • 🌱 育种项目变异鉴定
    • 💊 药物基因组学研究
    • 🧪 功能基因组学分析
    
    ⚙️ 过滤参数详解 | Parameter Details:
    
    \b
    📊 通用质量参数:
       QUAL (Quality Score):
       • 变异调用的质量分数
       • Phred-scaled质量值
       • 建议: ≥30 (研究级), ≥40 (临床级)
       
       DP (Depth):
       • 位点的测序深度
       • 影响变异检出的可信度
       • 建议: ≥10 (最低), ≥20 (高质量)
       
       MQ (Mapping Quality):
       • 比对质量均值
       • 反映比对的可靠性
       • 建议: ≥40 (标准), ≥50 (严格)
       
       QD (Quality by Depth):
       • 质量分数/深度比值
       • 归一化的质量指标
       • 建议: ≥2.0
    
    \b
    🔹 SNP特异参数:
       FS (FisherStrand):
       • 链偏好性检验
       • 检测测序偏差
       • SNP建议: ≤60.0
       
       SOR (StrandOddsRatio):
       • 链比值比
       • 更敏感的链偏检测
       • SNP建议: ≤3.0
       
       MQRankSum:
       • 比对质量秩和检验
       • 比较参考/变异reads的MQ
       • SNP建议: ≥-12.5
       
       ReadPosRankSum:
       • 读段位置秩和检验
       • 检测位置偏差
       • SNP建议: ≥-8.0
    
    \b
    🔸 INDEL特异参数:
       FS (FisherStrand):
       • INDEL更宽松: ≤200.0
       • INDEL易受链偏影响
       
       SOR (StrandOddsRatio):
       • INDEL建议: ≤10.0
       • 比SNP更宽松
       
       ReadPosRankSum:
       • INDEL建议: ≥-20.0
       • INDEL位置偏差容忍度更高
    
    \b
    参数调整建议:
       🏥 临床诊断:
       → 提高所有阈值(QUAL≥50, DP≥20)
       → 降低FS和SOR限制
       → 追求高特异性
       
       🔬 科研发现:
       → 使用标准阈值
       → 平衡敏感性和特异性
       → 后续验证重要变异
       
       📊 群体研究:
       → 可适当放宽阈值
       → 注重样本间一致性
       → 结合频率信息
    
    📈 过滤效果评估 | Filtering Effect Evaluation:
    
    \b
    统计报告内容:
       • 📊 原始变异数量
         - 总变异数
         - SNP数量
         - INDEL数量
       
       • ✅ 通过过滤的变异
         - 高质量SNP
         - 高质量INDEL
         - 合并后总数
       
       • ❌ 被过滤掉的变异
         - 低质量SNP
         - 低质量INDEL
         - 过滤原因统计
       
       • 📉 过滤率计算
         - SNP过滤率
         - INDEL过滤率
         - 总体过滤率
    
    \b
    质量评估指标:
       • Ti/Tv比值(SNP)
         - 转换/颠换比
         - 全基因组预期: ~2.0-2.1
         - 外显子组预期: ~3.0-3.3
       
       • INDEL长度分布
         - 短INDEL(<10bp)占比
         - 长INDEL(≥10bp)占比
       
       • 质量分数分布
         - QUAL分布
         - DP分布
         - 其他指标分布
    
    🛠️ 故障排除 | Troubleshooting:
    
    \b
    常见问题:
       1️⃣ 过滤后变异过少
          → 检查原始VCF质量
          → 适当放宽过滤阈值
          → 确认测序深度充足
       
       2️⃣ 缺少必需字段
          → 确认VCF包含所有INFO字段
          → 使用GATK HaplotypeCaller
          → 检查VCF格式规范性
       
       3️⃣ Ti/Tv比值异常
          → 检查参考基因组版本
          → 评估样本污染可能
          → 调整过滤参数
       
       4️⃣ BCFtools错误
          → 确认BCFtools已安装
          → 检查版本兼容性(1.10+)
          → 指定--bcftools-path
    
    \b
    性能优化:
       • 使用压缩VCF输入(.vcf.gz)
       • 合理设置线程数
       • 使用SSD存储临时文件
       • 大文件考虑分染色体处理
    
    🏆 最佳实践 | Best Practices:
    
    \b
    1️⃣ 过滤前准备:
       • 确认VCF格式完整性
       • 检查必需注释字段
       • 了解数据来源和特点
       • 评估测序深度分布
    
    \b
    2️⃣ 参数选择策略:
       # 首次过滤使用标准参数
       biopytools filter-snp-indel -i input.vcf -o standard/
       
       # 评估结果后调整
       # 如果过滤率过高(>50%)，放宽阈值
       # 如果Ti/Tv异常，检查质量标准
    
    \b
    3️⃣ 结果验证:
       • 检查过滤统计报告
       • 评估Ti/Tv比值合理性
       • 验证已知变异位点保留率
       • 人工检查关键位点
    
    \b
    4️⃣ 下游分析建议:
       • 使用merged.vcf.gz进行注释
       • 分别分析SNP和INDEL效应
       • 结合功能预测工具
       • 考虑群体频率信息
    
    \b
    5️⃣ 质量控制流程:
       输入VCF
       → vcf-filter过滤
       → 评估统计报告
       → 必要时调整参数
       → 验证关键变异
       → 下游分析
    
    📚 相关资源 | Related Resources:
    
    \b
    • GATK Best Practices
    • BCFtools官方文档
    • VCF格式规范
    • 变异过滤标准
    • 质量控制指南
    
    \b
    推荐工具:
    • SnpEff/VEP: 变异注释
    • SnpSift: VCF过滤和统计
    • vcftools: VCF处理工具集
    • IGV: 可视化验证
    
    ⚠️ 注意事项 | Important Notes:
    
    \b
    • 过滤标准应根据具体研究调整
    • 临床应用需更严格的质量控制
    • 保留原始VCF以便重新分析
    • 不同测序平台可能需要不同阈值
    • INDEL过滤通常比SNP更宽松
    • 低深度区域的变异需谨慎解释
    • 建议结合多种质量指标综合判断
    • 过滤后应进行下游验证
    
    🔬 高级应用 | Advanced Applications:
    
    \b
    多样本过滤策略:
    • 先进行样本级别过滤
    • 再进行队列级别过滤
    • 结合allele frequency
    • 考虑家系信息
    
    \b
    特定区域过滤:
    • 外显子区更严格
    • 非编码区可放宽
    • 重复区域特殊处理
    • 结合区域功能重要性
    """
    
    # 🚀 懒加载
    vcf_filter_main = _lazy_import_vcf_filter_main()
    
    # 构建参数列表
    args = ['filter_snp_indel.py']
    args.extend(['-i', input])
    
    # 输出控制
    if output != './filtered_vcf':
        args.extend(['-o', output])
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # SNP过滤参数
    if snp_qual != 30.0:
        args.extend(['--snp-qual', str(snp_qual)])
    if snp_dp != 10:
        args.extend(['--snp-dp', str(snp_dp)])
    if snp_mq != 40.0:
        args.extend(['--snp-mq', str(snp_mq)])
    if snp_qd != 2.0:
        args.extend(['--snp-qd', str(snp_qd)])
    if snp_fs != 60.0:
        args.extend(['--snp-fs', str(snp_fs)])
    if snp_sor != 3.0:
        args.extend(['--snp-sor', str(snp_sor)])
    if snp_mqrs != -12.5:
        args.extend(['--snp-mqrs', str(snp_mqrs)])
    if snp_rprs != -8.0:
        args.extend(['--snp-rprs', str(snp_rprs)])
    # [新增] 处理 MAF 参数
    if snp_maf != 0.05:
        args.extend(['--snp-maf', str(snp_maf)])
    
    # INDEL过滤参数
    if indel_qual != 30.0:
        args.extend(['--indel-qual', str(indel_qual)])
    if indel_dp != 10:
        args.extend(['--indel-dp', str(indel_dp)])
    if indel_mq != 40.0:
        args.extend(['--indel-mq', str(indel_mq)])
    if indel_qd != 2.0:
        args.extend(['--indel-qd', str(indel_qd)])
    if indel_fs != 200.0:
        args.extend(['--indel-fs', str(indel_fs)])
    if indel_sor != 10.0:
        args.extend(['--indel-sor', str(indel_sor)])
    if indel_rprs != -20.0:
        args.extend(['--indel-rprs', str(indel_rprs)])
    
    # 工具路径
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
    
    # 执行主程序
    original_argv = sys.argv
    sys.argv = args
    
    try:
        vcf_filter_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 VCF过滤被用户中断 | Filtering interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 VCF过滤执行失败 | Filtering execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv