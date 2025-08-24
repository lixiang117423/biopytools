"""
VCF PCA分析命令 | VCF PCA Analysis Command
"""

import click
import sys
from ...vcf_pca.main import main as vcf_pca_main


@click.command(short_help='VCF文件主成分分析(PCA)工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径 (支持压缩和未压缩) | Input VCF file path (supports compressed and uncompressed)')
@click.option('--output', '-o',
              default='./pca_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./pca_output)')
@click.option('--sample-info', '-s',
              type=click.Path(exists=True),
              help='样本信息文件 (制表符分隔) | Sample information file (tab-separated)')
@click.option('--components', '-c',
              default=10,
              type=int,
              help='主成分数量 | Number of principal components (default: 10)')
@click.option('--maf', '-m',
              default=0.05,
              type=float,
              help='最小等位基因频率阈值 | Minor allele frequency threshold (default: 0.05)')
@click.option('--missing',
              default=0.1,
              type=float,
              help='最大缺失率阈值 | Maximum missing rate threshold (default: 0.1)')
@click.option('--hwe',
              default=1e-6,
              type=float,
              help='Hardy-Weinberg平衡p值阈值 | Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)')
@click.option('--skip-qc',
              is_flag=True,
              help='跳过质量控制过滤，直接使用输入VCF文件 | Skip quality control filtering, use input VCF file directly')
@click.option('--plot', '-p',
              is_flag=True,
              help='生成PCA可视化图表 | Generate PCA visualization plots')
@click.option('--group-column', '-g',
              type=str,
              help='样本信息文件中用于分组的列名 | Column name for grouping in sample info file')
@click.option('--plink-path',
              default='plink',
              type=str,
              help='PLINK软件路径 | PLINK software path (default: plink)')
@click.option('--bcftools-path',
              default='bcftools',
              type=str,
              help='BCFtools软件路径 | BCFtools software path (default: bcftools)')
def vcf_pca(vcf, output, sample_info, components, maf, missing, hwe, skip_qc,
            plot, group_column, plink_path, bcftools_path):
    """
    VCF文件主成分分析(PCA)脚本 (模块化版本)
    
    一个完整的VCF基因型数据主成分分析流程工具，支持质量控制、PCA分析、
    样本信息整合和结果可视化，适用于群体遗传学和基因组学研究中的样本结构分析。
    
    功能特点 | Features:
    - 全自动VCF到PLINK格式转换
    - 多维度质量控制过滤
    - 高效主成分分析计算
    - 丰富的可视化图表生成
    - 样本群体信息整合
    - 详细的分析日志记录
    - 模块化可扩展架构
    
    分析流程 | Analysis Pipeline:
    1. VCF格式验证和转换
    2. 基因型质量控制过滤
    3. 主成分分析计算
    4. 结果可视化生成
    5. 样本信息整合
    6. 综合报告输出
    
    应用场景 | Use Cases:
    - 群体遗传结构分析
    - 样本质量评估和离群检测
    - 基因组关联研究(GWAS)前处理
    - 遗传多样性评估
    - 系统发育和进化分析
    - 育种和选择研究
    
    示例 | Examples:
    
    \b
    # 基本PCA分析
    biopytools vcf-pca -v variants.vcf -o pca_results
    
    \b
    # 指定主成分数量并生成可视化
    biopytools vcf-pca -v data.vcf.gz -o results -c 15 -p
    
    \b
    # 包含样本信息和群体分组
    biopytools vcf-pca -v variants.vcf -o pca_out \\
        -s samples.txt -g population -p
    
    \b
    # 跳过质控直接分析
    biopytools vcf-pca -v filtered_variants.vcf -o results --skip-qc -p
    
    \b
    # 自定义质控参数
    biopytools vcf-pca -v data.vcf -o results \\
        -m 0.01 --missing 0.05 --hwe 1e-5 -c 20 -p
    
    \b
    # 指定外部工具路径
    biopytools vcf-pca -v variants.vcf -o results \\
        --plink-path /usr/local/bin/plink \\
        --bcftools-path /opt/bcftools/bin/bcftools -p
    
    \b
    # 完整分析流程示例
    biopytools vcf-pca -v population_variants.vcf.gz \\
        -o comprehensive_pca \\
        -s sample_metadata.txt \\
        -g geographic_origin \\
        -c 20 -m 0.02 --missing 0.08 \\
        --plot
    
    输入文件格式 | Input File Formats:
    
    VCF文件要求:
    - 标准VCF 4.0+格式
    - 支持压缩文件(.vcf.gz)
    - 包含完整的基因型信息
    - 建议预先进行基本质量过滤
    
    样本信息文件格式 (可选):
    - 制表符分隔的文本文件
    - 第一列必须是样本ID (与VCF文件中样本名匹配)
    - 可包含多个注释列 (如群体、地理来源等)
    - 支持分组可视化分析
    
    样本信息文件示例:
    Sample_ID    Population    Geographic_Origin    Sex
    Sample001    EUR          Europe              M
    Sample002    EAS          East_Asia           F
    Sample003    AFR          Africa              M
    Sample004    EUR          Europe              F
    
    质量控制参数说明 | Quality Control Parameters:
    
    --maf (Minor Allele Frequency):
    - 过滤低频变异位点
    - 默认0.05 (5%)
    - 建议范围: 0.01-0.1
    - 用于去除可能的测序错误
    
    --missing (Missing Rate):
    - 过滤高缺失率位点
    - 默认0.1 (10%)
    - 建议范围: 0.05-0.2
    - 确保数据完整性
    
    --hwe (Hardy-Weinberg Equilibrium):
    - 过滤偏离HWE的位点
    - 默认1e-6
    - 建议范围: 1e-6 到 1e-3
    - 检测基因型错误和群体分层
    
    主成分分析参数 | PCA Parameters:
    
    --components:
    - 计算的主成分数量
    - 默认10个主成分
    - 建议范围: 5-50
    - 更多主成分提供更详细的结构信息
    
    可视化选项 | Visualization Options:
    
    生成的图表包括:
    - 碎石图 (Scree Plot): 显示各主成分解释的方差比例
    - PCA散点图: PC1 vs PC2, PC1 vs PC3, PC2 vs PC3
    - 成对散点图矩阵: 多个主成分的两两比较
    - 样本分组着色 (如果提供样本信息)
    
    输出文件说明 | Output Files:
    
    核心结果文件:
    - pca_eigenvalues.txt: 特征值文件
    - pca_eigenvectors.txt: 特征向量文件
    - pca_eigenvalues_formatted.txt: 格式化的特征值
    - pca_eigenvectors_formatted.txt: 格式化的特征向量
    
    整合结果文件 (如果提供样本信息):
    - pca_with_sample_info.txt: PCA结果与样本信息合并
    
    可视化文件 (如果启用--plot):
    - pca_scree_plot.png: 碎石图
    - pca_scatter_plots.png: PCA散点图
    - pca_pairs_plot.png: PCA成对图
    
    中间文件:
    - vcf_pca.bed/bim/fam: PLINK二进制格式文件
    - vcf_pca_filtered.*: 质量控制后的PLINK文件
    - analysis.log: 详细的分析日志
    - summary_report.txt: 分析总结报告
    
    性能和系统要求 | Performance & System Requirements:
    
    依赖软件:
    - PLINK (v1.9+): 用于基因型数据处理
    - BCFtools: 用于VCF文件操作
    - Python packages: pandas, numpy, matplotlib, seaborn
    
    系统建议:
    - RAM: 至少4GB，大文件建议16GB+
    - 存储: 输出目录需要足够空间存储中间文件
    - CPU: 多核处理器可提高处理速度
    
    文件大小估算:
    - 10K样本，100K位点: ~2GB内存，~5GB存储
    - 1K样本，1M位点: ~4GB内存，~10GB存储
    
    故障排除 | Troubleshooting:
    
    常见问题:
    1. "PLINK not found": 检查PLINK安装和路径设置
    2. "Memory error": 增加系统内存或减少数据规模
    3. "Empty result": 检查质控参数是否过严
    4. "Sample mismatch": 确保样本信息文件与VCF样本名匹配
    
    质量控制建议:
    - 大规模数据建议先进行预过滤
    - 检查样本和位点的质量分布
    - 适当调整质控参数避免过度过滤
    - 验证PCA结果的生物学合理性
    
    最佳实践 | Best Practices:
    
    1. 数据准备:
       - 使用高质量的VCF文件作为输入
       - 预先去除明显的低质量样本和位点
       - 确保VCF文件格式正确且完整
    
    2. 参数选择:
       - 根据研究目的调整质控参数严格程度
       - 对于群体结构分析，可适当放宽MAF阈值
       - 对于关联分析预处理，建议严格质控
    
    3. 结果解读:
       - 结合碎石图确定有意义的主成分数量
       - 注意PC1和PC2通常解释最多的遗传变异
       - 结合样本信息解释主成分的生物学意义
    
    4. 后续分析:
       - 可使用PCA结果作为协变量进行GWAS
       - 结合其他分析方法验证群体结构
       - 考虑进行更深入的群体遗传学分析
    
    引用和参考 | Citation & References:
    
    如果在学术研究中使用此工具，请引用相关的方法学文献:
    - PLINK: Purcell et al. (2007) AJHG
    - PCA in genetics: Patterson et al. (2006) PLoS Genet
    - Population structure: Pritchard et al. (2000) Genetics
    """
    
    # 构建参数列表传递给原始main函数
    args = ['vcf_pca.py']
    
    # 必需参数
    args.extend(['-v', vcf])
    
    # 可选参数（只在非默认值时添加）
    if output != './pca_output':
        args.extend(['-o', output])
    
    if sample_info:
        args.extend(['-s', sample_info])
    
    if components != 10:
        args.extend(['-c', str(components)])
    
    if maf != 0.05:
        args.extend(['-m', str(maf)])
    
    if missing != 0.1:
        args.extend(['--missing', str(missing)])
    
    if hwe != 1e-6:
        args.extend(['--hwe', str(hwe)])
    
    if group_column:
        args.extend(['-g', group_column])
    
    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])
    
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
    
    # 布尔选项
    if skip_qc:
        args.append('--skip-qc')
    
    if plot:
        args.append('-p')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        vcf_pca_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\nVCF PCA分析被用户中断 | VCF PCA analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"VCF PCA分析失败 | VCF PCA analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv