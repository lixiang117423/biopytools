"""
PLINK GWAS分析命令 | PLINK GWAS Analysis Command
"""

import click
import sys
from ...plink_gwas.main import main as plink_gwas_main


@click.command(short_help='PLINK全基因组关联分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf-file', '-v',
              required=True,
              type=click.Path(exists=True),
              help='VCF文件路径 | Path to VCF file')
@click.option('--phenotype-file', '-p',
              required=True,
              type=click.Path(exists=True),
              help='表型文件路径 | Path to phenotype file')
@click.option('--trait-type', '-t',
              default='qualitative',
              type=click.Choice(['qualitative', 'quantitative']),
              help='表型类型 | Trait type (default: qualitative)')
@click.option('--genetic-model', '-m',
              default='additive',
              type=click.Choice(['additive', 'dominant', 'recessive', 'all']),
              help='遗传模型 | Genetic model (default: additive)')
@click.option('--output-dir', '-o',
              default='gwas_results',
              type=click.Path(),
              help='输出目录 | Output directory (default: gwas_results)')
@click.option('--no-strat-corr',
              is_flag=True,
              help='禁用群体结构校正（跳过LD剪枝和PCA）| Disable population stratification correction (skip LD pruning and PCA)')
@click.option('--mind',
              default=0.05,
              type=float,
              help='个体缺失率阈值 | Individual missing rate threshold (default: 0.05)')
@click.option('--geno',
              default=0.05,
              type=float,
              help='SNP缺失率阈值 | SNP missing rate threshold (default: 0.05)')
@click.option('--maf',
              default=0.01,
              type=float,
              help='最小等位基因频率 | Minor allele frequency (default: 0.01)')
@click.option('--hwe',
              default=1e-6,
              type=float,
              help='Hardy-Weinberg平衡P值阈值 | HWE p-value threshold (default: 1e-6)')
@click.option('--ld-window-size',
              default=50,
              type=int,
              help='LD剪枝窗口大小(kb) | LD window size in kb (default: 50)')
@click.option('--ld-step-size',
              default=5,
              type=int,
              help='LD剪枝步长(SNP数) | LD step size in SNPs (default: 5)')
@click.option('--ld-r2-threshold',
              default=0.2,
              type=float,
              help='LD剪枝r²阈值 | LD r² threshold (default: 0.2)')
@click.option('--pca-components',
              default=10,
              type=int,
              help='计算的主成分数量 | Number of PCA components (default: 10)')
@click.option('--pca-use',
              default=5,
              type=int,
              help='关联分析中使用的主成分数量 | Number of PCs to use in association (default: 5)')
@click.option('--correction-method',
              default='all',
              type=click.Choice(['bonferroni', 'suggestive', 'fdr', 'all']),
              help='显著性校正方法 | Significance correction method (default: all)')
@click.option('--bonferroni-alpha',
              default=0.05,
              type=float,
              help='Bonferroni校正alpha水平 | Bonferroni alpha level (default: 0.05)')
@click.option('--suggestive-threshold',
              default=1e-5,
              type=float,
              help='提示性关联阈值 | Suggestive threshold (default: 1e-5)')
@click.option('--fdr-alpha',
              default=0.05,
              type=float,
              help='FDR校正q值阈值 | FDR q-value threshold (default: 0.05)')
@click.option('--threads',
              default=1,
              type=int,
              help='使用的线程数 | Number of threads (default: 1)')
def plinkgwas(vcf_file, phenotype_file, trait_type, genetic_model, output_dir,
               no_strat_corr, mind, geno, maf, hwe, ld_window_size, ld_step_size,
               ld_r2_threshold, pca_components, pca_use, correction_method,
               bonferroni_alpha, suggestive_threshold, fdr_alpha, threads):
    """
    PLINK GWAS分析工具
    
    使用PLINK进行全基因组关联分析，支持质量性状和数量性状，
    多种遗传模型，完整的数据质量控制和群体结构校正。
    
    功能特点 | Features:
    - 支持质量性状和数量性状分析
    - 多种遗传模型：加性、显性、隐性模型
    - 完整的数据质量控制流程
    - 群体结构分析和校正
    - 多种显著性校正方法
    - 丰富的结果可视化
    
    输出文件 | Output Files:
    - gwas_summary_report.txt: 总结报告
    - model_comparison_report.txt: 模型比较报告
    - gwas_results_*.txt: 详细关联分析结果
    - *.png: 各种可视化图表（曼哈顿图、QQ图等）
    
    示例 | Examples:
    
    \b
    # 基本质量性状分析（默认加性模型）
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt \\
        -t qualitative -o results
    
    \b
    # 使用显性模型分析
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt \\
        -t qualitative -m dominant -o results
    
    \b
    # 测试所有遗传模型
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt \\
        -t qualitative -m all -o results
    
    \b
    # 数量性状分析
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt \\
        -t quantitative -m additive -o results
    
    \b
    # 自定义质量控制参数
    biopytools plink-gwas -v genotypes.vcf.gz -p traits.txt \\
        --mind 0.02 --geno 0.02 --maf 0.05 --hwe 1e-8 \\
        -o strict_qc_results
    
    \b
    # 禁用群体结构校正的快速分析
    biopytools plink-gwas -v data.vcf.gz -p pheno.txt \\
        --no-strat-corr -o quick_analysis
    
    \b
    # 自定义LD剪枝和PCA参数
    biopytools plink-gwas -v large_dataset.vcf.gz -p phenotypes.txt \\
        --ld-window-size 25 --ld-r2-threshold 0.1 \\
        --pca-components 20 --pca-use 10 -o detailed_analysis
    
    \b
    # 多线程高性能分析
    biopytools plink-gwas -v genome.vcf.gz -p traits.txt \\
        --threads 32 --correction-method all \\
        -o high_performance_analysis
    
    表型文件格式 | Phenotype File Format:
    
    质量性状 (qualitative):
    FID    IID    pheno
    FAM1   IND1   1      # 对照
    FAM1   IND2   2      # 病例
    
    数量性状 (quantitative):
    FID    IID    pheno
    FAM1   IND1   1.25
    FAM1   IND2   2.67
    
    遗传模型说明 | Genetic Models:
    - additive: 加性模型，假设杂合子效应为纯合子效应的一半
    - dominant: 显性模型，假设至少一个变异等位基因产生效应
    - recessive: 隐性模型，假设需要两个变异等位基因才产生效应
    - all: 同时测试所有模型
    
    显著性校正方法 | Significance Correction:
    - bonferroni: Bonferroni多重检验校正
    - suggestive: 提示性关联阈值
    - fdr: 假发现率校正
    - all: 应用所有校正方法
    
    分析流程 | Analysis Pipeline:
    1. 输入文件格式转换和验证
    2. 表型数据处理和合并
    3. 数据质量控制（个体和SNP筛选）
    4. 群体结构分析（LD剪枝和PCA）
    5. 关联分析（多种遗传模型）
    6. 结果处理和显著性校正
    7. 报告生成和结果可视化
    
    质量控制参数 | Quality Control Parameters:
    - mind: 个体缺失率，超过此阈值的个体被移除
    - geno: SNP缺失率，超过此阈值的SNP被移除
    - maf: 最小等位基因频率，低于此频率的SNP被移除
    - hwe: Hardy-Weinberg平衡检验P值阈值
    
    性能建议 | Performance Tips:
    - 大数据集建议增加线程数提高计算速度
    - 严格的质量控制有助于提高结果可靠性
    - 群体结构校正对混合群体尤其重要
    - 选择合适的遗传模型基于先验生物学知识
    """
    
    # 构建参数列表传递给原始main函数
    args = ['plinkgwas.py']
    
    # 必需参数
    args.extend(['-v', vcf_file])
    args.extend(['-p', phenotype_file])
    
    # 分析参数（只在非默认值时添加）
    if trait_type != 'qualitative':
        args.extend(['-t', trait_type])
    
    if genetic_model != 'additive':
        args.extend(['-m', genetic_model])
    
    if output_dir != 'gwas_results':
        args.extend(['-o', output_dir])
    
    if no_strat_corr:
        args.append('--no-strat-corr')
    
    # 质量控制参数
    if mind != 0.05:
        args.extend(['--mind', str(mind)])
    
    if geno != 0.05:
        args.extend(['--geno', str(geno)])
    
    if maf != 0.01:
        args.extend(['--maf', str(maf)])
    
    if hwe != 1e-6:
        args.extend(['--hwe', str(hwe)])
    
    # LD剪枝参数
    if ld_window_size != 50:
        args.extend(['--ld-window-size', str(ld_window_size)])
    
    if ld_step_size != 5:
        args.extend(['--ld-step-size', str(ld_step_size)])
    
    if ld_r2_threshold != 0.2:
        args.extend(['--ld-r2-threshold', str(ld_r2_threshold)])
    
    # PCA参数
    if pca_components != 10:
        args.extend(['--pca-components', str(pca_components)])
    
    if pca_use != 5:
        args.extend(['--pca-use', str(pca_use)])
    
    # 显著性校正参数
    if correction_method != 'all':
        args.extend(['--correction-method', correction_method])
    
    if bonferroni_alpha != 0.05:
        args.extend(['--bonferroni-alpha', str(bonferroni_alpha)])
    
    if suggestive_threshold != 1e-5:
        args.extend(['--suggestive-threshold', str(suggestive_threshold)])
    
    if fdr_alpha != 0.05:
        args.extend(['--fdr-alpha', str(fdr_alpha)])
    
    # 计算参数
    if threads != 1:
        args.extend(['--threads', str(threads)])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        plink_gwas_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断分析 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"分析失败 | Analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv