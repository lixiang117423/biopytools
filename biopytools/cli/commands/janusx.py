"""
JanusX GWAS和基因组选择分析命令|JanusX GWAS and Genomic Selection Analysis Command
"""

import click
import sys
import os


def _lazy_import_janusx_main():
    """延迟加载janusx主函数|Lazy load janusx main function"""
    try:
        from ...janusx.main import main as janusx_main
        return janusx_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.group(
    short_help='JanusX GWAS和基因组选择分析|JanusX GWAS and Genomic Selection Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def janusx():
    """
    JanusX GWAS和基因组选择分析工具|JanusX GWAS and Genomic Selection Analysis Tool

    全基因组关联分析和基因组选择的统一接口|Unified interface for GWAS and Genomic Selection

    示例|Example: biopytools janusx gwas -i data.vcf.gz -p pheno.txt -o output
    """
    pass


@janusx.command(short_help='全基因组关联分析|Genome-Wide Association Study')
@click.option('-i', '--genotype',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
@click.option('-p', '--pheno',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型文件|Phenotype file')
@click.option('-t', '--type',
              default='vcf',
              type=click.Choice(['vcf', 'bfile']),
              help='基因型类型|Genotype type')
@click.option('-m', '--models',
              multiple=True,
              default=['lmm'],
              show_default=False,
              type=click.Choice(['lm', 'lmm', 'fastlmm', 'farmcpu']),
              help='GWAS模型|GWAS models')
@click.option('--maf',
              type=float,
              default=0.05,
              show_default=True,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--geno',
              type=float,
              default=0.0,
              show_default=True,
              help='缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)')
@click.option('-k', '--grm',
              default='1',
              show_default=True,
              help='亲缘关系矩阵方法|GRM method')
@click.option('-q', '--qcov',
              type=int,
              default=0,
              show_default=True,
              help='PCA数量或Q矩阵文件路径|Number of PCs or path to Q matrix file (0=不使用Q矩阵|no Q matrix)')
@click.option('-c', '--cov',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='协变量文件|Covariate file')
@click.option('-n', '--ncol',
              multiple=True,
              type=int,
              help='表型列索引(零基)|Phenotype column indices (zero-based)')
@click.option('--plot',
              is_flag=True,
              help='生成图表|Generate plots')
@click.option('--chunksize',
              type=int,
              default=100000,
              show_default=True,
              help='SNP分块大小|SNP chunk size')
@click.option('--mmap-limit',
              type=int,
              help='内存映射限制|Memory map limit')
@click.option('-th', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-o', '--output-dir',
              default='./janusx_gwas_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix',
              help='输出文件前缀|Output file prefix')
@click.option('--janusx-path',
              type=click.Path(),
              help='JanusX可执行文件路径|JanusX executable path')
def gwas(genotype, pheno, type, models, maf, geno, grm, qcov, cov, ncol,
         plot, chunksize, mmap_limit, threads, output_dir, prefix, janusx_path):
    """
    全基因组关联分析|Genome-Wide Association Study (GWAS)

    使用GLM、LMM、fastLMM或FarmCPU进行GWAS分析|Perform GWAS using GLM, LMM, fastLMM, or FarmCPU

    示例|Examples: biopytools janusx gwas -i data.vcf.gz -p pheno.txt -m lmm -o output

      # 基本LMM分析|Basic LMM analysis:

      # 多模型比较|Multi-model comparison:
      biopytools janusx gwas -i data.vcf.gz -p pheno.txt -m lmm -m farmcpu -o output
    """
    janusx_main = _lazy_import_janusx_main()

    args = _build_gwas_args({
        'genotype': genotype,
        'pheno': pheno,
        'type': type,
        'models': list(models),
        'maf': maf,
        'geno': geno,
        'grm': grm,
        'qcov': qcov,
        'cov': cov,
        'ncol': list(ncol) if ncol else None,
        'plot': plot,
        'chunksize': chunksize,
        'mmap_limit': mmap_limit,
        'threads': threads,
        'output_dir': output_dir,
        'prefix': prefix,
        'janusx_path': janusx_path
    })

    _execute_main(janusx_main, args)


@janusx.command(short_help='基因组选择分析|Genomic Selection')
@click.option('-i', '--genotype',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
@click.option('-p', '--pheno',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='表型文件|Phenotype file')
@click.option('-t', '--type',
              default='vcf',
              type=click.Choice(['vcf', 'bfile']),
              help='基因型类型|Genotype type')
@click.option('-m', '--models',
              multiple=True,
              default=['GBLUP'],
              show_default=False,
              type=click.Choice(['GBLUP', 'rrBLUP', 'BayesA', 'BayesB', 'BayesCpi']),
              help='GS模型|GS models')
@click.option('--maf',
              type=float,
              default=0.05,
              show_default=True,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--geno',
              type=float,
              default=0.0,
              show_default=True,
              help='缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)')
@click.option('--pcd',
              is_flag=True,
              help='启用PCA降维|Enable PCA-based dimensionality reduction')
@click.option('-n', '--ncol',
              type=int,
              help='表型列索引(零基)|Phenotype column index (zero-based)')
@click.option('--cv',
              type=int,
              help='交叉验证折数|K-fold cross-validation')
@click.option('--plot',
              is_flag=True,
              help='生成图表|Generate plots')
@click.option('-o', '--output-dir',
              default='./janusx_gs_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix',
              help='输出文件前缀|Output file prefix')
@click.option('--janusx-path',
              type=click.Path(),
              help='JanusX可执行文件路径|JanusX executable path')
def gs(genotype, pheno, type, models, maf, geno, pcd, ncol, cv, plot, output_dir, prefix, janusx_path):
    """
    基因组选择分析|Genomic Selection (GS)

    使用GBLUP、rrBLUP或贝叶斯方法进行基因组预测|Perform genomic prediction using GBLUP, rrBLUP, or Bayesian methods

    示例|Examples: biopytools janusx gs -i data.vcf.gz -p pheno.txt -m GBLUP -o output

      # 基本GBLUP分析|Basic GBLUP analysis:

      # 多模型比较|Multi-model comparison:
      biopytools janusx gs -i data.vcf.gz -p pheno.txt -m GBLUP -m rrBLUP -o output
    """
    janusx_main = _lazy_import_janusx_main()

    args = _build_gs_args({
        'genotype': genotype,
        'pheno': pheno,
        'type': type,
        'models': list(models),
        'maf': maf,
        'geno': geno,
        'pcd': pcd,
        'ncol': ncol,
        'cv': cv,
        'plot': plot,
        'output_dir': output_dir,
        'prefix': prefix,
        'janusx_path': janusx_path
    })

    _execute_main(janusx_main, args)


@janusx.command(short_help='主成分分析|Principal Component Analysis')
@click.option('-i', '--genotype',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
@click.option('-t', '--type',
              default='vcf',
              type=click.Choice(['vcf', 'bfile']),
              help='基因型类型|Genotype type')
@click.option('--grm',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='预计算的GRM前缀|Precomputed GRM prefix')
@click.option('--pcfile',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='预计算的PCA文件|Precomputed PCA file')
@click.option('-d', '--dim',
              type=int,
              default=3,
              show_default=True,
              help='输出主成分数量|Number of PCs to output')
@click.option('--plot',
              is_flag=True,
              help='生成2D散点图|Generate 2D scatter plot')
@click.option('--plot3d',
              is_flag=True,
              help='生成3D旋转GIF|Generate 3D rotating GIF')
@click.option('-g', '--group',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='分组文件路径|Group file path')
@click.option('--color',
              type=int,
              default=1,
              show_default=True,
              help='调色板索引(0-6)|Color palette index')
@click.option('-o', '--output-dir',
              default='./janusx_pca_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix',
              help='输出文件前缀|Output file prefix')
@click.option('--janusx-path',
              type=click.Path(),
              help='JanusX可执行文件路径|JanusX executable path')
def pca(genotype, type, grm, pcfile, dim, plot, plot3d, group, color, output_dir, prefix, janusx_path):
    """
    主成分分析|Principal Component Analysis (PCA)

    进行主成分分析并可视化|Perform PCA and generate visualizations

    示例|Examples: biopytools janusx pca -i data.vcf.gz -d 5 -o output

      # 基本PCA分析|Basic PCA analysis:

      # 带可视化的PCA|PCA with plots:
      biopytools janusx pca -i data.vcf.gz -d 5 --plot --plot3d -o output
    """
    janusx_main = _lazy_import_janusx_main()

    args = _build_pca_args({
        'genotype': genotype,
        'type': type,
        'grm': grm,
        'pcfile': pcfile,
        'dim': dim,
        'plot': plot,
        'plot3d': plot3d,
        'group': group,
        'color': color,
        'output_dir': output_dir,
        'prefix': prefix,
        'janusx_path': janusx_path
    })

    _execute_main(janusx_main, args)


@janusx.command(short_help='GWAS结果可视化和注释|GWAS result visualization and annotation')
@click.option('-f', '--files',
              multiple=True,
              required=True,
              callback=lambda ctx, param, value: [_validate_file_exists(f) for f in value] if value else None,
              help='GWAS结果文件列表|List of GWAS result files')
@click.option('--chr-col',
              default='#CHROM',
              show_default=True,
              help='染色体列名|Chromosome column name')
@click.option('--pos-col',
              default='POS',
              show_default=True,
              help='位置列名|Position column name')
@click.option('--pvalue-col',
              default='p',
              show_default=True,
              help='P值列名|P-value column name')
@click.option('--threshold',
              type=float,
              help='显著性阈值|Significance threshold')
@click.option('--noplot',
              is_flag=True,
              help='禁用绘图|Disable plotting')
@click.option('--color',
              type=int,
              default=0,
              show_default=True,
              help='颜色方案索引(0-6, -1为auto)|Color scheme index')
@click.option('--hl', '--highlight',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='高亮区域BED文件|Highlight regions BED file')
@click.option('--format',
              default='png',
              show_default=True,
              type=click.Choice(['pdf', 'png', 'svg', 'tif']),
              help='输出格式|Output format')
@click.option('-a', '--anno',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='注释文件路径|Annotation file path')
@click.option('--ab', '--anno-broaden',
              type=int,
              help='注释窗口|Annotation window (kb)')
@click.option('--desc-item',
              default='description',
              show_default=True,
              help='GFF描述键|GFF description key')
@click.option('-o', '--output-dir',
              default='./janusx_postgwas_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix',
              default='JanusX',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('-th', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--janusx-path',
              type=click.Path(),
              help='JanusX可执行文件路径|JanusX executable path')
def postgwas(files, chr_col, pos_col, pvalue_col, threshold, noplot, color,
              hl, format, anno, ab, desc_item, output_dir, prefix, threads, janusx_path):
    """
    GWAS结果可视化和注释|GWAS result visualization and annotation

    生成Manhattan图、QQ图并进行SNP注释|Generate Manhattan plots, QQ plots and perform SNP annotation

    示例|Examples: biopytools janusx postgwas -f result.lmm.tsv --threshold 1e-6 -o output

      # 基本可视化|Basic visualization:

      # 带注释的可视化|Visualization with annotation:
      biopytools janusx postgwas -f result.lmm.tsv -a genes.gff --threshold 1e-6 -o output
    """
    janusx_main = _lazy_import_janusx_main()

    args = _build_postgwas_args({
        'files': list(files),
        'chr_col': chr_col,
        'pos_col': pos_col,
        'pvalue_col': pvalue_col,
        'threshold': threshold,
        'noplot': noplot,
        'color': color,
        'highlight': hl,  # 将hl参数映射到highlight键
        'format': format,
        'anno': anno,
        'anno_broaden': ab,  # 将ab参数映射到anno_broaden键
        'desc_item': desc_item,
        'output_dir': output_dir,
        'prefix': prefix,
        'threads': threads,
        'janusx_path': janusx_path
    })

    _execute_main(janusx_main, args)


def _build_gwas_args(params):
    """构建GWAS参数列表|Build GWAS argument list"""
    args = ['janusx.py', 'gwas']

    # 必需参数|Required parameters
    args.extend(['-i', params['genotype']])
    args.extend(['-p', params['pheno']])

    # 可选参数|Optional parameters
    if params['type'] != 'vcf':
        args.extend(['-t', params['type']])

    if params['models']:
        for model in params['models']:
            args.extend(['-m', model])

    if params['maf'] != 0.05:
        args.extend(['--maf', str(params['maf'])])

    if params['geno'] != 0.0:
        args.extend(['--geno', str(params['geno'])])

    if params['grm'] != '1':
        args.extend(['-k', params['grm']])

    if params['qcov'] != 0:
        args.extend(['-q', str(params['qcov'])])

    if params['cov']:
        args.extend(['-c', params['cov']])

    if params['ncol']:
        for col in params['ncol']:
            args.extend(['-n', str(col)])

    if params['plot']:
        args.append('--plot')

    if params['chunksize'] != 100000:
        args.extend(['--chunksize', str(params['chunksize'])])

    if params['mmap_limit'] is not None:
        args.extend(['--mmap-limit', str(params['mmap_limit'])])

    if params['threads'] != 12:
        args.extend(['-th', str(params['threads'])])

    if params['output_dir'] != './janusx_gwas_output':
        args.extend(['-o', params['output_dir']])

    if params['prefix']:
        args.extend(['--prefix', params['prefix']])

    if params['janusx_path']:
        args.extend(['--janusx-path', params['janusx_path']])

    return args


def _build_gs_args(params):
    """构建GS参数列表|Build GS argument list"""
    args = ['janusx.py', 'gs']

    # 必需参数|Required parameters
    args.extend(['-i', params['genotype']])
    args.extend(['-p', params['pheno']])

    # 可选参数|Optional parameters
    if params['type'] != 'vcf':
        args.extend(['-t', params['type']])

    if params['models']:
        for model in params['models']:
            args.extend(['-m', model])

    if params['maf'] != 0.05:
        args.extend(['--maf', str(params['maf'])])

    if params['geno'] != 0.0:
        args.extend(['--geno', str(params['geno'])])

    if params['pcd']:
        args.append('--pcd')

    if params['ncol'] is not None:
        args.extend(['-n', str(params['ncol'])])

    if params['cv'] is not None:
        args.extend(['--cv', str(params['cv'])])

    if params['plot']:
        args.append('--plot')

    if params['output_dir'] != './janusx_gs_output':
        args.extend(['-o', params['output_dir']])

    if params['prefix']:
        args.extend(['--prefix', params['prefix']])

    if params['janusx_path']:
        args.extend(['--janusx-path', params['janusx_path']])

    return args


def _build_pca_args(params):
    """构建PCA参数列表|Build PCA argument list"""
    args = ['janusx.py', 'pca']

    # 输入源|Input source
    if params['grm']:
        args.extend(['--grm', params['grm']])
    elif params['pcfile']:
        args.extend(['--pcfile', params['pcfile']])
    elif params['genotype']:
        args.extend(['-i', params['genotype']])
        if params['type'] != 'vcf':
            args.extend(['-t', params['type']])

    # 可选参数|Optional parameters
    if params['dim'] != 3:
        args.extend(['-d', str(params['dim'])])

    if params['plot']:
        args.append('--plot')

    if params['plot3d']:
        args.append('--plot3d')

    if params['group']:
        args.extend(['-g', params['group']])

    if params['color'] != 1:
        args.extend(['--color', str(params['color'])])

    if params['output_dir'] != './janusx_pca_output':
        args.extend(['-o', params['output_dir']])

    if params['prefix']:
        args.extend(['--prefix', params['prefix']])

    if params['janusx_path']:
        args.extend(['--janusx-path', params['janusx_path']])

    return args


def _build_postgwas_args(params):
    """构建PostGWAS参数列表|Build PostGWAS argument list"""
    args = ['janusx.py', 'postgwas']

    # 必需参数|Required parameters
    for file in params['files']:
        args.extend(['-f', file])

    # 可选参数|Optional parameters
    if params['chr_col'] != '#CHROM':
        args.extend(['--chr-col', params['chr_col']])

    if params['pos_col'] != 'POS':
        args.extend(['--pos-col', params['pos_col']])

    if params['pvalue_col'] != 'p':
        args.extend(['--pvalue-col', params['pvalue_col']])

    if params['threshold'] is not None:
        args.extend(['--threshold', str(params['threshold'])])

    if params['noplot']:
        args.append('--noplot')

    if params['color'] != 0:
        args.extend(['--color', str(params['color'])])

    if params['highlight']:
        args.extend(['--hl', params['highlight']])

    if params['format'] != 'png':
        args.extend(['--format', params['format']])

    if params['anno']:
        args.extend(['-a', params['anno']])

    if params['anno_broaden'] is not None:
        args.extend(['--ab', str(params['anno_broaden'])])

    if params['desc_item'] != 'description':
        args.extend(['--desc-item', params['desc_item']])

    if params['output_dir'] != './janusx_postgwas_output':
        args.extend(['-o', params['output_dir']])

    if params['prefix'] != 'JanusX':
        args.extend(['--prefix', params['prefix']])

    if params['threads'] != 12:
        args.extend(['-th', str(params['threads'])])

    if params['janusx_path']:
        args.extend(['--janusx-path', params['janusx_path']])

    return args


def _execute_main(janusx_main, args):
    """执行主程序|Execute main program"""
    original_argv = sys.argv
    sys.argv = args

    try:
        janusx_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
