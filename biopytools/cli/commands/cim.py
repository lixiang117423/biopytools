"""
CIM分析命令|R/qtl CIM Analysis Command
"""

import click
import sys
import os


def _lazy_import_cim_main():
    """延迟加载cim主函数|Lazy load cim main function"""
    try:
        from ...cim.main import main as cim_main
        return cim_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="R/qtl复合区间作图(CIM)|R/qtl Composite Interval Mapping",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--input', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='输入VCF文件|Input VCF file (.vcf/.vcf.gz)')
@click.option('-p', '--pheno', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='表型文件|Phenotype file (TSV: sample, value)')
@click.option('-o', '--output', required=True,
              help='输出目录|Output directory')
@click.option('-t', '--type', 'cross_type',
              type=click.Choice(['f2', 'bc']), default='f2', show_default=True,
              help='群体类型|Cross type (f2/bc)')
@click.option('--map-mode',
              type=click.Choice(['physical', 'estimate', 'mstmap']), default='mstmap', show_default=True,
              help='cM位置来源|cM source (physical/estimate/mstmap)')
@click.option('--maf', type=float, default=0.05, show_default=True,
              help='最小等位基因频率|Minor allele frequency threshold')
@click.option('--missing', type=float, default=0.1, show_default=True,
              help='最大缺失率|Maximum missing rate')
@click.option('--n-marcovar', type=int, default=10, show_default=True,
              help='协因子数量|Number of marker covariates')
@click.option('--window', type=float, default=10.0, show_default=True,
              help='窗口大小(cM)|Window size in cM')
@click.option('--method',
              type=click.Choice(['hk', 'em', 'imp']), default='hk', show_default=True,
              help='扫描方法|Scanning method')
@click.option('--step', type=float, default=1.0, show_default=True,
              help='伪标记步长(cM)|Pseudomarker step in cM')
@click.option('--n-perm', type=int, default=1000, show_default=True,
              help='置换检验次数(0=跳过)|Permutation count (0=skip)')
@click.option('--ld-window', type=int, default=50, show_default=True,
              help='LD窗口(SNP数)|LD window in SNP count')
@click.option('--ld-step', type=int, default=5, show_default=True,
              help='LD步长(SNP数)|LD step in SNP count')
@click.option('--ld-r2', type=float, default=0.1, show_default=True,
              help='LD r2阈值|LD r2 threshold')
@click.option('--skip-ld', is_flag=True, default=True, show_default=True,
              help='跳过LD降维|Skip LD pruning (default: True)')
@click.option('--no-skip-ld', 'skip_ld', is_flag=True, default=None, flag_value=False,
              help='启用LD降维|Enable LD pruning')
@click.option('--mstmap-pvalue', type=float, default=1e-6, show_default=True,
              help='MSTmap聚类p值起始值(自动调优)|MSTmap clustering p-value start, auto-tuned')
@click.option('--mstmap-distfun',
              type=click.Choice(['kosambi', 'haldane']), default='kosambi', show_default=True,
              help='MSTmap距离函数|MSTmap distance function')
@click.option('--r-env', default='Rqtl', show_default=True,
              help='R conda环境名|R conda environment name')
def cim(input, pheno, output, cross_type, map_mode, maf, missing,
         n_marcovar, window, method, step, n_perm,
         ld_window, ld_step, ld_r2, skip_ld, mstmap_pvalue, mstmap_distfun, r_env):
    """R/qtl复合区间作图(CIM)分析|R/qtl Composite Interval Mapping (CIM)

    从VCF和表型文件自动完成CIM分析：标记过滤、LD降维、CIM扫描、置换检验、结果可视化
    Automated CIM from VCF + phenotype: marker filtering, LD pruning, CIM scan, permutation test, visualization

    示例|Examples: biopytools cim -i input.vcf.gz -p phe.txt -o output_dir
    """

    cim_main_func = _lazy_import_cim_main()

    args = ['cim.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-p', pheno])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    args.extend(['-t', cross_type])
    args.extend(['--map-mode', map_mode])
    args.extend(['--maf', str(maf)])
    args.extend(['--missing', str(missing)])
    args.extend(['--n-marcovar', str(n_marcovar)])
    args.extend(['--window', str(window)])
    args.extend(['--method', method])
    args.extend(['--step', str(step)])
    args.extend(['--n-perm', str(n_perm)])
    args.extend(['--ld-window', str(ld_window)])
    args.extend(['--ld-step', str(ld_step)])
    args.extend(['--ld-r2', str(ld_r2)])
    if skip_ld:
        args.append('--skip-ld')
    args.extend(['--mstmap-pvalue', str(mstmap_pvalue)])
    args.extend(['--mstmap-distfun', mstmap_distfun])
    args.extend(['--r-env', r_env])

    original_argv = sys.argv
    sys.argv = args

    try:
        cim_main_func()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
