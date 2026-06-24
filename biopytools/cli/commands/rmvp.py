"""
rMVP GWAS命令|rMVP GWAS Analysis Command
"""

import click
import sys
import os


def _lazy_import_rmvp_main():
    """延迟加载rmvp主函数|Lazy load rmvp main function"""
    try:
        from ...rmvp.cli import main as rmvp_main
        return rmvp_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


def _validate_directory_exists(dir_path):
    """验证目录存在(如不存在则创建)|Validate directory exists (create if not exists)"""
    if not _is_help_request() and dir_path:
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"无法创建输出目录|Cannot create output directory: {dir_path}")
    return dir_path


@click.command(short_help="rMVP GWAS批量分析工具|rMVP Batch GWAS Analysis Tool",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件|Input VCF file')
@click.option('--pheno', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入表型文件|Input phenotype file')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='输出目录|Output directory')
@click.option('--output-prefix',
              default='RMVP_Result',
              show_default=True,
              help='输出前缀|Output prefix')
@click.option('--models',
              multiple=True,
              type=click.Choice(['GLM', 'MLM', 'FarmCPU']),
              default=['GLM', 'MLM', 'FarmCPU'],
              show_default=True,
              help='分析模型|Analysis models (可多选|can select multiple)')
@click.option('--r-env',
              default='~/miniforge3/envs/rMVP',
              show_default=True,
              help='R conda环境路径|R conda environment path (e.g., ~/miniforge3/envs/rMVP)')
@click.option('--r-path',
              help='R可执行文件路径|R executable path')
@click.option('--ncpus',
              type=int,
              default=12,
              show_default=True,
              help='CPU核心数|Number of CPU cores')
@click.option('--maxLine',
              type=int,
              default=10000,
              show_default=True,
              help='每次读取的SNP数量|Number of SNPs to read at once (较小值减少内存|smaller uses less memory)')
@click.option('--n-pc-glm',
              type=int,
              default=3,
              show_default=True,
              help='GLM模型使用的PC数量|Number of PCs for GLM')
@click.option('--n-pc-mlm',
              type=int,
              default=3,
              show_default=True,
              help='MLM模型使用的PC数量|Number of PCs for MLM')
@click.option('--n-pc-farmcpu',
              type=int,
              default=3,
              show_default=True,
              help='FarmCPU模型使用的PC数量|Number of PCs for FarmCPU')
@click.option('--vc-method',
              type=click.Choice(['BRENT', 'EMMA', 'HE']),
              default='BRENT',
              show_default=True,
              help='MLM方差组分分析方法|MLM variance component method')
@click.option('--max-loop',
              type=int,
              default=10,
              show_default=True,
              help='FarmCPU最大迭代次数|FarmCPU max iterations')
@click.option('--method-bin',
              type=click.Choice(['static', 'fast-lmm']),
              default='static',
              show_default=True,
              help='FarmCPU bin方法|FarmCPU bin method')
@click.option('--maf',
              type=float,
              help='最小等位基因频率阈值|Minor allele frequency threshold')
@click.option('--miss',
              type=float,
              help='缺失率阈值|Missing rate threshold')
@click.option('--file-type',
              type=click.Choice(['jpg', 'pdf', 'tiff']),
              default='jpg',
              show_default=True,
              help='图片格式|Figure format')
@click.option('--dpi',
              type=int,
              default=300,
              show_default=True,
              help='图片分辨率|Figure DPI')
@click.option('--threshold',
              type=float,
              default=0.05,
              show_default=True,
              help='显著性阈值|Significance threshold')
@click.option('--ld-pruning/--no-ld-pruning',
              'ld_pruning',
              default=True,
              show_default=True,
              help='LD去连锁（默认开启）：kinship/PCA在去连锁SNP上计算，GWAS用全部SNP|LD pruning (default on): K/PCA on pruned SNPs, GWAS uses all SNPs')
@click.option('--ld-window',
              default='3000kb',
              show_default=True,
              help='LD修剪窗口|LD pruning window (e.g. 3000kb or 500)')
@click.option('--ld-step',
              type=int,
              default=1,
              show_default=True,
              help='LD修剪步长|LD pruning step size')
@click.option('--ld-r2',
              type=float,
              default=0.2,
              show_default=True,
              help='LD r2阈值|LD r2 threshold')
@click.option('--plink-path',
              default=None,
              help='PLINK可执行文件路径|PLINK executable path (default: conda env Population_genetics)')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARN', 'ERROR']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
def rmvp(**kwargs):
    """
    rMVP GWAS批量分析工具|rMVP Batch GWAS Analysis Tool

    使用rMVP进行全基因组关联分析(GWAS)，支持GLM、MLM、FarmCPU三个模型
    Perform GWAS using rMVP with GLM, MLM, and FarmCPU models

    示例|Examples: biopytools rmvp -i input.vcf.gz -p phenotype.txt -o output --r-env ~/miniforge3/envs/rMVP
    """
    # 获取主函数|Get main function
    rmvp_main = _lazy_import_rmvp_main()

    # 构建参数列表|Build arguments list
    argv = []

    # 必需参数|Required arguments
    argv.extend(['--vcf', kwargs['vcf']])
    argv.extend(['--pheno', kwargs['pheno']])
    argv.extend(['--output', kwargs['output']])

    # 可选参数|Optional arguments
    if kwargs['output_prefix'] != 'RMVP_Result':
        argv.extend(['--output-prefix', kwargs['output_prefix']])

    if kwargs['models']:
        argv.extend(['--models'] + list(kwargs['models']))

    if kwargs['r_env']:
        argv.extend(['--r-env', kwargs['r_env']])

    if kwargs['r_path']:
        argv.extend(['--r-path', kwargs['r_path']])

    argv.extend(['--ncpus', str(kwargs['ncpus'])])
    argv.extend(['--maxLine', str(kwargs['maxline'])])
    argv.extend(['--n-pc-glm', str(kwargs['n_pc_glm'])])
    argv.extend(['--n-pc-mlm', str(kwargs['n_pc_mlm'])])
    argv.extend(['--n-pc-farmcpu', str(kwargs['n_pc_farmcpu'])])
    argv.extend(['--vc-method', kwargs['vc_method']])
    argv.extend(['--max-loop', str(kwargs['max_loop'])])
    argv.extend(['--method-bin', kwargs['method_bin']])

    if kwargs['maf'] is not None:
        argv.extend(['--maf', str(kwargs['maf'])])

    if kwargs['miss'] is not None:
        argv.extend(['--miss', str(kwargs['miss'])])

    argv.extend(['--file-type', kwargs['file_type']])
    argv.extend(['--dpi', str(kwargs['dpi'])])
    argv.extend(['--threshold', str(kwargs['threshold'])])

    # LD去连锁参数（仅非默认值透传）|LD pruning params (non-defaults only)
    if not kwargs['ld_pruning']:
        argv.append('--no-ld-pruning')
    if kwargs['ld_window'] != '3000kb':
        argv.extend(['--ld-window', kwargs['ld_window']])
    if kwargs['ld_step'] != 1:
        argv.extend(['--ld-step', str(kwargs['ld_step'])])
    if kwargs['ld_r2'] != 0.2:
        argv.extend(['--ld-r2', str(kwargs['ld_r2'])])
    if kwargs['plink_path']:
        argv.extend(['--plink-path', kwargs['plink_path']])

    argv.extend(['--log-level', kwargs['log_level']])

    if kwargs['quiet']:
        argv.append('--quiet')

    # 修改sys.argv并调用主函数|Modify sys.argv and call main function
    sys.argv = ['rmvp'] + argv
    rmvp_main()
