"""
PopLD LD衰减分析命令 | PopLD LD Decay Analysis Command
"""

import click
import sys
import os


def _lazy_import_popld_main():
    """懒加载popld main函数 | Lazy load popld main function"""
    try:
        from ...popLD.main import main as popld_main
        return popld_main
    except ImportError as e:
        click.echo(f"导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="PopLD连锁不平衡衰减分析",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径 | Input VCF file path')
@click.option('--output-stat', '-o',
              required=True,
              help='输出统计文件前缀 | Output statistics file prefix')
@click.option('--poplddecay-path', '-p',
              default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay',
              help='PopLD软件路径 | PopLD software path')
@click.option('--sub-pop', '-s',
              help='亚群样本列表文件 | Subgroup sample list file')
@click.option('--max-dist', '-d',
              type=int,
              default=300,
              help='SNP间最大距离(kb) | Max distance between SNPs in kb (default: 300)')
@click.option('--maf', '-m',
              type=float,
              default=0.005,
              help='最小次等位基因频率 | Min minor allele frequency (default: 0.005)')
@click.option('--het',
              type=float,
              default=0.88,
              help='最大杂合率 | Max ratio of het allele (default: 0.88)')
@click.option('--miss',
              type=float,
              default=0.25,
              help='最大缺失率 | Max ratio of miss allele (default: 0.25)')
@click.option('--ehh',
              help='EHH起始位点 | EHH start site (NA for disabled)')
@click.option('--out-type',
              type=click.Choice(['1', '2', '3', '4', '5', '6', '7', '8']),
              default='1',
              help='输出类型 | Output type (default: 1)')
@click.option('--method',
              type=click.Choice(['1', '2']),
              default='1',
              help='算法方法(1=快速, 2=可能需要更多内存) | Algorithm method (default: 1)')
@click.option('--out-filter-snp',
              is_flag=True,
              help='输出最终用于计算的SNP | Output final SNP used for calculation')
def popld(input_vcf, output_stat, poplddecay_path, sub_pop, max_dist,
               maf, het, miss, ehh, out_type, method, out_filter_snp):
    """
    PopLD连锁不平衡衰减分析工具.

    对VCF格式的SNP数据进行连锁不平衡(LD)衰减分析，计算不同距离下
    SNP位点间的LD水平并生成统计结果。

    示例 | Examples:

    \b
    # 基本用法
    biopytools poplddecay -i snp.vcf.gz -o LD_result

    \b
    # 使用亚群样本
    biopytools poplddecay -i snp.vcf.gz -o LD_result -s subgroup.txt

    \b
    # 自定义参数
    biopytools poplddecay \\
        -i snp.vcf.gz \\
        -o LD_result \\
        -d 500 \\
        -m 0.01 \\
        --het 0.9

    \b
    # 指定软件路径
    biopytools poplddecay \\
        -i snp.vcf.gz \\
        -o LD_result \\
        -p /path/to/PopLDdecay
    """

    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    popld_main = _lazy_import_popld_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['popld.py']

    # 必需参数 | Required parameters
    args.extend(['-i', input_vcf])
    args.extend(['-o', output_stat])

    # 可选参数（只在非默认值时添加）| Optional parameters (add only when non-default)
    if poplddecay_path != '/share/org/YZWL/yzwl_lixg/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay':
        args.extend(['-p', poplddecay_path])

    if sub_pop:
        args.extend(['-s', sub_pop])

    if max_dist != 300:
        args.extend(['-d', str(max_dist)])

    if maf != 0.005:
        args.extend(['-m', str(maf)])

    if het != 0.88:
        args.extend(['--het', str(het)])

    if miss != 0.25:
        args.extend(['--miss', str(miss)])

    if ehh:
        args.extend(['--ehh', ehh])

    if out_type != '1':
        args.extend(['--out-type', out_type])

    if method != '1':
        args.extend(['--method', method])

    # 处理选项（布尔标志）| Processing options (boolean flags)
    if out_filter_snp:
        args.append('--out-filter-snp')

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        popld_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
