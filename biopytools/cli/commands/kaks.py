"""
Ka/Ks计算分析CLI包装器|Ka/Ks Calculation Analysis CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_kaks_main():
    """延迟加载Ka/Ks主函数|Lazy load Ka/Ks main function"""
    try:
        from ...kaks.main import main as kaks_main
        return kaks_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
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


@click.command(
    short_help='Ka/Ks计算分析|Ka/Ks Calculation Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--fasta1', '-1',
              required=True,
              type=click.Path(exists=True),
              help='第一个FASTA文件|First FASTA file')
@click.option('--fasta2', '-2',
              required=True,
              type=click.Path(exists=True),
              help='第二个FASTA文件|Second FASTA file')
@click.option('--pairs', '-p',
              required=True,
              type=click.Path(exists=True),
              help='序列配对文件|Sequence pair file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--method', '-m',
              default='GMYN',
              show_default=True,
              type=click.Choice([
                  'GMYN', 'MYN', 'YN', 'NG', 'LWL', 'LPB', 'MLWL', 'MLPB',
                  'GY', 'MS', 'MA', 'GNG', 'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN'
              ]),
              help='计算方法|Calculation method')
@click.option('--kaks-path',
              default='KaKs_Calculator',
              show_default=True,
              type=str,
              help='KaKs_Calculator软件路径|Path to KaKs_Calculator executable')
@click.option('--verbose', '-v',
              is_flag=True,
              help='启用详细日志|Enable verbose logging')
@click.option('--temp-dir',
              type=click.Path(),
              help='自定义临时目录|Custom temporary directory')
@click.option('--keep-temp',
              is_flag=True,
              help='保留临时文件|Keep temporary files')
def kaks(fasta1, fasta2, pairs, output, method, kaks_path, verbose, temp_dir, keep_temp):
    """
    Ka/Ks计算分析流程|Ka/Ks Calculation Analysis Pipeline

    计算同源基因对之间的非同义替换率(Ka)和同义替换率(Ks)
    Calculate non-synonymous (Ka) and synonymous (Ks) substitution rates between homologous gene pairs

    示例|Examples: biopytools kaks -1 species1.fasta -2 species2.fasta -p pairs.txt -o results/
    """

    # 延迟加载|Lazy load
    kaks_main = _lazy_import_kaks_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['kaks.py']

    # 必需参数|Required parameters
    args.extend(['-1', fasta1])
    args.extend(['-2', fasta2])
    args.extend(['-p', pairs])
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if method != 'GMYN':
        args.extend(['-m', method])

    if kaks_path != 'KaKs_Calculator':
        args.extend(['--kaks-path', kaks_path])

    if verbose:
        args.append('--verbose')

    if temp_dir:
        args.extend(['--temp-dir', temp_dir])

    if keep_temp:
        args.append('--keep-temp')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        kaks_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
