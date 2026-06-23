"""
亚基因组归属命令|Subgenome Assignment Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...subgenome_assign.main import main as subgenome_main
        return subgenome_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path existence"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


def _validate_parent(ctx, param, value):
    """验证 --parent 参数格式|Validate --parent format"""
    if _is_help_request():
        return value
    if not value:
        return value
    for spec in value:
        if ':' not in spec:
            raise click.BadParameter(
                f"--parent 格式错误|invalid format: {spec}. "
                f"应为|expected NAME:hap1.fa,hap2.fa"
            )
    return value


@click.command(
    short_help='亚基因组归属（基于亲本比对）|Subgenome assignment via parental alignment',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--target',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='目标多倍体基因组 FASTA|Target polyploid genome FASTA')
@click.option('--parent', 'parent_specs',
              required=True, multiple=True,
              callback=_validate_parent,
              metavar='NAME:hap1.fa,hap2.fa',
              help='亲本配置（可重复指定多个亲本）'
                   '|Parent spec (can be repeated for multiple parents). '
                   '格式|Format: NAME:hap1.fa,hap2.fa')
@click.option('-o', '--output-dir',
              default='./subgenome_assign_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--preset',
              type=click.Choice(['asm5', 'asm10', 'asm20', 'asm25']),
              default='asm10',
              show_default=True,
              help='minimap2 -x 预设|minimap2 preset (asm5=<1%%, asm10=1-5%%, asm20=5-15%%)')
@click.option('-t', '--threads',
              type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--minimap2-secondary',
              is_flag=True,
              help='保留次要比对（默认 --secondary=no）|Keep secondary alignments')
@click.option('--min-conf',
              type=float, default=0.65, show_default=True,
              help='置信度阈值|Confidence threshold for LOW_CONFIDENCE flag')
@click.option('--no-split',
              is_flag=True,
              help='不输出拆分的 FASTA|Do not output split FASTAs')
@click.option('--no-keep-unassigned',
              is_flag=True,
              help='不输出未归属染色体的 FASTA|Do not output unassigned FASTA')
@click.option('--minimap2-path',
              default='~/miniforge3/envs/cphasing/bin/minimap2', show_default=True,
              help='minimap2 二进制路径|minimap2 binary path')
@click.option('--samtools-path',
              default='~/.local/bin/samtools', show_default=True,
              help='samtools 二进制路径|samtools binary path')
def subgenome_assign(target, parent_specs, output_dir, preset, threads,
                     minimap2_secondary, min_conf,
                     no_split, no_keep_unassigned,
                     minimap2_path, samtools_path):
    """
    亚基因组归属工具|Subgenome Assignment Tool

    通过将目标多倍体基因组比对到各亲本参考，按比对覆盖度判定每条染色体的亚基因组来源
    |Assign each chromosome to a subgenome by alignment against parental references

    示例|Example:
        biopytools subgenome_assign -i Cf.chr.fa \\
            --parent Ca:Ca_hap1.fa,Ca_hap2.fa \\
            --parent Ch:Ch_hap1.fa,Ch_hap2.fa \\
            -o results -t 16
    """

    subgenome_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['subgenome_assign.py']
    args.extend(['-i', target])

    # --parent 必须全部传给 main.py 的 action='append'|Pass all --parent specs
    for spec in parent_specs:
        args.extend(['--parent', spec])

    # 非默认值才传|Pass non-defaults only
    if output_dir != './subgenome_assign_output':
        args.extend(['-o', output_dir])
    if preset != 'asm10':
        args.extend(['--preset', preset])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if minimap2_secondary:
        args.append('--minimap2-secondary')
    if min_conf != 0.65:
        args.extend(['--min-conf', str(min_conf)])
    if no_split:
        args.append('--no-split')
    if no_keep_unassigned:
        args.append('--no-keep-unassigned')
    if minimap2_path != '~/miniforge3/envs/cphasing/bin/minimap2':
        args.extend(['--minimap2-path', minimap2_path])
    if samtools_path != '~/.local/bin/samtools':
        args.extend(['--samtools-path', samtools_path])

    original_argv = sys.argv
    sys.argv = args

    try:
        subgenome_main()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
