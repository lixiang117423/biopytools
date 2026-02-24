"""
BUSCO质量评估分析命令|BUSCO Quality Assessment Analysis Command
"""

import click
import sys
import os


def _lazy_import_busco_main():
    """延迟加载BUSCO主函数|Lazy load BUSCO main function"""
    try:
        from ...busco.main import main as busco_main
        return busco_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='BUSCO质量评估|BUSCO quality assessment',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入文件或目录|Input file or directory')
@click.option('--lineage', '-l',
              required=True,
              type=str,
              help='BUSCO数据库谱系|BUSCO database lineage')
@click.option('--output-dir', '-o',
              default='./busco_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--mode', '-m',
              type=click.Choice(['genome', 'geno', 'transcriptome', 'tran', 'proteins', 'prot']),
              default='genome',
              show_default=True,
              help='BUSCO分析模式|BUSCO analysis mode')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='CPU线程数|Number of CPU threads')
@click.option('--sample-suffix',
              default='*.fa',
              show_default=True,
              help='样本名提取后缀|Sample name extraction suffix pattern')
@click.option('--output-format',
              type=click.Choice(['txt', 'csv', 'xlsx']),
              default='txt',
              show_default=True,
              help='输出文件格式|Output file format')
@click.option('--force', '-f',
              is_flag=True,
              help='强制重写现有文件|Force rewriting existing files')
@click.option('--augustus',
              is_flag=True,
              help='使用Augustus基因预测器|Use Augustus gene predictor')
@click.option('--augustus-parameters',
              type=str,
              help='Augustus额外参数|Additional Augustus parameters')
@click.option('--augustus-species',
              type=str,
              help='Augustus物种名|Augustus species name')
@click.option('--auto-lineage',
              is_flag=True,
              help='自动选择谱系|Automatically select lineage')
@click.option('--auto-lineage-euk',
              is_flag=True,
              help='自动选择真核生物谱系|Automatically select eukaryote lineage')
@click.option('--auto-lineage-prok',
              is_flag=True,
              help='自动选择原核生物谱系|Automatically select prokaryote lineage')
@click.option('--contig-break',
              type=int,
              default=10,
              show_default=True,
              help='Contig打断的N数量|Number of Ns for contig break')
@click.option('--datasets-version',
              default='odb12',
              show_default=True,
              help='数据集版本|Dataset version')
@click.option('--download-path',
              type=str,
              help='数据集下载路径|Dataset download path')
@click.option('--evalue', '-e',
              type=float,
              default=1e-3,
              show_default=True,
              help='BLAST E值阈值|BLAST E-value threshold')
@click.option('--limit',
              type=int,
              default=3,
              show_default=True,
              help='候选区域限制|Candidate region limit')
@click.option('--long',
              is_flag=True,
              help='启用Augustus长模式优化|Enable Augustus long mode optimization')
@click.option('--metaeuk',
              is_flag=True,
              help='使用Metaeuk基因预测器|Use Metaeuk gene predictor')
@click.option('--metaeuk-parameters',
              type=str,
              help='Metaeuk额外参数|Additional Metaeuk parameters')
@click.option('--metaeuk-rerun-parameters',
              type=str,
              help='Metaeuk重新运行参数|Metaeuk rerun parameters')
@click.option('--miniprot',
              is_flag=True,
              help='使用Miniprot基因预测器|Use Miniprot gene predictor')
@click.option('--skip-bbtools',
              is_flag=True,
              help='跳过BBTools统计|Skip BBTools statistics')
@click.option('--offline',
              is_flag=True,
              help='离线模式|Offline mode')
@click.option('--restart', '-r',
              is_flag=True,
              help='重启未完成的分析|Restart incomplete analysis')
@click.option('--quiet', '-q',
              is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--scaffold-composition',
              is_flag=True,
              help='生成scaffold组成文件|Generate scaffold composition file')
@click.option('--tar',
              is_flag=True,
              help='压缩子目录|Compress subdirectories')
@click.option('--busco-path',
              default='busco',
              show_default=True,
              help='BUSCO软件路径|BUSCO software path')
def busco(input, lineage, output_dir, mode, threads, sample_suffix, output_format,
                  force, augustus, augustus_parameters, augustus_species, auto_lineage,
                  auto_lineage_euk, auto_lineage_prok, contig_break, datasets_version,
                  download_path, evalue, limit, long, metaeuk, metaeuk_parameters,
                  metaeuk_rerun_parameters, miniprot, skip_bbtools, offline, restart,
                  quiet, scaffold_composition, tar, busco_path):
    """
    BUSCO质量评估分析工具|BUSCO Quality Assessment Analysis Tool

    基于单拷贝直系同源基因评估基因组完整性|Assess genome completeness using single-copy orthologs

    示例|Examples: biopytools busco -i genome.fa -l eukaryota_odb12
    """

    # Lazy loading: import only when actually called
    busco_main = _lazy_import_busco_main()

    # 构建参数列表|Build argument list for original main function
    args = ['busco_analysis.py']

    # Required parameters
    args.extend(['-i', input])
    args.extend(['-l', lineage])

    #Optional parameters (add only when non-default)
    if output_dir != './busco_output':
        args.extend(['-o', output_dir])

    if mode != 'genome':
        args.extend(['-m', mode])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if sample_suffix != '*.fa':
        args.extend(['--sample-suffix', sample_suffix])

    if output_format != 'txt':
        args.extend(['--output-format', output_format])

    if contig_break != 10:
        args.extend(['--contig-break', str(contig_break)])

    if datasets_version != 'odb12':
        args.extend(['--datasets-version', datasets_version])

    if evalue != 1e-3:
        args.extend(['-e', str(evalue)])

    if limit != 3:
        args.extend(['--limit', str(limit)])

    if busco_path != 'busco':
        args.extend(['--busco-path', busco_path])

    # String parameters
    if augustus_parameters:
        args.extend(['--augustus-parameters', augustus_parameters])

    if augustus_species:
        args.extend(['--augustus-species', augustus_species])

    if download_path:
        args.extend(['--download-path', download_path])

    if metaeuk_parameters:
        args.extend(['--metaeuk-parameters', metaeuk_parameters])

    if metaeuk_rerun_parameters:
        args.extend(['--metaeuk-rerun-parameters', metaeuk_rerun_parameters])

    # Boolean options
    if force:
        args.append('-f')

    if augustus:
        args.append('--augustus')

    if auto_lineage:
        args.append('--auto-lineage')

    if auto_lineage_euk:
        args.append('--auto-lineage-euk')

    if auto_lineage_prok:
        args.append('--auto-lineage-prok')

    if long:
        args.append('--long')

    if metaeuk:
        args.append('--metaeuk')

    if miniprot:
        args.append('--miniprot')

    if skip_bbtools:
        args.append('--skip-bbtools')

    if offline:
        args.append('--offline')

    if restart:
        args.append('-r')

    if quiet:
        args.append('-q')

    if scaffold_composition:
        args.append('--scaffold-composition')

    if tar:
        args.append('--tar')

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        busco_main()
    except SystemExit as e:
        # Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断BUSCO分析|BUSCO quality assessment analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"BUSCO分析失败|BUSCO quality assessment analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
