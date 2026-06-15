"""
HiFiasm基因组组装CLI包装器|HiFiasm Genome Assembly CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_hifiasm_main():
    """延迟加载hifiasm主函数|Lazy load hifiasm main function"""
    try:
        from ...hifiasm.main import main as hifiasm_main
        return hifiasm_main
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


def _validate_required_file_exists(file_path):
    """验证必需文件存在(仅在非帮助模式)|Validate required file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help='HiFiasm基因组组装|HiFiasm Genome Assembly',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-reads', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_required_file_exists(value) if value else None,
              help='HiFi测序数据文件|Input HiFi sequencing data file')
@click.option('--output-dir', '-o',
              default='./hifiasm_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              default='sample',
              type=str,
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--hg-size',
              default='auto',
              type=str,
              show_default=True,
              help='基因组大小估计(如1.4g, 2.1g)|Genome size estimation (e.g., 1.4g, 2.1g)')
@click.option('--purge-level', '-l',
              default=3,
              type=int,
              show_default=True,
              help='Purge级别(0-3)|Purge level (0-3)')
@click.option('--purge-max',
              default=65,
              type=int,
              show_default=True,
              help='最大purge覆盖度|Maximum purge coverage')
@click.option('--similarity-threshold', '-s',
              default=0.75,
              type=float,
              show_default=True,
              help='相似性阈值|Similarity threshold')
@click.option('--ont-reads',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='ONT长读长数据文件|ONT long-read data file')
@click.option('--hi-c-1',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C第一端数据文件|Hi-C first-end data file')
@click.option('--hi-c-2',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Hi-C第二端数据文件|Hi-C second-end data file')
@click.option('--extra-hifiasm-args',
              default='',
              type=str,
              show_default=True,
              help='额外的HiFiasm参数|Additional HiFiasm arguments')
@click.option('--skip-busco',
              is_flag=True,
              help='跳过BUSCO质量评估|Skip BUSCO quality assessment')
@click.option('--busco-lineage',
              default='auto',
              type=str,
              show_default=True,
              help='BUSCO谱系数据集|BUSCO lineage dataset')
@click.option('--busco-mode',
              default='genome',
              type=click.Choice(['genome', 'proteins', 'transcriptome']),
              show_default=True,
              help='BUSCO评估模式|BUSCO assessment mode')
@click.option('--skip-quast',
              is_flag=True,
              help='跳过QUAST质量评估|Skip QUAST quality assessment')
@click.option('--reference-genome',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组文件(用于QUAST)|Reference genome file (for QUAST)')
@click.option('--analyze-haplotypes',
              is_flag=True,
              help='分析单倍型差异|Analyze haplotype differences')
@click.option('--min-contig-length',
              default=1000,
              type=int,
              show_default=True,
              help='最小contig长度过滤|Minimum contig length filter')
@click.option('--generate-plots',
              is_flag=True,
              help='生成可视化图表|Generate visualization plots')
@click.option('--assembly-type',
              default='auto',
              type=click.Choice(['auto', 'diploid', 'triploid', 'polyploid']),
              show_default=True,
              help='组装类型|Assembly type')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--compress-output',
              is_flag=True,
              help='压缩输出文件|Compress output files')
@click.option('--output-formats',
              multiple=True,
              type=click.Choice(['fasta', 'gfa', 'both']),
              default=['both'],
              show_default=True,
              help='输出格式选择|Output format selection')
@click.option('--memory',
              default=100,
              type=int,
              show_default=True,
              help='内存大小(GB)|Memory size (GB)')
@click.option('--tmp-dir',
              default='/tmp',
              type=click.Path(),
              show_default=True,
              help='临时目录|Temporary directory')
@click.option('--max-runtime',
              default=48,
              type=int,
              show_default=True,
              help='最大运行时间(小时)|Maximum runtime (hours)')
@click.option('--resume',
              is_flag=True,
              help='恢复中断的分析|Resume interrupted analysis')
@click.option('--hifiasm-path',
              default='hifiasm',
              type=str,
              show_default=True,
              help='HiFiasm软件路径|HiFiasm software path')
@click.option('--busco-path',
              default='busco',
              type=str,
              show_default=True,
              help='BUSCO软件路径|BUSCO software path')
@click.option('--quast-path',
              default='quast',
              type=str,
              show_default=True,
              help='QUAST软件路径|QUAST software path')
@click.option('--python-path',
              default='python3',
              type=str,
              show_default=True,
              help='Python解释器路径|Python interpreter path')
@click.option('--samtools-path',
              default='samtools',
              type=str,
              show_default=True,
              help='Samtools软件路径|Samtools software path')
@click.option('--busco-db-path',
              type=click.Path(),
              help='BUSCO数据库路径|BUSCO database path')
@click.option('--busco-download-path',
              type=click.Path(),
              help='BUSCO数据集下载路径|BUSCO dataset download path')
@click.option('--debug',
              is_flag=True,
              help='启用调试模式|Enable debug mode')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-level',
              default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              show_default=True,
              help='日志级别|Log level')
@click.option('--config-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='配置文件路径|Configuration file path')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式(不执行实际命令)|Dry run mode (do not execute actual commands)')
def hifiasm(input_reads, output_dir, prefix, threads, hg_size, purge_level,
            purge_max, similarity_threshold, ont_reads, hi_c_1, hi_c_2,
            extra_hifiasm_args, skip_busco, busco_lineage, busco_mode,
            skip_quast, reference_genome, analyze_haplotypes, min_contig_length,
            generate_plots, assembly_type, keep_intermediate, compress_output,
            output_formats, memory, tmp_dir, max_runtime, resume, hifiasm_path,
            busco_path, quast_path, python_path, samtools_path, busco_db_path,
            busco_download_path, debug, verbose, log_level, config_file, dry_run):
    """
    HiFiasm基因组组装流程|HiFiasm Genome Assembly Pipeline

    基于HiFi reads的基因组组装自动化流程
    Automated genome assembly pipeline based on HiFi reads

    示例|Examples: biopytools hifiasm -i sample.hifi.fq.gz -o hifiasm_results -p sample_prefix
    """

    # 延迟加载|Lazy load
    hifiasm_main = _lazy_import_hifiasm_main()

    # 验证Hi-C数据|Validate Hi-C data
    if (hi_c_1 and not hi_c_2) or (hi_c_2 and not hi_c_1):
        raise click.ClickException("Hi-C数据需要同时提供两端文件|Hi-C data requires both end files")

    # 验证purge级别|Validate purge level
    if purge_level not in range(0, 4):
        raise click.ClickException("Purge级别必须在0-3之间|Purge level must be between 0-3")

    # 验证相似性阈值|Validate similarity threshold
    if similarity_threshold <= 0 or similarity_threshold > 1:
        raise click.ClickException("相似性阈值必须在0-1之间|Similarity threshold must be between 0-1")

    # 构建主函数参数列表|Build argument list for main function
    args = ['hifiasm.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_reads])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output_dir != './hifiasm_output':
        args.extend(['-o', output_dir])

    if prefix != 'sample':
        args.extend(['-p', prefix])

    if threads != 32:
        args.extend(['-t', str(threads)])

    # HiFiasm参数|HiFiasm parameters
    if hg_size != 'auto':
        args.extend(['--hg-size', hg_size])

    if purge_level != 3:
        args.extend(['-l', str(purge_level)])

    if purge_max != 65:
        args.extend(['--purge-max', str(purge_max)])

    if similarity_threshold != 0.75:
        args.extend(['-s', str(similarity_threshold)])

    if ont_reads:
        args.extend(['--ont-reads', ont_reads])

    if hi_c_1:
        args.extend(['--hi-c-1', hi_c_1])

    if hi_c_2:
        args.extend(['--hi-c-2', hi_c_2])

    if extra_hifiasm_args:
        args.extend(['--extra-hifiasm-args', extra_hifiasm_args])

    # 质量评估参数|Quality assessment parameters
    if skip_busco:
        args.append('--skip-busco')

    if busco_lineage != 'auto':
        args.extend(['--busco-lineage', busco_lineage])

    if busco_mode != 'genome':
        args.extend(['--busco-mode', busco_mode])

    if skip_quast:
        args.append('--skip-quast')

    if reference_genome:
        args.extend(['--reference-genome', reference_genome])

    # 分析参数|Analysis parameters
    if analyze_haplotypes:
        args.append('--analyze-haplotypes')

    if min_contig_length != 1000:
        args.extend(['--min-contig-length', str(min_contig_length)])

    if generate_plots:
        args.append('--generate-plots')

    if assembly_type != 'auto':
        args.extend(['--assembly-type', assembly_type])

    # 输出控制参数|Output control parameters
    if keep_intermediate:
        args.append('--keep-intermediate')

    if compress_output:
        args.append('--compress-output')

    if output_formats and set(output_formats) != {'both'}:
        args.extend(['--output-formats'] + list(output_formats))

    # 系统参数|System parameters
    if memory != 64:
        args.extend(['--memory', str(memory)])

    if tmp_dir != '/tmp':
        args.extend(['--tmp-dir', tmp_dir])

    if max_runtime != 48:
        args.extend(['--max-runtime', str(max_runtime)])

    if resume:
        args.append('--resume')

    # 工具路径参数|Tool path parameters
    if hifiasm_path != 'hifiasm':
        args.extend(['--hifiasm-path', hifiasm_path])

    if busco_path != 'busco':
        args.extend(['--busco-path', busco_path])

    if quast_path != 'quast':
        args.extend(['--quast-path', quast_path])

    if python_path != 'python3':
        args.extend(['--python-path', python_path])

    if samtools_path != 'samtools':
        args.extend(['--samtools-path', samtools_path])

    # 数据库路径参数|Database path parameters
    if busco_db_path:
        args.extend(['--busco-db-path', busco_db_path])

    if busco_download_path:
        args.extend(['--busco-download-path', busco_download_path])

    # 高级参数|Advanced parameters
    if debug:
        args.append('--debug')

    if verbose > 0:
        args.extend(['-' + 'v' * verbose])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    if config_file:
        args.extend(['--config-file', config_file])

    if dry_run:
        args.append('--dry-run')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        hifiasm_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"意外错误|Unexpected error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
