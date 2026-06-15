"""
叶绿体基因组组装命令|Plastome Assembly Command
"""

import click
import sys
import os


def _lazy_import_plastome_classes():
    """延迟加载叶绿体组装类|Lazy load plastome classes"""
    try:
        from ...get_plastome.main import PlastomeAssembler
        return PlastomeAssembler
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_dir_exists(dir_path):
    """验证目录存在性(仅在非帮助模式下)|Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.isdir(dir_path):
        raise click.BadParameter(f"Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='叶绿体基因组组装|Plastome Assembly',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='输入目录(包含reads文件)|Input directory containing reads files')
@click.option('--output-dir', '-o',
              default='./plastome_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              help='输出前缀|Output prefix')
@click.option('--organelle-type',
              default='embplant_pt',
              show_default=True,
              type=click.Choice(['embplant_pt', 'embplant_mt', 'embplant_nr', 'other_pt',
                                'animal_mt', 'fungus_mt', 'fungus_nr']),
              help='细胞器类型|Organelle type')
@click.option('--max-rounds', '-R',
              type=int,
              default=15,
              show_default=True,
              help='最大扩展轮数|Maximum extension rounds')
@click.option('--kmer-list', '-k',
              default='21,45,65,85,105',
              show_default=True,
              help='Kmer列表(逗号分隔)|Kmer list comma-separated')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Threads')
@click.option('--getorganelle-path',
              default='~/miniforge3/envs/getorganelle_v.1.7.71/bin/get_organelle_from_reads.py',
              show_default=True,
              help='GetOrganelle脚本路径|GetOrganelle script path')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-file',
              help='日志文件路径|Log file path')
@click.option('--single-mode',
              is_flag=True,
              hidden=True,
              help='单样品模式|Single sample mode')
def get_plastome(input, output_dir, prefix, organelle_type, max_rounds,
                 kmer_list, threads, getorganelle_path, verbose, log_file, single_mode):
    """
    叶绿体基因组组装工具|Plastome Assembly Tool

    基于GetOrganelle从二代全基因组测序数据中自动组装叶绿体基因组|Assemble plastome from WGS data using GetOrganelle

    示例|Examples: biopytools get-plastome -i fastq_folder -o plastome_output
    """

    # 延迟加载|Lazy loading: import only when actually called
    PlastomeAssembler = _lazy_import_plastome_classes()

    # 默认使用批量模式，除非明确指定--single-mode|Use batch mode by default, unless --single-mode is explicitly specified
    if not single_mode:
        # 批量模式|Batch mode
        kwargs = {
            'getorganelle_path': getorganelle_path,
            'organelle_type': organelle_type,
            'max_rounds': max_rounds,
            'kmer_list': kmer_list,
            'threads': threads,
        }

        results = PlastomeAssembler.run_batch(input, output_dir, **kwargs)

        # 根据结果设置退出码|Set exit code based on results
        if results:
            successful = sum(1 for v in results.values() if v)
            failed = len(results) - successful
            if failed > 0:
                sys.exit(1)  # 有失败的样品|Some samples failed
            else:
                sys.exit(0)  # 全部成功|All successful
        else:
            sys.exit(1)  # 没有样品|No samples
    else:
        # 单个样品模式|Single sample mode
        # 延迟加载|Lazy loading: import only when actually called
        from ...get_plastome.main import main as plastome_main

        # 构建参数列表|Build argument list
        args = ['plastome_assembler.py']

        # 必需参数|Required parameters
        args.extend(['--input', input])

        # 可选参数|Optional parameters (add only when non-default)
        if output_dir != './plastome_output':
            args.extend(['--output-dir', output_dir])

        if prefix:
            args.extend(['--prefix', prefix])

        # GetOrganelle参数|GetOrganelle parameters
        if organelle_type != 'embplant_pt':
            args.extend(['--organelle-type', organelle_type])

        if max_rounds != 15:
            args.extend(['--max-rounds', str(max_rounds)])

        if kmer_list != '21,45,65,85,105':
            args.extend(['--kmer-list', kmer_list])

        if threads != 12:
            args.extend(['--threads', str(threads)])

        # 软件配置|Software configuration
        if getorganelle_path != '~/miniforge3/envs/getorganelle_v.1.7.71/bin/get_organelle_from_reads.py':
            args.extend(['--getorganelle-path', getorganelle_path])

        # 日志参数|Logging parameters
        if verbose:
            args.append('--verbose')

        if log_file:
            args.extend(['--log-file', log_file])

        # 保存和恢复sys.argv|Save and restore sys.argv
        original_argv = sys.argv
        sys.argv = args

        try:
            # 调用原始主函数|Call original main function
            plastome_main()
        except SystemExit as e:
            # 处理程序正常退出|Handle normal program exit
            if e.code != 0:
                sys.exit(e.code)
        except KeyboardInterrupt:
            click.echo("\n用户中断|Interrupted by user", err=True)
            sys.exit(1)
        except Exception as e:
            click.echo(f"执行失败|Execution failed: {e}", err=True)
            sys.exit(1)
        finally:
            sys.argv = original_argv
