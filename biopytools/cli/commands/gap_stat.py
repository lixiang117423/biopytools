"""
基因组Gap统计命令|Genome Gap Statistics Command
"""

import click
import sys
import os


def _lazy_import_gapstat_main():
    """延迟加载gapstat主函数|Lazy load gapstat main function"""
    try:
        from ...gapstat.main import main as gapstat_main
        return gapstat_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_input_exists(input_path):
    """验证输入存在(文件或文件夹，仅在非帮助模式)|Validate input exists (file or directory, only in non-help mode)"""
    help_flags = {'-h', '--help'}
    if not any(arg in help_flags for arg in sys.argv) and input_path and not os.path.exists(input_path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {input_path}")
    return input_path


def _get_fasta_files_from_dir(input_dir):
    """从文件夹获取所有FASTA文件|Get all FASTA files from directory"""
    fasta_files = []
    fasta_extensions = ['.fa', '.fasta', '.FA', '.FASTA',
                       '.fna', '.fnasta', '.FNA', '.FNASTA']

    for filename in os.listdir(input_dir):
        for ext in fasta_extensions:
            if filename.endswith(ext):
                fasta_files.append(os.path.join(input_dir, filename))
                break

    return sorted(fasta_files)


@click.command(
    short_help='基因组Gap统计工具|Genome gap statistics tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_input_exists(value),
              help='输入FASTA文件或文件夹|Input FASTA file or directory')
@click.option('-o', '--output',
              help='输出文件路径|Output file path')
@click.option('--min-n',
              type=int,
              default=1,
              show_default=True,
              help='最少N数量|Minimum consecutive N count')
def gap_stat(input, output, min_n):
    """
    基因组Gap统计工具|Genome Gap Statistics Tool

    统计基因组FASTA文件中Gap的位置和长度，支持批量处理|Statistics gap positions and lengths in genome FASTA files with batch processing support

    示例|Examples: biopytools gap-stat -i genome.fa -o gaps.txt
    """

    # 检查输入是文件还是文件夹|Check if input is file or directory
    if os.path.isfile(input):
        # 单文件模式|Single file mode
        click.echo(f"处理单个FASTA文件|Processing single FASTA file: {input}")

        # 延迟加载|Lazy loading
        gapstat_main = _lazy_import_gapstat_main()

        # 构建参数列表|Build argument list
        args = ['gapstat.py']
        args.extend([input])

        if output:
            args.extend(['-o', output])

        if min_n != 1:
            args.extend(['--min-n', str(min_n)])

        # 执行主程序|Execute main program
        original_argv = sys.argv
        sys.argv = args

        try:
            gapstat_main()
        except SystemExit as e:
            sys.exit(e.code)
        except Exception as e:
            click.echo(f"错误|Error: {e}", err=True)
            sys.exit(1)
        finally:
            sys.argv = original_argv

    elif os.path.isdir(input):
        # 批量处理模式 - 合并所有FASTA文件到一个输出文件|Batch processing mode
        click.echo(f"批量处理模式|Batch processing mode")
        click.echo(f"输入文件夹|Input directory: {input}")
        click.echo(f"输出文件|Output file: {output}")

        # 获取所有FASTA文件|Get all FASTA files
        fasta_files = _get_fasta_files_from_dir(input)

        if not fasta_files:
            click.echo(f"错误|Error: 在文件夹中未找到FASTA文件|No FASTA files found in directory: {input}", err=True)
            sys.exit(1)

        click.echo(f"找到|Found {len(fasta_files)} 个FASTA文件|FASTA files")
        click.echo(f"合并输出到一个文件|Merging output to single file")

        # 使用批量合并处理|Use batch merge processing
        from ...gapstat.main import batch_process_gapstat

        try:
            success = batch_process_gapstat(
                fasta_files=fasta_files,
                output_file=output,
                min_n=min_n
            )

            if success:
                click.echo(f"\n✓ 批量处理成功|Batch processing completed successfully")
                click.echo(f"输出文件|Output file: {output}")
                sys.exit(0)
            else:
                click.echo(f"\n✗ 批量处理失败|Batch processing failed", err=True)
                sys.exit(1)

        except Exception as e:
            click.echo(f"✗ 错误|Error: {str(e)}", err=True)
            import traceback
            click.echo(traceback.format_exc(), err=True)
            sys.exit(1)

    else:
        click.echo(f"错误|Error: 输入路径既不是文件也不是文件夹|Input path is neither a file nor a directory: {input}", err=True)
        sys.exit(1)
