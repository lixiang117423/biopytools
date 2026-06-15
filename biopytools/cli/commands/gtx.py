"""
GTX WGS分析命令|GTX WGS Analysis Command
"""

import click
import sys
import os


def _lazy_import_gtx_main():
    """延迟加载gtx主函数|Lazy load gtx main function"""
    try:
        from ...gtx.main import main as gtx_main
        return gtx_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_dir(dir_path):
    """验证输入目录存在性(仅在非帮助模式)|Validate input directory existence (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(dir_path):
            raise click.BadParameter(f"输入目录不存在|Input directory does not exist: {dir_path}")
        if not os.path.isdir(dir_path):
            raise click.BadParameter(f"输入路径不是目录|Input path is not a directory: {dir_path}")
    return dir_path


def _validate_reference_file(file_path):
    """验证参考基因组文件存在性(仅在非帮助模式)|Validate reference file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"参考基因组文件不存在|Reference genome file does not exist: {file_path}")
    return file_path


def _validate_joint_output(joint_output):
    """验证Joint输出文件名格式|Validate joint output filename format"""
    if joint_output and not joint_output.endswith('.vcf.gz'):
        raise click.BadParameter(f"Joint输出文件必须以.vcf.gz结尾|Joint output file must end with .vcf.gz: {joint_output}")
    return joint_output


@click.command(short_help='GTX WGS分析|GTX WGS Analysis',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_dir(value) if value else None,
              help='输入目录(包含clean FASTQ文件)|Input directory path (containing clean FASTQ files)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录路径|Output directory path')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_reference_file(value) if value else None,
              help='参考基因组文件路径|Reference genome file path')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--gtx-path',
              default='~/software/gtx/GTX.CAT_2.2.1/bin/gtx',
              type=click.Path(),
              show_default=True,
              help='GTX程序路径|GTX program path')
@click.option('--tmp-dir',
              type=click.Path(),
              help='临时目录路径|Temporary directory path')
@click.option('--enable-joint',
              is_flag=True,
              help='启用Joint Calling|Enable Joint Calling for multi-sample variant detection')
@click.option('--joint-output',
              default='merged_gtx.vcf.gz',
              show_default=True,
              callback=lambda ctx, param, value: _validate_joint_output(value) if value else None,
              help='Joint calling输出VCF文件名|Joint calling output VCF filename')
@click.option('--joint-threads',
              default=88,
              type=int,
              show_default=True,
              help='Joint calling使用的线程数|Number of threads for joint calling')
@click.option('--min-confidence',
              default=30,
              type=int,
              show_default=True,
              help='最小置信度阈值|Minimum confidence threshold')
@click.option('--min-base-quality',
              default=20,
              type=int,
              show_default=True,
              help='最小碱基质量阈值|Minimum base quality threshold')
@click.option('--ploidy',
              default=2,
              type=int,
              show_default=True,
              help='倍性|Ploidy')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              show_default=True,
              type=click.Choice(['CONSERVATIVE', 'AGGRESSIVE']),
              help='PCR indel模型|PCR indel model')
@click.option('--read1-pattern',
              default='*_1.clean.fq.gz',
              show_default=True,
              type=str,
              help='R1文件匹配模式|R1 file pattern')
@click.option('--read2-pattern',
              default='*_2.clean.fq.gz',
              show_default=True,
              type=str,
              help='R2文件匹配模式|R2 file pattern')
def gtx(input_dir, output_dir, reference, threads, gtx_path, tmp_dir,
            enable_joint, joint_output, joint_threads, min_confidence, min_base_quality,
            ploidy, pcr_indel_model, read1_pattern, read2_pattern):
    """
    GTX WGS分析工具|GTX WGS Analysis Tool

    基于GTX进行WGS数据分析，支持单样品和Joint Calling模式
    Perform WGS data analysis using GTX, supports individual and joint calling modes

    示例|Examples: biopytools gtx -i /data/clean_fastq -o /results -r /genome/reference.fa --enable-joint
    """

    # 延迟加载|Lazy load: import only when actually called
    gtx_main = _lazy_import_gtx_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['gtx_wgs.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])

    # 可选参数|Optional parameters
    if threads != 88:
        args.extend(['-t', str(threads)])

    if gtx_path != '~/software/gtx/GTX.CAT_2.2.1/bin/gtx':
        args.extend(['--gtx-path', gtx_path])

    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])

    # Joint calling参数|Joint calling parameters
    if enable_joint:
        args.append('--enable-joint')

    if joint_output != 'merged_gtx.vcf.gz':
        args.extend(['--joint-output', joint_output])

    if joint_threads != 88:
        args.extend(['--joint-threads', str(joint_threads)])

    # 质量控制参数|Quality control parameters
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])

    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])

    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])

    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])

    # 文件模式参数|File pattern parameters
    if read1_pattern != '*_1.clean.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])

    if read2_pattern != '*_2.clean.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gtx_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
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
