"""
YaHS Hi-C Scaffolding命令|YaHS Hi-C Scaffolding Command
"""

import click
import sys
import os


def _lazy_import_yahs_main():
    """延迟加载yahs主函数|Lazy load yahs main function"""
    try:
        from ...yahs.main import main as yahs_main
        return yahs_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
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
    short_help='YaHS Hi-C scaffolding流程|YaHS Hi-C scaffolding pipeline',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# 必需参数|Required parameters
@click.option('-r', '--ref',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='参考基因组FASTA|Reference genome FASTA')
@click.option('-1', '--hic-r1',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='Hi-C R1文件|Hi-C R1 file')
@click.option('-2', '--hic-r2',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='Hi-C R2文件|Hi-C R2 file')

# 输出配置|Output configuration
@click.option('-o', '--output-dir',
              default='./yahs_output',
              show_default=True,
              help='输出目录|Output directory')

# 资源配置|Resource configuration
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--java-ram',
              default='300G',
              show_default=True,
              help='Java内存|Java memory')
@click.option('--sam-ram',
              default='300G',
              show_default=True,
              help='Samtools排序内存|Samtools sort memory')

# YaHS 核心参数|YaHS core parameters
@click.option('-e', '--enzyme',
              default='GATC',
              show_default=True,
              help='限制性酶切位点|Restriction enzyme sequence')
@click.option('--min-len',
              type=int,
              default=10000,
              show_default=True,
              help='最小contig长度|Minimum contig length')
@click.option('--min-mapq',
              type=int,
              default=30,
              show_default=True,
              help='最小MAPQ值|Minimum MAPQ value')
@click.option('--no-contig-ec',
              is_flag=True,
              help='跳过contig错误校正|Skip contig error correction')
@click.option('--no-scaffold-ec',
              is_flag=True,
              help='跳过scaffold错误校正|Skip scaffold error correction')
@click.option('--resolutions',
              help='分辨率列表(逗号分隔)|Resolution list (comma-separated)')
@click.option('--rounds',
              type=int,
              default=1,
              show_default=True,
              help='每分辨率运行轮数|Rounds per resolution')
@click.option('--telo-motif',
              help='端粒序列模体|Telomeric sequence motif')

# 工具路径|Tool paths
@click.option('--yahs-bin',
              default='~/miniforge3/envs/yahs_v.1.2.2/bin/yahs',
              show_default=True,
              help='YaHS可执行文件|YaHS executable')
@click.option('--juicer-bin',
              default='~/miniforge3/envs/yahs_v.1.2.2/bin/juicer',
              show_default=True,
              help='juicer可执行文件|juicer executable')
@click.option('--juicer-jar',
              default='~/software/juicer/scripts/juicer_tools.jar',
              show_default=True,
              help='juicer_tools.jar文件|juicer_tools.jar file')
@click.option('--bwa-bin',
              default='~/miniforge3/envs/Population_genetics/bin/bwa',
              show_default=True,
              help='BWA可执行文件|BWA executable')
@click.option('--samtools-bin',
              default='~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools',
              show_default=True,
              help='samtools可执行文件|samtools executable')
@click.option('--java-cmd',
              default='java',
              show_default=True,
              help='Java可执行文件|Java executable')

# 执行控制|Execution control
@click.option('-s', '--step',
              type=click.Choice(['1', '2', '3', '4', '5', '6']),
              help='运行指定步骤|Run specified step')
@click.option('--force-rerun',
              is_flag=True,
              help='强制重新运行|Force rerun all steps')
@click.option('--keep-temp',
              is_flag=True,
              help='保留临时文件|Keep temporary files')

def yahs(ref, hic_r1, hic_r2, output_dir, threads, java_ram, sam_ram,
         enzyme, min_len, min_mapq, no_contig_ec, no_scaffold_ec,
         resolutions, rounds, telo_motif,
         yahs_bin, juicer_bin, juicer_jar, bwa_bin, samtools_bin, java_cmd,
         step, force_rerun, keep_temp):
    """
    YaHS Hi-C scaffolding流程|YaHS Hi-C scaffolding pipeline

    使用YaHS进行Hi-C scaffolding分析，支持完整流程和单步执行
    Perform Hi-C scaffolding using YaHS, supports complete pipeline and single-step execution

    示例|Examples: biopytools yahs -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz
    """

    # 延迟加载|Lazy loading
    yahs_main = _lazy_import_yahs_main()

    # 构建参数列表|Build argument list
    args = ['yahs.py']

    # 必需参数|Required parameters
    args.extend(['-r', ref])
    args.extend(['-1', hic_r1])
    args.extend(['-2', hic_r2])

    # 输出配置|Output configuration
    if output_dir != './yahs_output':
        args.extend(['-o', output_dir])

    # 资源配置|Resource configuration
    if threads != 12:
        args.extend(['-t', str(threads)])
    if java_ram != '300G':
        args.extend(['--java-ram', java_ram])
    if sam_ram != '300G':
        args.extend(['--sam-ram', sam_ram])

    # YaHS参数|YaHS parameters
    if enzyme != 'GATC':
        args.extend(['-e', enzyme])
    if min_len != 10000:
        args.extend(['--min-len', str(min_len)])
    if min_mapq != 30:
        args.extend(['--min-mapq', str(min_mapq)])
    if no_contig_ec:
        args.append('--no-contig-ec')
    if no_scaffold_ec:
        args.append('--no-scaffold-ec')
    if resolutions:
        args.extend(['--resolutions', resolutions])
    if rounds != 1:
        args.extend(['--rounds', str(rounds)])
    if telo_motif:
        args.extend(['--telo-motif', telo_motif])

    # 工具路径|Tool paths
    if yahs_bin != '~/miniforge3/envs/yahs_v.1.2.2/bin/yahs':
        args.extend(['--yahs-bin', yahs_bin])
    if juicer_bin != '~/miniforge3/envs/yahs_v.1.2.2/bin/juicer':
        args.extend(['--juicer-bin', juicer_bin])
    if juicer_jar != '~/software/juicer/scripts/juicer_tools.jar':
        args.extend(['--juicer-jar', juicer_jar])
    if bwa_bin != '~/miniforge3/envs/Population_genetics/bin/bwa':
        args.extend(['--bwa-bin', bwa_bin])
    if samtools_bin != '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools':
        args.extend(['--samtools-bin', samtools_bin])
    if java_cmd != 'java':
        args.extend(['--java-cmd', java_cmd])

    # 执行控制|Execution control
    if step:
        args.extend(['-s', step])
    if force_rerun:
        args.append('--force-rerun')
    if keep_temp:
        args.append('--keep-temp')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        yahs_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
