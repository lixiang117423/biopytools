"""
BRAKER3基因组注释|BRAKER3 Genome Annotation Command
"""

import click
import sys
import os


def _lazy_import_braker_main():
    """延迟加载braker主函数|Lazy load braker main function"""
    try:
        # 修复：从 main.py 导入 parse_arguments 和其他需要的函数
        import sys
        from ...braker import main as braker_module
        return braker_module.main
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
    short_help='BRAKER3基因组注释工具|BRAKER3 genome annotation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# 必需参数|Required parameters
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')

@click.option('--species', '-s',
              required=True,
              help='物种名称|Species name (for BRAKER output naming)')

# 输入数据|Input data
@click.option('--prot_seq', '-p',
              help='近缘物种蛋白质序列文件或文件夹|Protein sequences file or directory')

@click.option('--isoseq', '-l',
              help='三代全长转录本文件夹|Long-read transcript directory')

@click.option('--rnaseq_dirs',
              help='二代RNA-seq目录列表(逗号分隔)|Comma-separated RNA-seq directories')

@click.option('--read1_pattern',
              default='_1.clean.fq.gz',
              show_default=True,
              help='R1文件模式|R1 file pattern')

@click.option('--read2_pattern',
              default='_2.clean.fq.gz',
              show_default=True,
              help='R2文件模式|R2 file pattern')

# 输出配置|Output configuration
@click.option('--output_dir', '-o',
              default='./braker_output',
              show_default=True,
              help='输出目录|Output directory')

# 流程参数|Pipeline parameters
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')

@click.option('--fungus',
              is_flag=True,
              help='使用真菌模式|Use fungus mode (suitable for oomycetes)')

# Singularity配置|Singularity configuration
@click.option('--singularity_image',
              default="~/software/singularity/braker3_devel.sif",
              show_default=True,
              help='Singularity镜像路径|Singularity image path')

@click.option('--no_singularity',
              is_flag=True,
              help='不使用Singularity镜像|Do not use Singularity image')

# 步骤控制|Step control
@click.option('--skip_repeat',
              is_flag=True,
              help='跳过重复序列屏蔽|Skip repeat masking')

@click.option('--skip_long_reads',
              is_flag=True,
              help='跳过三代转录本处理|Skip long-read processing')

@click.option('--skip_short_reads',
              is_flag=True,
              help='跳过二代RNA-seq处理|Skip short-read processing')

# BRAKER3特定参数|BRAKER3 specific parameters
@click.option('--busco_lineage',
              help='BUSCO谱系|BUSCO lineage')

@click.option('--utr',
              is_flag=True,
              help='预测UTR|Predict UTR')

@click.option('--training_genes',
              help='训练基因集|Training gene set file')

@click.option('--use_existing',
              is_flag=True,
              help='使用已有参数|Use existing parameters')

def braker(genome, species, prot_seq, isoseq, rnaseq_dirs,
           output_dir, threads, fungus, singularity_image, no_singularity,
           skip_repeat, skip_long_reads, skip_short_reads,
           busco_lineage, utr, training_genes, use_existing,
           read1_pattern, read2_pattern):
    """
    BRAKER3基因组注释工具|BRAKER3 Genome Annotation Tool

    基于BRAKER3、GeneMark-ETP和AUGUSTUS的真核生物基因组基因结构注释流程|Eukaryotic genome annotation pipeline based on BRAKER3, GeneMark-ETP and AUGUSTUS

    示例|Example: biopytools braker -g genome.fa -s my_oomycete -p proteins.fa
    """

    # 延迟加载|Lazy loading
    braker_main = _lazy_import_braker_main()

    # 构建参数列表|Build argument list
    args = ['braker.py']

    # 必需参数|Required parameters
    args.extend(['--genome', genome])
    args.extend(['--species', species])

    # 输入数据|Input data
    if prot_seq:
        args.extend(['--prot_seq', prot_seq])

    if isoseq:
        args.extend(['--isoseq', isoseq])

    if rnaseq_dirs:
        args.extend(['--rnaseq_dirs', rnaseq_dirs])

    if read1_pattern != '_1.clean.fq.gz':
        args.extend(['--read1_pattern', read1_pattern])

    if read2_pattern != '_2.clean.fq.gz':
        args.extend(['--read2_pattern', read2_pattern])

    # 输出配置|Output configuration
    if output_dir != './braker_output':
        args.extend(['--output_dir', output_dir])

    # 流程参数|Pipeline parameters
    if threads != 40:
        args.extend(['--threads', str(threads)])

    if fungus:
        args.append('--fungus')

    # Singularity配置|Singularity configuration
    if singularity_image != "~/software/singularity/braker3_devel.sif":
        args.extend(['--singularity_image', singularity_image])

    if no_singularity:
        args.append('--no_singularity')

    # 步骤控制|Step control
    if skip_repeat:
        args.append('--skip_repeat')

    if skip_long_reads:
        args.append('--skip_long_reads')

    if skip_short_reads:
        args.append('--skip_short_reads')

    # BRAKER3特定参数|BRAKER3 specific parameters
    if busco_lineage:
        args.extend(['--busco_lineage', busco_lineage])

    if utr:
        args.append('--utr')

    if training_genes:
        args.extend(['--training_genes', training_genes])

    if use_existing:
        args.append('--use_existing')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        braker_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
