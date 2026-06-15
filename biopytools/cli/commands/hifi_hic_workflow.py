"""
HiFi+Hi-C工作流CLI命令|HiFi+Hi-C Workflow CLI Command

biopytools hifi-hic-workflow命令的Click包装器
Click wrapper for biopytools hifi-hic-workflow command
"""

import click
import sys
import os
from pathlib import Path


@click.command(
    name='hifi-hic-workflow',
    help='HiFi+Hi-C基因组组装与挂载流程|HiFi+Hi-C Genome Assembly and Scaffolding Workflow'
)
@click.option(
    '--hifi',
    'hifi_reads',
    required=True,
    type=click.Path(exists=True),
    help='HiFi reads文件|HiFi reads file'
)
@click.option(
    '--hic-r1',
    required=True,
    type=click.Path(exists=True),
    help='Hi-C R1文件|Hi-C R1 file'
)
@click.option(
    '--hic-r2',
    required=True,
    type=click.Path(exists=True),
    help='Hi-C R2文件|Hi-C R2 file'
)
@click.option(
    '--ref', '--reference',
    'reference_genome',
    required=True,
    type=click.Path(exists=True),
    help='参考基因组FASTA文件（仅用于命名）|Reference genome FASTA file (for naming only)'
)
@click.option(
    '-o', '--output',
    required=True,
    type=click.Path(),
    help='输出目录|Output directory'
)
@click.option(
    '-p', '--prefix',
    default='genome_sample',
    help='样本前缀|Sample prefix (default: genome_sample)'
)
@click.option(
    '-t', '--threads',
    default=64,
    type=int,
    help='线程数|Number of threads (default: 64)'
)
@click.option(
    '--use-ngs-polish',
    is_flag=True,
    help='启用NGS polish|Enable NGS polish'
)
@click.option(
    '--ngs-data',
    type=click.Path(exists=True),
    help='NGS二代数据目录|NGS second-generation data directory'
)
@click.option(
    '--nchrs',
    type=int,
    help='染色体数量（如不指定，从reference统计）|Number of chromosomes (count from reference if not specified)'
)
@click.option(
    '--skip-hifi-hic',
    is_flag=True,
    help='跳过HiFi组装|Skip HiFi assembly'
)
@click.option(
    '--skip-haphic',
    is_flag=True,
    help='跳过Hi-C挂载|Skip Hi-C scaffolding'
)
@click.option(
    '--skip-rename',
    is_flag=True,
    help='跳过重命名|Skip renaming'
)
@click.option(
    '--skip-heatmap',
    is_flag=True,
    help='跳过热图|Skip heatmap'
)
@click.option(
    '--no-resume',
    is_flag=True,
    help='禁用断点续传|Disable resume mode'
)
@click.option(
    '--force',
    is_flag=True,
    help='强制重新运行所有步骤|Force rerun all steps'
)
@click.option(
    '-v', '--verbose',
    is_flag=True,
    help='显示详细日志|Show verbose logs'
)
def hifi_hic_workflow(
    hifi_reads, hic_r1, hic_r2, reference_genome, output, prefix,
    threads, use_ngs_polish, ngs_data, nchrs,
    skip_hifi_hic, skip_haphic, skip_rename, skip_heatmap,
    no_resume, force, verbose
):
    """
    HiFi+Hi-C工作流命令|HiFi+Hi-C Workflow Command

    完整的植物基因组组装流程：HiFi组装 → Hi-C挂载 → 染色体重命名 → Hi-C热图
    Complete plant genome assembly workflow: HiFi assembly → Hi-C scaffolding → Chromosome rename → Hi-C heatmap
    """
    # 延迟导入主模块|Lazy import main module
    try:
        from biopytools.hifi_hic_workflow.main import main as workflow_main

        # 构建参数列表|Build argument list
        sys.argv = [
            'hifi_hic_workflow',
            '--hifi', hifi_reads,
            '--hic-r1', hic_r1,
            '--hic-r2', hic_r2,
            '--ref', reference_genome,
            '-o', output,
            '-p', prefix,
            '-t', str(threads),
        ]

        # NGS polish|NGS polish
        if use_ngs_polish:
            sys.argv.append('--use-ngs-polish')
            if ngs_data:
                sys.argv.extend(['--ngs-data', ngs_data])

        # 染色体数|Chromosomes
        if nchrs:
            sys.argv.extend(['--nchrs', str(nchrs)])

        # 流程控制|Workflow control
        if skip_hifi_hic:
            sys.argv.append('--skip-hifi-hic')
        if skip_haphic:
            sys.argv.append('--skip-haphic')
        if skip_rename:
            sys.argv.append('--skip-rename')
        if skip_heatmap:
            sys.argv.append('--skip-heatmap')
        if no_resume:
            sys.argv.append('--no-resume')
        if force:
            sys.argv.append('--force')

        # 详细日志|Verbose logging
        if verbose:
            sys.argv.append('--verbose')

        # 运行主程序|Run main program
        exit_code = workflow_main()
        sys.exit(exit_code)

    except ImportError as e:
        click.echo(f"错误|Error: 无法导入hifi_hic_workflow模块|Cannot import hifi_hic_workflow module: {e}", err=True)
        click.echo("请确保biopytools已正确安装|Please ensure biopytools is properly installed", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
