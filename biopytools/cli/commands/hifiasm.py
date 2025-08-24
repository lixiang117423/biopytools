"""
HiFiasm基因组组装分析命令 | HiFiasm Genome Assembly Analysis Command
"""

import click
import sys
from ...hifiasm.main import main as hifiasm_main


@click.command(short_help = '运行hifiasm组装流程',
               context_settings=dict(help_option_names=['-h', '--help'],max_content_width=120))
@click.option('--input-reads', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入HiFi测序数据文件 | Input HiFi sequencing data file')
@click.option('--output-dir', '-o',
              default='./hifiasm_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./hifiasm_output)')
@click.option('--prefix', '-p',
              default='sample',
              type=str,
              help='输出文件前缀 | Output file prefix (default: sample)')
@click.option('--threads', '-t',
              default=32,
              type=int,
              help='线程数 | Number of threads (default: 32)')
@click.option('--hg-size',
              default='auto',
              type=str,
              help='基因组大小估计 (如 1.4g, 2.1g) | Genome size estimation (e.g., 1.4g, 2.1g)')
@click.option('--purge-level', '-l',
              default=3,
              type=int,
              help='purge级别 (0-3) | Purge level (0-3) (default: 3)')
@click.option('--purge-max',
              default=65,
              type=int,
              help='最大purge覆盖度 | Maximum purge coverage (default: 65)')
@click.option('--similarity-threshold', '-s',
              default=0.75,
              type=float,
              help='相似性阈值 | Similarity threshold (default: 0.75)')
@click.option('--ont-reads',
              type=click.Path(exists=True),
              help='ONT长读长数据文件 | ONT long-read data file')
@click.option('--hi-c-1',
              type=click.Path(exists=True),
              help='Hi-C第一端数据文件 | Hi-C first-end data file')
@click.option('--hi-c-2',
              type=click.Path(exists=True),
              help='Hi-C第二端数据文件 | Hi-C second-end data file')
@click.option('--extra-hifiasm-args',
              default='',
              type=str,
              help='额外的HiFiasm参数 | Additional HiFiasm arguments')
@click.option('--skip-busco',
              is_flag=True,
              help='跳过BUSCO质量评估 | Skip BUSCO quality assessment')
@click.option('--busco-lineage',
              default='auto',
              type=str,
              help='BUSCO谱系数据集 (如 embryophyta_odb10) | BUSCO lineage dataset (default: auto)')
@click.option('--busco-mode',
              default='genome',
              type=click.Choice(['genome', 'proteins', 'transcriptome']),
              help='BUSCO评估模式 | BUSCO assessment mode (default: genome)')
@click.option('--skip-quast',
              is_flag=True,
              help='跳过QUAST质量评估 | Skip QUAST quality assessment')
@click.option('--reference-genome',
              type=click.Path(exists=True),
              help='参考基因组文件 (用于QUAST) | Reference genome file (for QUAST)')
@click.option('--analyze-haplotypes',
              is_flag=True,
              help='分析单倍型差异 | Analyze haplotype differences')
@click.option('--min-contig-length',
              default=1000,
              type=int,
              help='最小contig长度过滤 | Minimum contig length filter (default: 1000)')
@click.option('--generate-plots',
              is_flag=True,
              help='生成可视化图表 | Generate visualization plots')
@click.option('--assembly-type',
              default='auto',
              type=click.Choice(['auto', 'diploid', 'triploid', 'polyploid']),
              help='组装类型 | Assembly type (default: auto)')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件 | Keep intermediate files')
@click.option('--compress-output',
              is_flag=True,
              help='压缩输出文件 | Compress output files')
@click.option('--output-formats',
              multiple=True,
              type=click.Choice(['fasta', 'gfa', 'both']),
              default=['both'],
              help='输出格式选择 | Output format selection (default: both)')
@click.option('--memory',
              default=64,
              type=int,
              help='内存大小(GB) | Memory size (GB) (default: 64)')
@click.option('--tmp-dir',
              default='/tmp',
              type=click.Path(),
              help='临时目录 | Temporary directory (default: /tmp)')
@click.option('--max-runtime',
              default=48,
              type=int,
              help='最大运行时间(小时) | Maximum runtime (hours) (default: 48)')
@click.option('--resume',
              is_flag=True,
              help='恢复中断的分析 | Resume interrupted analysis')
@click.option('--hifiasm-path',
              default='hifiasm',
              type=str,
              help='HiFiasm软件路径 | HiFiasm software path (default: hifiasm)')
@click.option('--busco-path',
              default='busco',
              type=str,
              help='BUSCO软件路径 | BUSCO software path (default: busco)')
@click.option('--quast-path',
              default='quast',
              type=str,
              help='QUAST软件路径 | QUAST software path (default: quast)')
@click.option('--python-path',
              default='python3',
              type=str,
              help='Python解释器路径 | Python interpreter path (default: python3)')
@click.option('--samtools-path',
              default='samtools',
              type=str,
              help='Samtools软件路径 | Samtools software path (default: samtools)')
@click.option('--busco-db-path',
              type=click.Path(),
              help='BUSCO数据库路径 | BUSCO database path')
@click.option('--busco-download-path',
              type=click.Path(),
              help='BUSCO数据集下载路径 | BUSCO dataset download path')
@click.option('--debug',
              is_flag=True,
              help='启用调试模式 | Enable debug mode')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式 (-v, -vv, -vvv) | Verbose output mode')
@click.option('--log-level',
              default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='日志级别 | Log level (default: INFO)')
@click.option('--config-file',
              type=click.Path(exists=True),
              help='配置文件路径 | Configuration file path')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式 (不执行实际命令) | Dry run mode (do not execute actual commands)')
def hifiasm(input_reads, output_dir, prefix, threads, hg_size, purge_level, 
            purge_max, similarity_threshold, ont_reads, hi_c_1, hi_c_2, 
            extra_hifiasm_args, skip_busco, busco_lineage, busco_mode, 
            skip_quast, reference_genome, analyze_haplotypes, min_contig_length,
            generate_plots, assembly_type, keep_intermediate, compress_output,
            output_formats, memory, tmp_dir, max_runtime, resume, hifiasm_path,
            busco_path, quast_path, python_path, samtools_path, busco_db_path,
            busco_download_path, debug, verbose, log_level, config_file, dry_run):
    """
    HiFiasm基因组组装完整流水线
    
    使用HiFiasm进行高质量基因组组装，包括组装、质量评估、
    单倍型分析等完整流程。支持二倍体和多倍体基因组组装。
    
    示例 | Examples:
    
    \b
    # 基本用法
    biopytools hifiasm -i sample.hifi.fq.gz -o hifiasm_results -p sample_prefix
    
    \b
    # 二倍体组装
    biopytools hifiasm -i reads.hifi.fq.gz -o output \\
        --hg-size 1.4g --purge-max 65 -s 0.75
    
    \b
    # 包含ONT数据的组装
    biopytools hifiasm -i hifi.fq.gz --ont-reads ont.fq.gz \\
        -o output -p sample
    
    \b
    # Hi-C辅助组装
    biopytools hifiasm -i hifi.fq.gz \\
        --hi-c-1 hic_R1.fq.gz --hi-c-2 hic_R2.fq.gz \\
        -o output -p sample
    
    \b
    # 完整质量评估
    biopytools hifiasm -i sample.hifi.fq.gz -o output \\
        --busco-lineage embryophyta_odb10 \\
        --reference-genome ref.fa --analyze-haplotypes
    
    \b
    # 高内存高线程组装
    biopytools hifiasm -i large.hifi.fq.gz -o output \\
        -t 128 --memory 256 --hg-size 5g
    
    \b
    # 试运行模式检查参数
    biopytools hifiasm -i sample.hifi.fq.gz -o output \\
        --dry-run --debug -vvv
    
    \b
    # 自定义工具路径
    biopytools hifiasm -i sample.hifi.fq.gz -o output \\
        --hifiasm-path /opt/hifiasm/hifiasm \\
        --busco-path /opt/busco/bin/busco \\
        --quast-path /opt/quast/quast.py
    
    \b
    # 多倍体组装
    biopytools hifiasm -i polyploid.hifi.fq.gz -o output \\
        --assembly-type polyploid --purge-level 1 \\
        --min-contig-length 10000
    """
    
    # 验证Hi-C数据配对
    if (hi_c_1 and not hi_c_2) or (hi_c_2 and not hi_c_1):
        raise click.ClickException("Hi-C数据需要同时提供两端数据 | Hi-C data requires both end files")
    
    # 验证purge级别
    if purge_level not in range(0, 4):
        raise click.ClickException("purge级别必须在0-3之间 | Purge level must be between 0-3")
    
    # 验证相似性阈值
    if similarity_threshold <= 0 or similarity_threshold > 1:
        raise click.ClickException("相似性阈值必须在0-1之间 | Similarity threshold must be between 0-1")
    
    # 构建参数列表传递给原始main函数
    args = ['hifiasm.py']
    
    # 必需参数
    args.extend(['-i', input_reads])
    
    # 基本参数（只在非默认值时添加）
    if output_dir != './hifiasm_output':
        args.extend(['-o', output_dir])
    
    if prefix != 'sample':
        args.extend(['-p', prefix])
    
    if threads != 32:
        args.extend(['-t', str(threads)])
    
    # HiFiasm组装参数
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
    
    # 质量评估参数
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
    
    # 分析参数
    if analyze_haplotypes:
        args.append('--analyze-haplotypes')
    
    if min_contig_length != 1000:
        args.extend(['--min-contig-length', str(min_contig_length)])
    
    if generate_plots:
        args.append('--generate-plots')
    
    if assembly_type != 'auto':
        args.extend(['--assembly-type', assembly_type])
    
    # 输出控制参数
    if keep_intermediate:
        args.append('--keep-intermediate')
    
    if compress_output:
        args.append('--compress-output')
    
    if output_formats and set(output_formats) != {'both'}:
        args.extend(['--output-formats'] + list(output_formats))
    
    # 系统参数
    if memory != 64:
        args.extend(['--memory', str(memory)])
    
    if tmp_dir != '/tmp':
        args.extend(['--tmp-dir', tmp_dir])
    
    if max_runtime != 48:
        args.extend(['--max-runtime', str(max_runtime)])
    
    if resume:
        args.append('--resume')
    
    # 工具路径参数
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
    
    # 数据库路径参数
    if busco_db_path:
        args.extend(['--busco-db-path', busco_db_path])
    
    if busco_download_path:
        args.extend(['--busco-download-path', busco_download_path])
    
    # 高级参数
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
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        hifiasm_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"意外错误 | Unexpected error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv