"""
NGenomeSyn Command
NGenomeSyn基因组共线性分析命令
"""

import click
import sys
import os


def _lazy_import_ngenomesyn_main():
    """懒加载ngenomesyn main函数 | Lazy load ngenomesyn main function"""
    try:
        from ...ngenomesyn.main import main as ngenomesyn_main
        return ngenomesyn_main
    except ImportError as e:
        click.echo(f"[ERROR] 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="NGenomeSyn基因组共线性分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--sample-map', '-s',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='[FILE] 样本映射文件 | Sample mapping file (genome_file\\tgenome_name)')
@click.option('--config', '-c',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='[FILE] 配置文件 | Configuration file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='[DIR] 输出目录 | Output directory')
@click.option('--aligner', '-a',
              default='minimap2',
              type=click.Choice(['minimap2', 'mummer']),
              help='[STR] 比对器类型 (默认: minimap2) | Aligner type (default: minimap2)')
@click.option('--threads', '-t',
              default=16,
              type=int,
              help='[INT] 线程数 (默认: 16) | Number of threads (default: 16)')
@click.option('--min-length',
              default=5000,
              type=int,
              help='[INT] 最小比对长度 (默认: 5000) | Minimum alignment length (default: 5000)')
@click.option('--minimap-preset',
              default='asm5',
              help='[STR] Minimap2预设模式 (默认: asm5) | Minimap2 preset (default: asm5)')
@click.option('--mummer-match-type',
              default='mumreference',
              type=click.Choice(['mum', 'mumreference', 'maxmatch']),
              help='[STR] MUMmer匹配类型 (默认: mumreference) | MUMmer match type (default: mumreference)')
@click.option('--mummer-min-match',
              default=20,
              type=int,
              help='[INT] MUMmer最小匹配长度 (默认: 20) | MUMmer min match (default: 20)')
@click.option('--chromosome',
              type=str,
              help='[STR] 指定染色体 | Specify chromosomes (e.g., "1,2,3" or "1-5")')
@click.option('--output-formats',
              multiple=True,
              default=['svg', 'png'],
              type=click.Choice(['svg', 'png']),
              help='[STR] 输出格式 (默认: svg png) | Output formats (default: svg png)')
@click.option('--ngenomesyn-bin',
              help='[FILE] NGenomeSyn二进制文件路径 | NGenomeSyn binary path')
@click.option('--use-syri',
              is_flag=True,
              default=False,
              help='[FLAG] 使用SyRI进行结构变异分析 | Use SyRI for structural variation analysis')
@click.option('--syri-bin',
              help='[FILE] SyRI二进制文件路径 | SyRI binary path')
def ngenomesyn(sample_map, config, output, aligner, threads, min_length,
               minimap_preset, mummer_match_type, mummer_min_match,
               chromosome, output_formats, ngenomesyn_bin, use_syri, syri_bin):
    """
    NGenomeSyn基因组共线性分析工具

    从基因组FASTA文件进行共线性分析并生成可视化图形。
    支持Minimap2和MUMmer两种比对器，自动完成比对和可视化。
    可选使用SyRI进行结构变异分析。

    示例 | Examples:

    基本用法:
      biopytools ngenomesyn -s samples.txt -o output_dir

    使用MUMmer比对器:
      biopytools ngenomesyn -s samples.txt -o output_dir --aligner mummer

    使用SyRI进行结构变异分析:
      biopytools ngenomesyn -s samples.txt -o output_dir --use-syri

    分析特定染色体:
      biopytools ngenomesyn -s samples.txt -o output_dir --chromosome "1,2,3"

    自定义参数:
      biopytools ngenomesyn -s samples.txt -o output_dir \\
          --threads 32 --min-length 10000

    sample_map文件格式 | sample_map File Format:

    格式为两列，用Tab分隔：
      genome_file.fa    GenomeName
      reference.fa     Ref
      query.fa         Query

    输出说明 | Output Description:

    工具会自动生成以下文件:
      - .len文件: 基因组长度信息
      - .paf/.link文件: 比对结果
      - .svg/.png文件: 可视化图形
      - .conf文件: NGenomeSyn配置文件
      - 使用--use-syri时还会生成:
        - .syri.out: SyRI分析结果
        - .syri.vcf: 结构变异VCF文件

    注意事项 | Notes:

      - 至少需要2个基因组文件
      - 支持FASTA和GZIP压缩的FASTA格式
      - Minimap2速度较快，MUMmer更精确
      - 染色体编号从1开始，连续编号
      - 使用--use-syri需要安装SyRI程序
    """
    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    ngenomesyn_main = _lazy_import_ngenomesyn_main()

    # 验证输入参数的互斥性
    if not sample_map and not config:
        raise click.ClickException("必须指定 --sample-map 或 --config 参数之一 | Must specify either --sample-map or --config")

    if sample_map and config:
        raise click.ClickException("--sample-map 和 --config 参数不能同时使用 | Cannot use both --sample-map and --config")

    # 构建参数列表传递给原始main函数
    args = ['ngenomesyn.py']

    # 输入文件参数（互斥）
    if sample_map:
        args.extend(['-s', sample_map])
    if config:
        args.extend(['-c', config])

    # 输出参数
    args.extend(['-o', output])

    # 比对参数（只在非默认值时添加）
    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if min_length != 5000:
        args.extend(['--min-length', str(min_length)])

    # 比对器特定参数
    if minimap_preset != 'asm5':
        args.extend(['--minimap-preset', minimap_preset])

    if mummer_match_type != 'mumreference':
        args.extend(['--mummer-match-type', mummer_match_type])

    if mummer_min_match != 20:
        args.extend(['--mummer-min-match', str(mummer_min_match)])

    # 染色体过滤参数
    if chromosome:
        args.extend(['--chromosome', chromosome])

    # 可视化参数
    if output_formats and set(output_formats) != {'svg', 'png'}:
        args.extend(['--output-formats'] + list(output_formats))

    if ngenomesyn_bin:
        args.extend(['--ngenomesyn-bin', ngenomesyn_bin])

    # SyRI参数
    if use_syri:
        args.append('--use-syri')

    if syri_bin:
        args.extend(['--syri-bin', syri_bin])

    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数
        ngenomesyn_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n[WARNING] 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"[ERROR] 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
