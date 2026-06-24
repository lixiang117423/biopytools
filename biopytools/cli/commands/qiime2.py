"""
QIIME2微生物组多样性分析命令|QIIME2 Microbiome Diversity Analysis Command
"""

import click
import sys
import os


def _lazy_import_qiime2_main():
    """延迟加载qiime2主函数|Lazy load qiime2 main function"""
    try:
        from ...qiime2.main import main as qiime2_main
        return qiime2_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.isdir(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='QIIME2微生物组多样性分析|QIIME2 Microbiome Diversity Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-dir',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='双端FASTQ输入目录|Input directory of paired-end FASTQ')
@click.option('-o', '--output-dir',
              default='./qiime2_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--amplicon',
              type=click.Choice(['16s', 'its']),
              default='16s',
              show_default=True,
              help='扩增子类型|Amplicon type')
@click.option('--method',
              type=click.Choice(['asv', 'otu']),
              default='asv',
              show_default=True,
              help='聚类方法ASV(DADA2)或OTU(vsearch)|Method: ASV or OTU')
@click.option('--fwd-primer',
              default='CCTACGGGNGGCWGCAG',
              show_default=True,
              help='正向引物序列(IUPAC)|Forward primer')
@click.option('--rev-primer',
              default='GACTACHVGGGTATCTAATCC',
              show_default=True,
              help='反向引物序列(IUPAC)|Reverse primer')
@click.option('--trunc-len-f', type=int, default=0, show_default=True,
              help='R1截断长度(0=不截断)|R1 truncation length (0=none)')
@click.option('--trunc-len-r', type=int, default=0, show_default=True,
              help='R2截断长度(0=不截断)|R2 truncation length (0=none)')
@click.option('--trim-left-f', type=int, default=0, show_default=True,
              help='R1左侧裁剪|R1 trim left')
@click.option('--trim-left-r', type=int, default=0, show_default=True,
              help='R2左侧裁剪|R2 trim left')
@click.option('--sampling-depth', type=int, default=0, show_default=True,
              help='抽平深度(0=自动)|Rarefaction depth (0=auto)')
@click.option('--perc-identity', type=float, default=0.97, show_default=True,
              help='OTU聚类相似度|OTU identity')
@click.option('--confidence', type=float, default=0.7, show_default=True,
              help='分类置信度|Classification confidence')
@click.option('--min-length', type=int, default=50, show_default=True,
              help='extract-reads最小长度|extract-reads min length')
@click.option('--max-length', type=int, default=0, show_default=True,
              help='extract-reads最大长度(0=不限)|extract-reads max length (0=none)')
@click.option('-t', '--threads', type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--validate-level',
              type=click.Choice(['min', 'max']),
              default='min', show_default=True,
              help='tools import校验级别|import validate level')
@click.option('--classifier',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='预训练分类器(.qza),省略则自动训练|Pre-trained classifier (.qza)')
@click.option('--database-dir',
              help='原始参考库目录(SILVA/UNITE)|Raw reference DB directory')
@click.option('--qiime-path',
              help='qiime可执行文件路径|qiime executable path')
@click.option('--classifier-cache-dir',
              help='分类器缓存目录|Classifier cache directory')
@click.option('--r1-suffix', default='_1.clean.fq.gz', show_default=True,
              help='R1文件后缀|R1 file suffix')
@click.option('--r2-suffix', default='_2.clean.fq.gz', show_default=True,
              help='R2文件后缀|R2 file suffix')
@click.option('--skip-cutadapt', is_flag=True,
              help='跳过引物切除|Skip primer trimming')
@click.option('--skip-phylogeny', is_flag=True,
              help='跳过系统发育建树(ITS自动跳过)|Skip phylogeny (auto for ITS)')
@click.option('--metadata-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='样品元数据TSV(可选)|Sample metadata TSV (optional)')
@click.option('-f', '--force', is_flag=True,
              help='覆盖已有输出|Overwrite existing outputs')
@click.option('-v', '--verbose', is_flag=True,
              help='详细输出|Verbose output')
def qiime2(input_dir, output_dir, amplicon, method, fwd_primer, rev_primer,
           trunc_len_f, trunc_len_r, trim_left_f, trim_left_r, sampling_depth,
           perc_identity, confidence, min_length, max_length, threads,
           validate_level, classifier, database_dir, qiime_path,
           classifier_cache_dir, r1_suffix, r2_suffix, skip_cutadapt,
           skip_phylogeny, metadata_file, force, verbose):
    """
    QIIME2微生物组多样性分析|QIIME2 Microbiome Diversity Analysis

    自动识别双端样品并完成QIIME2全流程(ASV/OTU、丰度表、分类注释、多样性、抽平)|
    Auto-detect paired reads and run full QIIME2 pipeline

    示例|Example: biopytools qiime2 -i raw_reads/ -o qiime2_output
    """

    qiime2_main = _lazy_import_qiime2_main()

    args = ['qiime2.py']
    args.extend(['-i', input_dir])

    if output_dir != './qiime2_output':
        args.extend(['-o', output_dir])
    if amplicon != '16s':
        args.extend(['--amplicon', amplicon])
    if method != 'asv':
        args.extend(['--method', method])
    if fwd_primer != 'CCTACGGGNGGCWGCAG':
        args.extend(['--fwd-primer', fwd_primer])
    if rev_primer != 'GACTACHVGGGTATCTAATCC':
        args.extend(['--rev-primer', rev_primer])
    if trunc_len_f != 0:
        args.extend(['--trunc-len-f', str(trunc_len_f)])
    if trunc_len_r != 0:
        args.extend(['--trunc-len-r', str(trunc_len_r)])
    if trim_left_f != 0:
        args.extend(['--trim-left-f', str(trim_left_f)])
    if trim_left_r != 0:
        args.extend(['--trim-left-r', str(trim_left_r)])
    if sampling_depth != 0:
        args.extend(['--sampling-depth', str(sampling_depth)])
    if perc_identity != 0.97:
        args.extend(['--perc-identity', str(perc_identity)])
    if confidence != 0.7:
        args.extend(['--confidence', str(confidence)])
    if min_length != 50:
        args.extend(['--min-length', str(min_length)])
    if max_length != 0:
        args.extend(['--max-length', str(max_length)])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if validate_level != 'min':
        args.extend(['--validate-level', validate_level])
    if classifier:
        args.extend(['--classifier', classifier])
    if database_dir:
        args.extend(['--database-dir', database_dir])
    if qiime_path:
        args.extend(['--qiime-path', qiime_path])
    if classifier_cache_dir:
        args.extend(['--classifier-cache-dir', classifier_cache_dir])
    if r1_suffix != '_1.clean.fq.gz':
        args.extend(['--r1-suffix', r1_suffix])
    if r2_suffix != '_2.clean.fq.gz':
        args.extend(['--r2-suffix', r2_suffix])
    if skip_cutadapt:
        args.append('--skip-cutadapt')
    if skip_phylogeny:
        args.append('--skip-phylogeny')
    if metadata_file:
        args.extend(['--metadata-file', metadata_file])
    if force:
        args.append('-f')
    if verbose:
        args.append('-v')

    original_argv = sys.argv
    sys.argv = args

    try:
        qiime2_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
