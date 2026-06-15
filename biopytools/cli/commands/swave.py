"""
Swave结构变异检测命令|Swave Structural Variant Detection Command
"""

import click
import sys
import os
import tempfile


def _lazy_import_main():
    """延迟加载swave主函数|Lazy load swave main function"""
    try:
        from ...swave.main import main as swave_main
        return swave_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    help_flags = {'-h', '--help'}
    is_help = any(arg in help_flags for arg in sys.argv)

    if not is_help and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


def _absolutize_tsv_paths(tsv_path):
    """将TSV文件中的相对路径转为绝对路径|Convert relative paths in TSV to absolute paths

    Swave内部通过os.system调用外部工具，cwd为swave_path而非用户工作目录，
    TSV中的相对路径会因cwd不同而解析失败。此函数基于TSV文件所在目录将所有
    相对路径转为绝对路径，并返回处理后的临时文件路径。

    Args:
        tsv_path: TSV文件路径|TSV file path

    Returns:
        处理后的TSV文件路径（临时文件）|Processed TSV file path (temp file)
    """
    tsv_dir = os.path.dirname(os.path.abspath(tsv_path))

    with open(tsv_path, 'r') as f:
        lines = f.readlines()

    need_convert = False
    processed_lines = []
    for line in lines:
        if 'NAME' in line:
            processed_lines.append(line)
            continue
        parts = line.strip().split()
        if len(parts) < 2:
            processed_lines.append(line)
            continue
        new_parts = [parts[0]]
        for p in parts[1:]:
            if not os.path.isabs(p):
                new_parts.append(os.path.abspath(os.path.join(tsv_dir, p)))
                need_convert = True
            else:
                new_parts.append(p)
        processed_lines.append('\t'.join(new_parts) + '\n')

    if not need_convert:
        return tsv_path

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False, dir=tsv_dir)
    tmp.writelines(processed_lines)
    tmp.close()
    return tmp.name


@click.group(
    short_help='Swave结构变异检测工具|Swave structural variant detection tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def swave():
    """
    Swave结构变异检测工具|Swave Structural Variant Detection Tool

    从泛基因组图检测结构变异和复杂SV|Detect structural variants and complex SVs from pangenome graphs
    """
    pass


@swave.command(short_help='检测结构变异|Call structural variants')
@click.option('-i', '--assemblies-tsv',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='样本组装TSV文件|Assemblies TSV file')
@click.option('-r', '--ref-fasta',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('-g', '--gfa-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='泛基因组图GFA文件|Pangenome graph GFA file')
@click.option('-s', '--gfa-source',
              required=True,
              type=click.Choice(['minigraph', 'cactus', 'pggb']),
              help='GFA文件来源|GFA file source (minigraph/cactus/pggb)')
@click.option('--swave-path',
              default='~/software/swave/Swave-main',
              show_default=True,
              help='Swave软件路径|Swave software path')
@click.option('-o', '--output-dir',
              default='./swave_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--decomposed-vcf',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Decomposed VCF文件(cactus/pggb必需)|Decomposed VCF file (required for cactus/pggb)')
@click.option('--output-mode',
              default='auto',
              type=click.Choice(['auto', 'population', 'single']),
              show_default=True,
              help='输出模式|Output mode')
@click.option('--spec-samples',
              multiple=True,
              help='指定样本|Specify samples')
@click.option('--min-sv-size',
              type=int,
              default=50,
              show_default=True,
              help='最小SV大小|Minimum SV size')
@click.option('--max-sv-size',
              type=int,
              default=1000000,
              show_default=True,
              help='最大SV大小|Maximum SV size')
@click.option('--max-sv-comps',
              type=int,
              default=5,
              show_default=True,
              help='最大SV组件数|Maximum number of SV components')
@click.option('--dup-to-ins',
              is_flag=True,
              help='将duplication报告为insertion|Report duplications as insertions')
@click.option('--remove-small',
              is_flag=True,
              help='移除小于min_sv_size的节点|Remove nodes smaller than min_sv_size')
@click.option('--force-reverse',
              is_flag=True,
              help='强制调用反向映射snarls|Force call reversed mapping snarls')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--minigraph-path',
              default='minigraph',
              show_default=True,
              help='minigraph工具路径|minigraph tool path')
@click.option('--gfatools-path',
              default='gfatools',
              show_default=True,
              help='gfatools工具路径|gfatools tool path')
@click.option('--spec-snarl',
              help='只调用特定snarl|Only call specific snarl')
@click.option('--spec-path',
              help='只调用特定path|Only call specific path')
def call(assemblies_tsv, ref_fasta, gfa_file, gfa_source, swave_path,
         output_dir, decomposed_vcf, output_mode, spec_samples,
         min_sv_size, max_sv_size, max_sv_comps,
         dup_to_ins, remove_small, force_reverse, threads,
         minigraph_path, gfatools_path, spec_snarl, spec_path):
    """
    检测结构变异|Call structural variants from inputs

    示例|Example: biopytools swave call -i assemblies.tsv -r ref.fa -g graph.gfa -s minigraph
    """
    swave_main = _lazy_import_main()

    # 自动将TSV中的相对路径转为绝对路径（基于TSV所在目录）
    # Auto-convert relative paths in TSV to absolute paths (based on TSV location)
    assemblies_tsv = _absolutize_tsv_paths(assemblies_tsv)

    # 构建参数列表|Build argument list
    args = ['swave.py', 'call']

    # 必需参数|Required parameters
    args.extend(['-i', assemblies_tsv])
    args.extend(['-r', ref_fasta])
    args.extend(['-g', gfa_file])
    args.extend(['-s', gfa_source])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if swave_path != '~/software/swave/Swave-main':
        args.extend(['--swave-path', swave_path])

    if decomposed_vcf:
        args.extend(['--decomposed-vcf', decomposed_vcf])

    if output_mode != 'auto':
        args.extend(['--output-mode', output_mode])

    if spec_samples:
        args.extend(['--spec-samples'] + list(spec_samples))

    if min_sv_size != 50:
        args.extend(['--min-sv-size', str(min_sv_size)])

    if max_sv_size != 1000000:
        args.extend(['--max-sv-size', str(max_sv_size)])

    if max_sv_comps != 5:
        args.extend(['--max-sv-comps', str(max_sv_comps)])

    # 布尔选项|Boolean options
    if dup_to_ins:
        args.append('--dup-to-ins')

    if remove_small:
        args.append('--remove-small')

    if force_reverse:
        args.append('--force-reverse')

    if threads != 12:
        args.extend(['-t', str(threads)])

    if minigraph_path != 'minigraph':
        args.extend(['--minigraph-path', minigraph_path])

    if gfatools_path != 'gfatools':
        args.extend(['--gfatools-path', gfatools_path])

    if spec_snarl:
        args.extend(['--spec-snarl', spec_snarl])

    if spec_path:
        args.extend(['--spec-path', spec_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        swave_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@swave.command(short_help='转换图路径为序列|Convert graph paths to sequences')
@click.option('--vcf-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('--gfa-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFA文件路径|GFA file path')
@click.option('--ref-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('--swave-path',
              default='~/software/swave/Swave-main',
              show_default=True,
              help='Swave软件路径|Swave software path')
@click.option('--output-path',
              help='输出路径|Output path')
@click.option('--force-pangenie',
              is_flag=True,
              help='强制输出满足pangenie要求的序列|Force output sequences to meet pangenie requirements')
def convert_seq(vcf_path, gfa_path, ref_path, swave_path, output_path, force_pangenie):
    """
    转换图路径为VCF的REF和ALT序列|Convert graph paths to sequence in VCF REF and ALT columns

    示例|Example: biopytools swave convert-seq --vcf-path input.vcf --gfa-path graph.gfa --ref-path ref.fa
    """
    swave_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['swave.py', 'convert_seq']

    args.extend(['--vcf-path', vcf_path])
    args.extend(['--gfa-path', gfa_path])
    args.extend(['--ref-path', ref_path])

    if swave_path != '~/software/swave/Swave-main':
        args.extend(['--swave-path', swave_path])

    if output_path:
        args.extend(['--output-path', output_path])

    if force_pangenie:
        args.append('--force-pangenie')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        swave_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@swave.command(short_help='转换VCF为GFA P lines|Convert VCF to GFA P lines')
@click.option('--gfa-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFA文件路径|GFA file path')
@click.option('--vcf-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('--ref-vcf-path',
              help='参考VCF文件路径|Reference VCF file path')
@click.option('--swave-path',
              default='~/software/swave/Swave-main',
              show_default=True,
              help='Swave软件路径|Swave software path')
@click.option('--output-path',
              help='输出路径|Output path')
@click.option('--force-vg',
              is_flag=True,
              help='强制输出满足vg要求的序列|Force output sequences to meet vg requirements')
def convert_plines(gfa_path, vcf_path, ref_vcf_path, swave_path, output_path, force_vg):
    """
    转换minigraph call VCF为GFA P lines|Convert minigraph call VCF into GFA P lines

    示例|Example: biopytools swave convert-plines --gfa-path graph.gfa --vcf-path input.vcf
    """
    swave_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['swave.py', 'convert_Plines']

    args.extend(['--gfa-path', gfa_path])
    args.extend(['--vcf-path', vcf_path])

    if ref_vcf_path:
        args.extend(['--ref-vcf-path', ref_vcf_path])

    if swave_path != '~/software/swave/Swave-main':
        args.extend(['--swave-path', swave_path])

    if output_path:
        args.extend(['--output-path', output_path])

    if force_vg:
        args.append('--force-vg')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        swave_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@swave.command(short_help='从VCF提取CSV|Extract CSV from VCF')
@click.option('--vcf-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('--swave-path',
              default='~/software/swave/Swave-main',
              show_default=True,
              help='Swave软件路径|Swave software path')
@click.option('--spec-csv',
              type=click.Choice(['INV', 'DUP', 'All']),
              help='特定CSV类型|Specific CSV type (INV/DUP/All)')
@click.option('--output-path',
              help='输出路径|Output path')
def extract_csv(vcf_path, swave_path, spec_csv, output_path):
    """
    从VCF提取CSV文件|Extract CSV from VCF

    示例|Example: biopytools swave extract-csv --vcf-path input.vcf --spec-csv INV
    """
    swave_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['swave.py', 'extract_csv']

    args.extend(['--vcf-path', vcf_path])

    if swave_path != '~/software/swave/Swave-main':
        args.extend(['--swave-path', swave_path])

    if spec_csv:
        args.extend(['--spec-csv', spec_csv])

    if output_path:
        args.extend(['--output-path', output_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        swave_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@swave.command(short_help='提取特定样本的SV|Extract SVs for specific samples')
@click.option('--vcf-path',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('--spec-samples',
              required=True,
              multiple=True,
              help='指定样本|Specify samples')
@click.option('--swave-path',
              default='~/software/swave/Swave-main',
              show_default=True,
              help='Swave软件路径|Swave software path')
@click.option('--output-path',
              help='输出路径|Output path')
def extract_sample(vcf_path, spec_samples, swave_path, output_path):
    """
    从VCF提取特定样本的SV|Extract SVs for specific samples from VCF

    示例|Example: biopytools swave extract-sample --vcf-path input.vcf --spec-samples sample1 sample2
    """
    swave_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['swave.py', 'extract_sample']

    args.extend(['--vcf-path', vcf_path])
    args.extend(['--spec-samples'] + list(spec_samples))

    if swave_path != '~/software/swave/Swave-main':
        args.extend(['--swave-path', swave_path])

    if output_path:
        args.extend(['--output-path', output_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        swave_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@swave.command(short_help='提取PAV矩阵|Extract PAV (Presence/Absence Variation) matrix')
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='swave converted VCF文件|swave converted VCF file')
@click.option('-o', '--output-file',
              default='pav_matrix.tsv',
              show_default=True,
              help='输出TSV文件|Output TSV file')
@click.option('--min-ac',
              type=int,
              default=1,
              show_default=True,
              help='最小等位基因数|Minimum allele count')
@click.option('--no-strip-prefix',
              is_flag=True,
              help='保留CHROM中的样本前缀|Keep sample prefix in CHROM')
@click.option('--svtype',
              multiple=True,
              help='仅保留指定SV类型|Only keep specified SV types (e.g., DUP INS DEL)')
def pav(vcf_file, output_file, min_ac, no_strip_prefix, svtype):
    """
    从swave converted VCF提取PAV矩阵|Extract PAV matrix from swave converted VCF

    示例|Example: biopytools swave pav -i converted.vcf -o pav_matrix.tsv
    """
    try:
        from ...swave.pav_extractor import PAVExtractor
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    extractor = PAVExtractor()
    try:
        result = extractor.extract(
            vcf_file=vcf_file,
            output_file=output_file,
            min_ac=min_ac,
            strip_prefix=not no_strip_prefix,
            svtype_only=list(svtype) if svtype else None
        )
        click.echo(f"PAV矩阵已保存|PAV matrix saved: {result}")
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
