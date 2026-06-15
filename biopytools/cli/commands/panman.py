"""
PanMAN泛基因组分析命令|PanMAN Pangenome Analysis Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载panman主函数|Lazy load panman main function"""
    try:
        from ...panman.main import main as panman_main
        return panman_main
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


@click.group()
@click.version_option(version="1.0.0")
def panman():
    """PanMAN泛基因组分析工具|PanMAN Pangenome Analysis Tool"""
    pass


@panman.command(
    short_help='构建PanMAN泛基因组网络|Build PanMAN pangenome network',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-p', '--pangraph',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='PanGraph JSON文件|PanGraph JSON file path')
@click.option('-g', '--gfa',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFA文件路径|GFA file path')
@click.option('-m', '--msa',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='MSA文件路径 (FASTA格式)|MSA file path (FASTA format)')
@click.option('-n', '--newick',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='Newick树文件路径|Newick tree file path')
@click.option('-o', '--output-prefix',
              default='output',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir',
              default='./panman_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--backend',
              type=click.Choice(['conda', 'docker', 'singularity']),
              default='conda',
              show_default=True,
              help='后端选择|Backend selection')
@click.option('--conda-env',
              default='panman_v.0.1.4',
              show_default=True,
              help='Conda环境名称|Conda environment name')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--sif-image',
              help='PanMAN Singularity SIF镜像路径|PanMAN Singularity SIF image path')
@click.option('--singularity-path',
              help='Singularity可执行文件路径|Singularity executable path')
def build(pangraph, gfa, msa, newick, output_prefix, output_dir, backend, conda_env, threads, sif_image, singularity_path):
    """
    构建PanMAN泛基因组网络|Build PanMAN pangenome network

    从PanGraph/GFA/MSA文件构建PanMAN网络|Build PanMAN network from PanGraph/GFA/MSA files

    示例|Example: biopytools panman build -P input.json -N tree.nwk -o my_panman
    """
    panman_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['panman.py', '--build']

    if pangraph:
        args.extend(['--pangraph', pangraph])
    if gfa:
        args.extend(['--gfa', gfa])
    if msa:
        args.extend(['--msa', msa])

    args.extend(['--newick', newick])

    if output_prefix != 'output':
        args.extend(['--output-prefix', output_prefix])

    if output_dir != './panman_output':
        args.extend(['--output-dir', output_dir])

    # 总是传递backend参数（因为默认值可能是singularity）|Always pass backend (default might be singularity)
    args.extend(['--backend', backend])

    if conda_env != 'panman_v.0.1.4':
        args.extend(['--conda-env', conda_env])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if sif_image:
        args.extend(['--sif-image', sif_image])

    if singularity_path:
        args.extend(['--singularity-path', singularity_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        panman_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@panman.command(
    short_help='从PanMAN提取数据|Extract data from PanMAN',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-panman',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='PanMAN文件路径|PanMAN file path')
@click.option('-o', '--output-prefix',
              default='output',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir',
              default='./panman_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-r', '--reference',
              help='参考序列名称 (VCF提取/重新扎根需要)|Reference sequence name (required for VCF/reroot)')
@click.option('--summary',
              is_flag=True,
              help='提取摘要统计|Extract summary statistics')
@click.option('--extract-fasta',
              is_flag=True,
              help='提取FASTA序列|Extract FASTA sequences')
@click.option('--extract-msa',
              is_flag=True,
              help='提取MSA比对|Extract MSA alignment')
@click.option('--vcf',
              is_flag=True,
              help='提取VCF变异|Extract VCF variants')
@click.option('--extract-gfa',
              is_flag=True,
              help='提取GFA格式|Extract GFA format')
@click.option('--extract-newick',
              is_flag=True,
              help='提取Newick树|Extract Newick tree')
@click.option('--extended-newick',
              is_flag=True,
              help='提取扩展Newick格式|Extract extended Newick format')
@click.option('--maf',
              is_flag=True,
              help='提取MAF格式|Extract MAF format')
@click.option('--aa',
              is_flag=True,
              help='提取氨基酸翻译|Extract amino acid translations')
@click.option('--subnet',
              is_flag=True,
              help='提取子网络|Extract subnet (requires --input-file)')
@click.option('--annotate',
              is_flag=True,
              help='注释节点|Annotate nodes (requires --input-file)')
@click.option('--reroot',
              is_flag=True,
              help='重新扎根树|Reroot tree (requires --reference)')
@click.option('--create-network',
              is_flag=True,
              help='创建网络|Create network (requires --input-file)')
@click.option('--print-mutations',
              is_flag=True,
              help='打印突变信息|Print mutations')
@click.option('--backend',
              type=click.Choice(['conda', 'docker', 'singularity']),
              default='conda',
              show_default=True,
              help='后端选择|Backend selection')
@click.option('--conda-env',
              default='panman_v.0.1.4',
              show_default=True,
              help='Conda环境名称|Conda environment name')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--sif-image',
              help='PanMAN Singularity SIF镜像路径|PanMAN Singularity SIF image path')
@click.option('--singularity-path',
              help='Singularity可执行文件路径|Singularity executable path')
@click.option('--acr',
              type=click.Choice(['fitch', 'mppa']),
              default='fitch',
              show_default=True,
              help='ACR方法|ACR method (fitch/mppa)')
@click.option('--tree-id',
              help='树ID (VCF提取可选)|Tree ID (optional for VCF extraction)')
@click.option('--input-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入文件路径 (subnet/annotate/create-network需要)|Input file path (required for subnet/annotate/create-network)')
@click.option('--range-index',
              help='范围查询index参数|Range query index parameter')
@click.option('--range-start',
              type=int,
              help='范围查询起始坐标|Range query start coordinate')
@click.option('--range-end',
              type=int,
              help='范围查询结束坐标|Range query end coordinate')
def extract(input_panman, output_prefix, output_dir, reference, summary, extract_fasta, extract_msa, vcf,
            extract_gfa, extract_newick, extended_newick, maf, aa, subnet, annotate, reroot, create_network,
            print_mutations, backend, conda_env, threads, sif_image, singularity_path, acr, tree_id, input_file, range_index,
            range_start, range_end):
    """
    从PanMAN提取数据|Extract data from PanMAN

    支持提取多种格式的数据|Support extracting data in multiple formats

    示例|Example: biopytools panman extract -I data.panman --summary --vcf --extract-fasta
    """
    panman_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['panman.py', '--extract']

    args.extend(['--input-panman', input_panman])

    if output_prefix != 'output':
        args.extend(['--output-prefix', output_prefix])

    if output_dir != './panman_output':
        args.extend(['--output-dir', output_dir])

    if reference:
        args.extend(['--reference', reference])

    if summary:
        args.append('--summary')
    if extract_fasta:
        args.append('--extract-fasta')
    if extract_msa:
        args.append('--extract-msa')
    if vcf:
        args.append('--vcf')
    if extract_gfa:
        args.append('--extract-gfa')
    if extract_newick:
        args.append('--extract-newick')
    if extended_newick:
        args.append('--extended-newick')
    if maf:
        args.append('--maf')
    if aa:
        args.append('--aa')
    if subnet:
        args.append('--subnet')
    if annotate:
        args.append('--annotate')
    if reroot:
        args.append('--reroot')
    if create_network:
        args.append('--create-network')
    if print_mutations:
        args.append('--print-mutations')

    # 总是传递backend参数（因为默认值可能是singularity）|Always pass backend (default might be singularity)
    args.extend(['--backend', backend])

    if conda_env != 'panman_v.0.1.4':
        args.extend(['--conda-env', conda_env])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if sif_image:
        args.extend(['--sif-image', sif_image])

    if singularity_path:
        args.extend(['--singularity-path', singularity_path])

    if acr != 'fitch':
        args.extend(['--acr', acr])

    if tree_id:
        args.extend(['--tree-id', tree_id])

    if input_file:
        args.extend(['--input-file', input_file])

    if range_index:
        args.extend(['--range-index', range_index])
    if range_start is not None:
        args.extend(['--range-start', str(range_start)])
    if range_end is not None:
        args.extend(['--range-end', str(range_end)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        panman_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@panman.command(
    short_help='从FASTA生成PanGraph|Generate PanGraph from FASTA',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--fasta',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件路径|Input FASTA file path')
@click.option('-o', '--output-prefix',
              default='output',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir',
              default='./panman_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--pangraph-path',
              help='PanGraph可执行文件路径|PanGraph executable path')
@click.option('--pangraph-sif',
              help='PanGraph Singularity SIF镜像路径|PanGraph Singularity SIF image path')
@click.option('--singularity-path',
              help='Singularity可执行文件路径|Singularity executable path')
@click.option('--backend',
              type=click.Choice(['conda', 'docker', 'singularity']),
              default='singularity',
              show_default=True,
              help='后端选择|Backend selection')
def generate_pangraph(fasta, output_prefix, output_dir, threads, pangraph_path, pangraph_sif, singularity_path, backend):
    """
    从FASTA文件生成PanGraph JSON|Generate PanGraph JSON from FASTA file

    从FASTA序列生成PanGraph JSON和Newick树文件|Generate PanGraph JSON and Newick tree from FASTA sequences

    示例|Example: biopytools panman generate-pangraph -i input.fa -o my_pangraph
    """
    panman_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['panman.py', '--generate-pangraph']

    args.extend(['--fasta', fasta])

    if output_prefix != 'output':
        args.extend(['--output-prefix', output_prefix])

    if output_dir != './panman_output':
        args.extend(['--output-dir', output_dir])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if pangraph_path:
        args.extend(['--pangraph-path', pangraph_path])

    if pangraph_sif:
        args.extend(['--pangraph-sif', pangraph_sif])

    if singularity_path:
        args.extend(['--singularity-path', singularity_path])

    if backend != 'singularity':
        args.extend(['--backend', backend])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        panman_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
