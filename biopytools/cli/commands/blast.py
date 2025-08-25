"""
🧬 BLAST序列比对分析命令 | BLAST Sequence Alignment Analysis Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_blast_main():
    """懒加载blast main函数 | Lazy load blast main function"""
    try:
        from ...blast_analysis.main import main as blast_main
        return blast_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


def _validate_required_file_exists(file_path):
    """验证必需文件是否存在（仅在非帮助模式下）| Validate required file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="BLAST序列比对分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📁 输入文件或目录路径 | Input file or directory path')
@click.option('--sample-map-file', '-s',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧪 样品映射文件，格式：文件路径<TAB>样品名称 | Sample mapping file, format: file_path<TAB>sample_name')
@click.option('--target-file', '-t',
              required=True,
              callback=lambda ctx, param, value: _validate_required_file_exists(value) if value else None,
              help='🎯 目标基因序列文件 | Target gene sequence file')
@click.option('--output-dir', '-o',
              default='./blast_output',
              type=click.Path(),
              help='📁 输出目录 (默认: ./blast_output) | Output directory (default: ./blast_output)')
@click.option('--blast-type',
              type=click.Choice(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']),
              default='blastn',
              help='🔬 BLAST程序类型 (默认: blastn) | BLAST program type (default: blastn)')
@click.option('--evalue', '-e',
              type=float,
              default=1e-5,
              help='🔢 E-value阈值 (默认: 1e-5) | E-value threshold (default: 1e-5)')
@click.option('--max-target-seqs',
              type=int,
              default=10,
              help='🎯 最大目标序列数 (默认: 10) | Maximum target sequences (default: 10)')
@click.option('--word-size',
              type=int,
              default=11,
              help='📏 词大小 (默认: 11) | Word size (default: 11)')
@click.option('--threads', '-j',
              type=int,
              default=88,
              help='⚡ 线程数 (默认: 88) | Number of threads (default: 88)')
@click.option('--input-suffix',
              default='*.fa',
              help='📁 输入文件后缀模式 (默认: *.fa) | Input file suffix pattern (default: *.fa)')
@click.option('--target-db-type',
              type=click.Choice(['nucl', 'prot']),
              default='nucl',
              help='🗄️ 目标数据库类型 (默认: nucl) | Target database type (default: nucl)')
@click.option('--min-identity',
              type=float,
              default=70.0,
              help='📊 最小序列相似度 %% (默认: 70.0) | Minimum sequence identity %% (default: 70.0)')
@click.option('--min-coverage',
              type=float,
              default=50.0,
              help='📐 最小覆盖度 %% (默认: 50.0) | Minimum coverage %% (default: 50.0)')
@click.option('--high-quality-evalue',
              type=float,
              default=1e-10,
              help='🌟 高质量比对E-value阈值 (默认: 1e-10) | High quality E-value threshold (default: 1e-10)')
@click.option('--auto-detect-samples',
              is_flag=True,
              default=True,
              help='🔍 自动检测样品名称 (默认: True) | Auto detect sample names (default: True)')
@click.option('--sample-name-pattern',
              default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
              help='🔍 样品名称提取正则表达式 | Sample name extraction regex pattern')
@click.option('--makeblastdb-path',
              default='makeblastdb',
              help='🗄️ makeblastdb程序路径 (默认: makeblastdb) | makeblastdb program path (default: makeblastdb)')
@click.option('--blastn-path',
              default='blastn',
              help='🔬 blastn程序路径 (默认: blastn) | blastn program path (default: blastn)')
@click.option('--blastp-path',
              default='blastp',
              help='🔬 blastp程序路径 (默认: blastp) | blastp program path (default: blastp)')
@click.option('--blastx-path',
              default='blastx',
              help='🔬 blastx程序路径 (默认: blastx) | blastx program path (default: blastx)')
@click.option('--tblastn-path',
              default='tblastn',
              help='🔬 tblastn程序路径 (默认: tblastn) | tblastn program path (default: tblastn)')
@click.option('--tblastx-path',
              default='tblastx',
              help='🔬 tblastx程序路径 (默认: tblastx) | tblastx program path (default: tblastx)')
def blast(input, sample_map_file, target_file, output_dir, blast_type, evalue, 
         max_target_seqs, word_size, threads, input_suffix, target_db_type,
         min_identity, min_coverage, high_quality_evalue, auto_detect_samples,
         sample_name_pattern, makeblastdb_path, blastn_path, blastp_path,
         blastx_path, tblastn_path, tblastx_path):
    """
    BLAST序列比对分析工具.
    
    对DNA/蛋白质序列进行BLAST比对分析，支持多种BLAST程序类型，
    自动处理样品映射，生成详细的比对结果和统计报告。
    
    示例 | Examples:
    
    \b
    # 🎯 基本DNA序列比对
    biopytools blast -i sequences/ -t target_genes.fa -o blast_results
    
    \b
    # 🔬 蛋白质序列比对
    biopytools blast -i proteins/ -t targets.fa -o results \\
        --blast-type blastp --target-db-type prot -j 32
    
    \b
    # 📊 高质量比对筛选
    biopytools blast -i sequences/ -t nlr.fa -o results \\
        --min-identity 80 --min-coverage 70 --high-quality-evalue 1e-15
    
    \b
    # 🧪 使用现有样品映射文件
    biopytools blast -s sample_map.txt -t targets.fa -o results
    
    \b
    # 🔧 自定义BLAST参数
    biopytools blast -i sequences/ -t targets.fa -o results \\
        --evalue 1e-10 --max-target-seqs 20 --word-size 15 \\
        --blast-type blastx -j 64
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    blast_main = _lazy_import_blast_main()
    
    # 验证输入参数 📋 | Validate input parameters
    if not input and not sample_map_file:
        click.echo("❌ 错误：必须指定输入路径(-i)或样品映射文件(-s)中的一个", err=True)
        sys.exit(1)
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'blast']
    
    # 输入数据参数 📁 | Input data parameters
    if input:
        args.extend(['-i', input])
    
    if sample_map_file:
        args.extend(['-s', sample_map_file])
    
    # 必需参数 📋 | Required parameters
    args.extend(['-t', target_file])
    
    # 基本参数 ⚙️ | Basic parameters
    if output_dir != './blast_output':
        args.extend(['-o', output_dir])
    
    if blast_type != 'blastn':
        args.extend(['--blast-type', blast_type])
    
    # BLAST参数 🔬 | BLAST parameters
    if evalue != 1e-5:
        args.extend(['-e', str(evalue)])
    
    if max_target_seqs != 10:
        args.extend(['--max-target-seqs', str(max_target_seqs)])
    
    if word_size != 11:
        args.extend(['--word-size', str(word_size)])
    
    if threads != 88:
        args.extend(['--threads', str(threads)])
    
    # 文件格式参数 📄 | File format parameters
    if input_suffix != '*.fa':
        args.extend(['--input-suffix', input_suffix])
    
    if target_db_type != 'nucl':
        args.extend(['--target-db-type', target_db_type])
    
    # 过滤参数 📊 | Filtering parameters
    if min_identity != 70.0:
        args.extend(['--min-identity', str(min_identity)])
    
    if min_coverage != 50.0:
        args.extend(['--min-coverage', str(min_coverage)])
    
    if high_quality_evalue != 1e-10:
        args.extend(['--high-quality-evalue', str(high_quality_evalue)])
    
    # 样品检测参数 🔍 | Sample detection parameters
    if not auto_detect_samples:  # 只在False时添加，因为默认是True
        # 注意：原始脚本使用store_true，所以我们不需要显式添加True的情况
        pass
    
    if sample_name_pattern != r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$':
        args.extend(['--sample-name-pattern', sample_name_pattern])
    
    # 工具路径参数 🛠️ | Tool path parameters
    if makeblastdb_path != 'makeblastdb':
        args.extend(['--makeblastdb-path', makeblastdb_path])
    
    if blastn_path != 'blastn':
        args.extend(['--blastn-path', blastn_path])
    
    if blastp_path != 'blastp':
        args.extend(['--blastp-path', blastp_path])
    
    if blastx_path != 'blastx':
        args.extend(['--blastx-path', blastx_path])
    
    if tblastn_path != 'tblastn':
        args.extend(['--tblastn-path', tblastn_path])
    
    if tblastx_path != 'tblastx':
        args.extend(['--tblastx-path', tblastx_path])
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        blast_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv