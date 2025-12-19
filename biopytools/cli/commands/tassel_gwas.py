"""
🌾 TASSEL GWAS分析命令 | TASSEL GWAS Analysis Command
🚀 高级优化版本：解决--help响应速度问题 ⚡️
"""

import click
import sys
import os


def _lazy_import_tassel_gwas_main():
    """😴 懒加载tassel_gwas main函数 | Lazy load tassel_gwas main function"""
    try:
        from ...tassel_gwas.main import main as tassel_gwas_main
        return tassel_gwas_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """❓ 检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """🔍 验证文件是否存在（仅在非帮助模式下）| Validate file exists (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"📂❌ 文件不存在 | File does not exist: {file_path}")
    return file_path


def _validate_directory_exists(dir_path):
    """🔍 验证目录是否存在（仅在非帮助模式下）| Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(dir_path):
        try:
            os.makedirs(dir_path, exist_ok=True)
        except OSError:
            raise click.BadParameter(f"📂❌ 无法创建输出目录 | Cannot create output directory: {dir_path}")
    return dir_path


@click.command(short_help="🌾 TASSEL GWAS分析工具")
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📊 VCF基因型文件路径 | VCF genotype file path')
@click.option('--pheno', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📋 表型文件路径 | Phenotype file path (可包含多个表型)')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='📁 输出目录路径 | Output directory path')
@click.option('--model', '-m',
              type=click.Choice(['GLM', 'MLM', 'BOTH']),
              default='MLM',
              help='🧬 GWAS模型 | GWAS model (GLM/MLM/BOTH) (default: MLM)')
@click.option('--memory',
              default='100g',
              help='💾 Java最大内存 | Java maximum memory (default: 100g)')
@click.option('--threads', '-t',
              type=int,
              default=4,
              help='🔢 并行线程数 | Number of parallel threads (default: 4)')
@click.option('--maf',
              type=float,
              help='🎯 最小等位基因频率过滤 | Minimum allele frequency filter')
@click.option('--miss',
              type=float,
              help='❓ 最大缺失率过滤 | Maximum missing rate filter')
@click.option('--pca-components',
              type=int,
              default=5,
              help='🧮 PCA主成分数量（用作MLM协变量）| Number of PCA components for MLM covariates (default: 5)')
@click.option('--q-matrix', '-q',
              type=click.Path(exists=True),
              help='🏛️ 群体结构Q矩阵文件 | Population structure Q matrix file')
@click.option('--kinship', '-k',
              type=click.Path(exists=True),
              help='🧬 亲缘关系K矩阵文件 | Kinship matrix file')
@click.option('--tassel-path',
              type=click.Path(exists=True),
              help='🔧 TASSEL安装路径 | TASSEL installation path')
@click.option('--skip-sort',
              is_flag=True,
              help='⏭️ 跳过VCF排序 | Skip VCF sorting')
@click.option('--keep-temp',
              is_flag=True,
              help='💾 保留临时文件 | Keep temporary files')
@click.option('--parallel',
              is_flag=True,
              help='🚀 并行处理多个表型 | Parallel process multiple traits')
@click.option('--workers',
              type=int,
              default=4,
              help='👥 并行工作线程数 | Number of parallel workers (default: 4)')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARN', 'ERROR']),
              default='INFO',
              help='📝 日志级别 | Log level (default: INFO)')
def tassel_gwas(vcf, pheno, output, model, memory, threads, maf, miss, pca_components, q_matrix,
                kinship, tassel_path, skip_sort, keep_temp, parallel, workers, log_level):
    """
    🌾 TASSEL GWAS分析工具 | TASSEL GWAS Analysis Tool

    基于TASSEL软件的全基因组关联分析工具，支持自动识别表型数量并批量处理。

    💡 功能特性 | Features:

    • 📋 自动识别表型数量 | Automatic detection of trait count
    • 🚀 支持并行处理 | Support parallel processing
    • 🧬 支持GLM/MLM/BOTH模型 | Support GLM/MLM/BOTH models
    • 📊 自动计算Kinship矩阵 | Automatic Kinship matrix calculation
    • 🧮 可配置PCA主成分数量 | Configurable number of PCA components (default: 5)
    • 🏛️ 支持群体结构校正 | Support population structure correction
    • ⚙️ 灵活的质控参数 | Flexible quality control parameters
    • 📈 生成曼哈顿图输入文件 | Generate Manhattan plot input files
    • 📊 生成处理报告 | Generate processing report

    📚 示例 | Examples:

    \b
    # 🎯 基本分析 - 自动识别并处理所有表型
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results

    \b
    # ⚙️ 指定模型和参数
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \\
        --model MLM --memory 200g --maf 0.05 --miss 0.1

    \b
    # 🧮 指定PCA主成分数量（默认5个）
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \\
        --pca-components 10

    \b
    # 🏛️ 使用Q矩阵校正群体结构
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \\
        --q-matrix Q.txt

    \b
    # 🚀 并行处理多个表型
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \\
        --parallel --workers 8

    \b
    # 🧬 运行GLM和MLM两种模型
    biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \\
        --model BOTH --keep-temp
    """

    # 😴 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    tassel_gwas_main = _lazy_import_tassel_gwas_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['tassel_gwas.py']

    # 必需参数 📋 | Required parameters
    args.extend(['--vcf', vcf])
    args.extend(['--pheno', pheno])
    args.extend(['--output', output])

    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    args.extend(['--model', model])

    if memory != '100g':
        args.extend(['--memory', memory])

    if threads != 4:
        args.extend(['--threads', str(threads)])

    if maf is not None:
        args.extend(['--maf', str(maf)])

    if miss is not None:
        args.extend(['--miss', str(miss)])

    if pca_components != 5:
        args.extend(['--pca-components', str(pca_components)])

    if q_matrix:
        args.extend(['--q-matrix', q_matrix])

    if kinship:
        args.extend(['--kinship', kinship])

    if tassel_path:
        args.extend(['--tassel-path', tassel_path])

    if skip_sort:
        args.append('--skip-sort')

    if keep_temp:
        args.append('--keep-temp')

    if parallel:
        args.append('--parallel')

    if workers != 4:
        args.extend(['--workers', str(workers)])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        tassel_gwas_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 运行错误 | Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv