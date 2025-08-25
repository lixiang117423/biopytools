"""
🧬 ADMIXTURE群体结构分析命令 | ADMIXTURE Population Structure Analysis Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_admixture_main():
    """懒加载admixture main函数 | Lazy load admixture main function"""
    try:
        from ...admixture.main import main as admixture_main
        return admixture_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="ADMIXTURE群体结构分析",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 输入VCF文件路径 | Input VCF file path')
@click.option('--output', '-o',
              default='admixture_results',
              type=click.Path(),
              help='📁 输出目录 | Output directory')
@click.option('--min-k', '-k',
              type=int,
              default=2,
              help='🔢 最小K值 (默认: 2) | Minimum K value (default: 2)')
@click.option('--max-k', '-K',
              type=int,
              default=10,
              help='🔢 最大K值 (默认: 10) | Maximum K value (default: 10)')
@click.option('--cv-folds', '-c',
              type=int,
              default=5,
              help='🔄 交叉验证折数 (默认: 5) | Cross-validation folds (default: 5)')
@click.option('--threads', '-t',
              type=int,
              default=4,
              help='⚡ 线程数 (默认: 4) | Number of threads (default: 4)')
@click.option('--maf', '-m',
              type=float,
              default=0.01,
              help='📊 MAF阈值 (默认: 0.01) | MAF threshold (default: 0.01)')
@click.option('--missing', '-M',
              type=float,
              default=0.1,
              help='❓ 缺失率阈值 (默认: 0.1) | Missing rate threshold (default: 0.1)')
@click.option('--hwe', '-H',
              type=float,
              default=1e-6,
              help='⚖️ HWE p值阈值 (默认: 1e-6) | HWE p-value threshold (default: 1e-6)')
@click.option('--skip-preprocessing', '-s',
              is_flag=True,
              help='⏩ 跳过VCF预处理 | Skip VCF preprocessing')
@click.option('--keep-intermediate', '-i',
              is_flag=True,
              help='💾 保留中间文件 | Keep intermediate files')
def admixture(vcf, output, min_k, max_k, cv_folds, threads, maf, missing,
              hwe, skip_preprocessing, keep_intermediate):
    """
    ADMIXTURE群体结构分析.

    示例 | Examples:
    
    \b
    # 🎯 基本分析
    biopytools admixture -v input.vcf -o results
    
    \b
    # 🔧 指定K值范围和线程数
    biopytools admixture -v input.vcf -o results -k 3 -K 8 -t 8
    
    \b
    # 📊 自定义质控参数
    biopytools admixture -v input.vcf -o results \\
        --maf 0.05 --missing 0.05 --hwe 1e-5
    
    \b
    # ⏩ 跳过预处理并保留中间文件
    biopytools admixture -v clean.vcf -o results \\
        --skip-preprocessing --keep-intermediate
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    admixture_main = _lazy_import_admixture_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'admixture']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-v', vcf])
    
    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    args.extend(['-o', output])  # 总是添加输出目录 | Always add output directory
    
    if min_k != 2:
        args.extend(['-k', str(min_k)])
    
    if max_k != 10:
        args.extend(['-K', str(max_k)])
    
    if cv_folds != 5:
        args.extend(['-c', str(cv_folds)])
    
    if threads != 4:
        args.extend(['-t', str(threads)])
    
    # 质控参数 📊 | Quality control parameters
    if maf != 0.01:
        args.extend(['-m', str(maf)])
    
    if missing != 0.1:
        args.extend(['-M', str(missing)])
    
    if hwe != 1e-6:
        args.extend(['-H', str(hwe)])
    
    # 布尔选项 🚩 | Boolean options
    if skip_preprocessing:
        args.append('--skip-preprocessing')
    
    if keep_intermediate:
        args.append('--keep-intermediate')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        admixture_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


# 🔥 额外优化：为其他命令提供通用的懒加载装饰器 | Additional optimization: generic lazy loading decorator
def lazy_command(module_path, main_func_name='main'):
    """
    通用懒加载装饰器 | Generic lazy loading decorator
    
    Args:
        module_path: 模块路径，如 'biopytools.admixture.main'
        main_func_name: main函数名称，默认为'main'
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # 懒加载模块
            try:
                module = __import__(module_path, fromlist=[main_func_name])
                main_func = getattr(module, main_func_name)
            except (ImportError, AttributeError) as e:
                click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
                sys.exit(1)
            
            # 调用原始函数，传入main_func
            return func(main_func, *args, **kwargs)
        return wrapper
    return decoratora