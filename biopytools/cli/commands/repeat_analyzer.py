"""
🧬 RepeatAnalyzerCLI命令 | RepeatAnalyzer CLI Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_repeat_main():
    """懒加载repeat analyzer main函数 | Lazy load repeat analyzer main function"""
    try:
        from ...repeat_analyzer.main import main as repeat_main
        return repeat_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_file(file_path):
    """验证输入基因组文件是否存在（仅在非帮助模式下）| Validate input genome file existence (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(file_path):
            raise click.BadParameter(f"输入基因组文件不存在 | Input genome file does not exist: {file_path}")
        if not os.path.isfile(file_path):
            raise click.BadParameter(f"输入路径不是文件 | Input path is not a file: {file_path}")
    return file_path


def _validate_output_dir(dir_path):
    """验证输出目录路径（仅在非帮助模式下）| Validate output directory path (only in non-help mode)"""
    if not _is_help_request():
        # 输出目录可以不存在，程序会创建
        parent_dir = os.path.dirname(os.path.abspath(dir_path))
        if parent_dir and not os.path.exists(parent_dir):
            raise click.BadParameter(f"输出目录的父目录不存在 | Parent directory of output does not exist: {parent_dir}")
    return dir_path


@click.command(short_help='植物基因组重复序列分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_file(value) if value else None,
              help='🧬 输入基因组FASTA文件路径 | Input genome FASTA file path')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_output_dir(value) if value else None,
              help='📂 输出目录 | Output directory')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 (默认: 88) | Number of threads (default: 88)')
@click.option('--skip-modeler',
              is_flag=True,
              help='⏭️ 跳过RepeatModeler步骤 | Skip RepeatModeler step')
@click.option('--skip-ltr',
              is_flag=True,
              help='⏭️ 跳过LTR分析步骤 | Skip LTR analysis step')
@click.option('--repeatmodeler-path',
              default='RepeatModeler',
              help='🔧 RepeatModeler程序路径 (默认: RepeatModeler) | RepeatModeler program path (default: RepeatModeler)')
@click.option('--ltr-finder-path',
              default='ltr_finder',
              help='🔧 LTR_FINDER程序路径 (默认: ltr_finder) | LTR_FINDER program path (default: ltr_finder)')
@click.option('--ltrharvest-path',
              default='gt ltrharvest',
              help='🔧 LTRharvest程序路径 (默认: gt ltrharvest) | LTRharvest program path (default: gt ltrharvest)')
@click.option('--ltr-retriever-path',
              default='LTR_retriever',
              help='🔧 LTR_retriever程序路径 (默认: LTR_retriever) | LTR_retriever program path (default: LTR_retriever)')
@click.option('--repeatmasker-path',
              default='RepeatMasker',
              help='🔧 RepeatMasker程序路径 (默认: RepeatMasker) | RepeatMasker program path (default: RepeatMasker)')
@click.option('--tesorter-path',
              default='TEsorter',
              help='🔧 TEsorter程序路径 (默认: TEsorter) | TEsorter program path (default: TEsorter)')
def repeat_analyzer(input, output, threads, skip_modeler, skip_ltr,
                   repeatmodeler_path, ltr_finder_path, ltrharvest_path,
                   ltr_retriever_path, repeatmasker_path, tesorter_path):
    """
    🧬 植物基因组重复序列和转座元件分析工具
    
    使用RepeatModeler、LTR工具、RepeatMasker和TEsorter对植物基因组进行
    全面的重复序列识别和分类，支持EDTA下游分析格式。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本重复序列分析
    biopytools repeat-analyzer -i genome.fasta -o repeat_results
    
    \b
    # ⚡ 大基因组高线程分析
    biopytools repeat-analyzer -i large_genome.fa -o results -t 64
    
    \b
    # ⏩ 跳过RepeatModeler加速分析
    biopytools repeat-analyzer -i genome.fasta -o results --skip-modeler
    
    \b
    # 🔧 只进行RepeatMasker和TEsorter分析
    biopytools repeat-analyzer -i genome.fa -o results \\
        --skip-modeler --skip-ltr
    
    \b
    # 🛠️ 自定义工具路径
    biopytools repeat-analyzer -i genome.fasta -o results \\
        --repeatmasker-path ~/.local/bin/RepeatMasker \\
        --tesorter-path ~/.local/bin/TEsorter -t 32
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    repeat_main = _lazy_import_repeat_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['repeat_analyzer.py']  # 模拟脚本名，避免传递biopytools子命令
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）⚙️ | Optional parameters (add only when non-default)
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # 布尔选项 🚩 | Boolean options
    if skip_modeler:
        args.append('--skip-modeler')
    
    if skip_ltr:
        args.append('--skip-ltr')
    
    # 工具路径参数（只在非默认值时添加）🔧 | Tool path parameters (add only when non-default)
    if repeatmodeler_path != 'RepeatModeler':
        args.extend(['--repeatmodeler-path', repeatmodeler_path])
    
    if ltr_finder_path != 'ltr_finder':
        args.extend(['--ltr-finder-path', ltr_finder_path])
        
    if ltrharvest_path != 'gt ltrharvest':
        args.extend(['--ltrharvest-path', ltrharvest_path])
        
    if ltr_retriever_path != 'LTR_retriever':
        args.extend(['--ltr-retriever-path', ltr_retriever_path])
    
    if repeatmasker_path != 'RepeatMasker':
        args.extend(['--repeatmasker-path', repeatmasker_path])
    
    if tesorter_path != 'TEsorter':
        args.extend(['--tesorter-path', tesorter_path])
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        repeat_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv