"""
EGAPx批量运行配置生成工具 | EGAPx Batch Config Generator
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_egapx_batch_main():
    """懒加载egapx-batch main函数 | Lazy load egapx-batch main function"""
    try:
        from ...egapx_batch.main import main as egapx_batch_main
        return egapx_batch_main
    except ImportError as e:
        click.echo(f"[ERROR] 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_genome(file_path):
    """验证基因组文件是否存在（仅在非帮助模式下）| Validate genome file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("基因组文件路径不能为空 | Genome file path cannot be empty")
    if not os.path.exists(file_path):
        raise click.BadParameter(f"基因组文件不存在 | Genome file not found: {file_path}")
    if not file_path.endswith(('.fa', '.fa.gz', '.fasta', '.fasta.gz')):
        raise click.BadParameter(f"基因组文件格式不正确 | Invalid genome file format: {file_path}")
    return file_path


def _validate_yaml(file_path):
    """验证YAML模板文件是否存在（仅在非帮助模式下）| Validate YAML template (only in non-help mode)"""
    if _is_help_request():
        return file_path
    # 如果为None，表示使用内置默认模板 | If None, use built-in default template
    if file_path is None:
        return None
    if not os.path.exists(file_path):
        raise click.BadParameter(f"YAML模板不存在 | YAML template not found: {file_path}")
    if not file_path.endswith('.yaml'):
        raise click.BadParameter(f"YAML模板格式不正确 | Invalid YAML template format: {file_path}")
    return file_path


def _validate_script(file_path):
    """验证脚本模板文件是否存在（仅在非帮助模式下）| Validate script template (only in non-help mode)"""
    if _is_help_request():
        return file_path
    # 如果为None，表示使用内置默认模板 | If None, use built-in default template
    if file_path is None:
        return None
    if not os.path.exists(file_path):
        raise click.BadParameter(f"脚本模板不存在 | Script template not found: {file_path}")
    if not file_path.endswith('.sh'):
        raise click.BadParameter(f"脚本模板格式不正确 | Invalid script template format: {file_path}")
    return file_path


@click.command(short_help="EGAPx批量运行配置生成工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_genome(value) if value else None,
              help='[FILE] 基因组FASTA文件路径 | Genome FASTA file path')
@click.option('--yaml', '-y',
              required=False,
              default=None,
              callback=lambda ctx, param, value: _validate_yaml(value) if value else None,
              help='[FILE] YAML模板文件路径 (可选，内置默认模板) | YAML template file path (optional, built-in default template)')
@click.option('--script', '-s',
              required=False,
              default=None,
              callback=lambda ctx, param, value: _validate_script(value) if value else None,
              help='[FILE] Shell脚本模板路径 (可选，内置默认模板) | Shell script template path (optional, built-in default template)')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='[DIR] 输出目录路径 | Output directory path')
@click.option('--egapx', '-e',
              default='/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx',
              type=click.Path(),
              help='[PATH] EGAPx安装路径 (默认: /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx) | EGAPx installation path')
@click.option('--chr-prefix', '-p',
              help='[STR] 染色体前缀过滤 | Chromosome prefix filter')
@click.option('--locus-prefix',
              default='',
              help='[STR] locus标签前缀 (默认: 空) | Locus tag prefix (default: empty)')
@click.option('--report-name',
              default='EGAPx',
              help='[STR] 报告名称 (默认: EGAPx) | Report name (default: EGAPx)')
@click.option('--short-reads',
              default='',
              help='[FILE] 短读长测序数据文件路径 | Short reads file path')
@click.option('--long-reads',
              default='',
              help='[FILE] 长读长测序数据文件路径 | Long reads file path')
def egapx_batch(genome, yaml, script, output, egapx, chr_prefix, locus_prefix, report_name, short_reads, long_reads):
    """
    EGAPx批量运行配置生成工具

    按染色体拆分基因组，为每个染色体生成独立的EGAPx运行配置文件
    和执行脚本，支持并行化批量处理。

    核心功能 | Core Features:
    - 按染色体/scaffold拆分基因组
    - 为每个序列生成独立的YAML配置
    - 为每个序列生成运行脚本
    - 自动创建EGAPx软链接
    - 生成任务列表和并行执行脚本
    - 支持自定义locus标签前缀
    - 支持自定义报告名称
    - 支持染色体前缀过滤

    使用场景 | Use Cases:
    - 大型基因组需要分染色体注释
    - 并行化EGAPx提高处理速度
    - 批量基因预测和注释
    - 多样本基因组比较分析

    示例 | Examples:

    \b
    # 基本用法 (使用内置模板)
    biopytools egapx-batch \\
        -g genome.fa \\
        -o output_dir

    \b
    # 自定义locus标签前缀和报告名
    biopytools egapx-batch \\
        -g genome.fa \\
        -o output_dir \\
        --locus-prefix Gene \\
        --report-name MyEGAPx

    \b
    # 使用自定义模板
    biopytools egapx-batch \\
        -g genome.fa \\
        -o output_dir \\
        -y template.yaml \\
        -s template.sh

    \b
    # 只处理主染色体（Chr开头）
    biopytools egapx-batch \\
        -g genome.fa \\
        -o output_dir \\
        --chr-prefix Chr

    \b
    # 自定义EGAPx路径
    biopytools egapx-batch \\
        -g genome.fa \\
        -o output_dir \\
        --egapx /custom/path/to/egapx

    输出说明 | Output Description:

    为每个染色体/序列生成独立的目录结构：
    ```
    output_dir/
    ├── Chr01/
    │   ├── Chr01.fa          # 染色体序列
    │   ├── Chr01.yaml        # YAML配置
    │   ├── egapx_Chr01.sh    # 运行脚本
    │   ├── work/             # 工作目录
    │   ├── output/           # 输出目录
    │   └── [EGAPx软链接]     # EGAPx symlinks
    ├── Chr02/
    │   └── ...
    ├── all_jobs_submit.list.sh    # 任务列表
    └── run_all_parallel.sh         # 并行执行脚本
    ```

    执行方式 | Execution Methods:

    \b
    # 1. 顺序执行
    bash output_dir/all_jobs_submit.list.sh

    \b
    # 2. 并行执行（推荐）
    bash output_dir/run_all_parallel.sh 4

    \b
    # 3. 使用GNU parallel
    cat output_dir/all_jobs_submit.list.sh | parallel -j 4

    模板文件要求 | Template File Requirements:

    工具内置了默认模板，无需额外指定。如需自定义，可以提供以下模板：

    **YAML模板** 需要包含以下字段：
    ```yaml
    genome: /path/to/genome.fa
    taxid: 71234
    short_reads: /path/to/short_reads.txt
    long_reads: /path/to/long_reads.txt
    locus_tag_prefix: Target
    ```

    **Shell脚本模板** 需要包含EGAPx运行命令：
    ```bash
    python3 ui/egapx.py config.yaml \\
        -e singularity \\
        -w work_dir \\
        -o output_dir \\
        -r report_name
    ```

    如果不提供-y和-s参数，将使用内置默认模板。

    故障排除 | Troubleshooting:

    1. "基因组文件不存在"
       - 检查-g参数指定的文件路径
       - 确认文件格式为.fasta/.fa

    2. "未拆分出任何序列"
       - 确认基因组文件格式正确
       - 检查是否使用了错误的--chr-prefix
       - 尝试不使用--chr-prefix参数

    3. "YAML模板不存在"
       - 如果使用-y参数，确认文件路径正确
       - 或者不使用-y参数，使用内置默认模板

    4. "脚本模板不存在"
       - 如果使用-s参数，确认文件路径正确
       - 或者不使用-s参数，使用内置默认模板

    5. "EGAPx路径不存在"
       - 使用--egapx参数指定正确的EGAPx安装路径
       - 或使用默认路径: /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx

    最佳实践 | Best Practices:

    1. 使用内置模板
       - 直接使用内置默认模板，无需额外准备
       - 如需自定义，可使用-y和-s参数指定模板文件

    2. 命名规范
       - locus_tag_prefix使用有意义的名称（如Gene, Target等）
       - report_name描述项目名称

    3. 染色体过滤
       - 使用--chr-prefix只处理需要的染色体
       - 例如：--chr-prefix Chr 只处理Chr01, Chr02...

    4. 并行执行
       - 根据可用资源选择合适的并行度
       - 推荐使用run_all_parallel.sh脚本

    5. 结果整理
       - 每个任务独立运行，互不干扰
       - 完成后可单独检查每个结果

    后续处理 | Post-processing:

    EGAPx运行完成后，每个染色体的输出在对应的output目录中。
    可能需要合并各个染色体的结果：

    ```bash
    # 合并GFF文件
    cat Chr*/output/*.gff > all_genes.gff

    # 合并统计信息
    cat Chr*/output/*summary.txt > all_summary.txt
    ```
    """

    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    egapx_batch_main = _lazy_import_egapx_batch_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['egapx_batch.py']

    # 必需参数 | Required parameters
    args.extend(['-g', genome])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加，减少命令行长度）| Optional parameters (add only when non-default)
    if yaml is not None:
        args.extend(['-y', yaml])

    if script is not None:
        args.extend(['-s', script])

    if egapx != '/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx':
        args.extend(['-e', egapx])

    if chr_prefix:
        args.extend(['-p', chr_prefix])

    if locus_prefix:
        args.extend(['--locus-prefix', locus_prefix])

    if report_name != 'EGAPx':
        args.extend(['--report-name', report_name])

    if short_reads:
        args.extend(['--short-reads', short_reads])

    if long_reads:
        args.extend(['--long-reads', long_reads])

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        egapx_batch_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n[WARNING] 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"[ERROR] 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
