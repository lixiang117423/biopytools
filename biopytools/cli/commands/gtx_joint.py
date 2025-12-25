"""
🧬 GTX Joint Calling命令生成器 | GTX Joint Calling Command Generator
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_gtx_joint_main():
    """懒加载gtx-joint main函数 | Lazy load gtx-joint main function"""
    try:
        from ...gtx_joint.main import main as gtx_joint_main
        return gtx_joint_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_gtx_exec(file_path):
    """验证GTX可执行文件是否存在（仅在非帮助模式下）| Validate GTX executable (only in non-help mode)"""
    if _is_help_request():
        return file_path
    # 使用默认值（如果未提供）
    if not file_path:
        file_path = '/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx'
    if not os.path.exists(file_path):
        raise click.BadParameter(f"GTX可执行文件不存在 | GTX executable not found: {file_path}")
    if not os.access(file_path, os.X_OK):
        raise click.BadParameter(f"GTX文件不可执行 | GTX file is not executable: {file_path}")
    return file_path


def _validate_reference_file(file_path):
    """验证参考基因组文件是否存在（仅在非帮助模式下）| Validate reference file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("参考基因组文件路径不能为空 | Reference file path cannot be empty")
    if not os.path.exists(file_path):
        raise click.BadParameter(f"参考基因组文件不存在 | Reference file not found: {file_path}")
    # 检查索引文件
    if not os.path.exists(f"{file_path}.fai"):
        raise click.BadParameter(
            f"参考基因组索引不存在 | Reference index not found: {file_path}.fai\n"
            f"请运行 | Please run: samtools faidx {file_path}"
        )
    return file_path


def _validate_gvcf_dir(dir_path):
    """验证GVCF目录是否存在（仅在非帮助模式下）| Validate GVCF directory (only in non-help mode)"""
    if _is_help_request():
        return dir_path
    if not dir_path:
        raise click.BadParameter("GVCF目录路径不能为空 | GVCF directory path cannot be empty")
    if not os.path.exists(dir_path):
        raise click.BadParameter(f"GVCF目录不存在 | GVCF directory not found: {dir_path}")
    if not os.path.isdir(dir_path):
        raise click.BadParameter(f"GVCF路径不是目录 | GVCF path is not a directory: {dir_path}")
    return dir_path


@click.command(short_help="GTX Joint Calling命令生成工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--gtx', '-g',
              default='/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx',
              callback=lambda ctx, param, value: _validate_gtx_exec(value) if value else None,
              help='📂 GTX可执行文件路径 (默认: /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx) | GTX executable path')
@click.option('--ref', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_reference_file(value) if value else None,
              help='🧬 参考基因组文件路径 | Reference genome file path')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_gvcf_dir(value) if value else None,
              help='📂 GVCF文件所在目录 | GVCF files directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📤 输出结果目录 | Output directory')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--tmp-dir', '-T',
              default='./tmp',
              type=click.Path(),
              help='🗑️ 临时目录 | Temporary directory (default: ./tmp)')
@click.option('--script', '-s',
              default='run_gtx_joint.sh',
              help='📄 输出脚本文件名 | Output script filename (default: run_gtx_joint.sh)')
@click.option('--faketime', '-f',
              default='2020-10-20 00:00:00',
              help='⏰ faketime时间 | faketime time (default: 2020-10-20 00:00:00)')
@click.option('--pattern', '-p',
              help='🔍 染色体过滤正则 | Chromosome filter pattern (e.g., "^Chr[0-9]+$")')
@click.option('--window', '-w',
              type=int,
              help='📏 区间大小(bp) | Window size in bp (e.g., 10000000 for 10M)')
def gtx_joint(gtx, ref, input, output, threads, tmp_dir, script, faketime, pattern, window):
    """
    🧬 GTX Joint Calling命令生成工具

    按染色体或区间自动生成GTX joint calling命令脚本，支持大规模多样品
    联合变异检测任务的并行化处理。

    核心功能 | Core Features:
    - 📂 自动扫描GVCF文件并构建输入参数
    - 🧬 从参考基因组索引读取染色体信息
    - 🔄 支持按染色体或按区间拆分任务
    - 📝 生成可直接执行的shell脚本
    - 🔍 支持染色体过滤正则表达式
    - ⏰ 可选的faketime支持
    - 💡 提供多种并行执行建议

    示例 | Examples:

    \b
    # 🚀 基本用法 - 按染色体生成（使用默认GTX路径）
    biopytools gtx-joint \\
        -r genome.fa \\
        -i ./gvcf \\
        -o ./output

    \b
    # 📏 按10M区间拆分（防止内存溢出）
    biopytools gtx-joint \\
        -r genome.fa \\
        -i ./gvcf \\
        -o ./output \\
        -w 10000000

    \b
    # 🔍 只处理主染色体（Chr01-Chr99）
    biopytools gtx-joint \\
        -r genome.fa \\
        -i ./gvcf \\
        -o ./output \\
        -p "^Chr[0-9]+$"

    \b
    # ⚡ 指定自定义GTX路径
    biopytools gtx-joint \\
        -g /custom/path/to/gtx \\
        -r reference.fa \\
        -i ./gvcf_files \\
        -o ./joint_output \\
        -t 24 \\
        -T /tmp \\
        -s my_joint_commands.sh \\
        -f "2020-01-01 00:00:00"

    输出说明 | Output Description:

    生成的脚本文件包含所有GTX joint calling命令，每条命令处理一个
    染色体或区间。脚本示例：
    ```
    faketime '2020-10-20 00:00:00' /path/to/gtx joint \\
        -r genome.fa \\
        -o output/Chr01.joint.vcf.gz \\
        -L Chr01 \\
        --tmp-dir ./tmp \\
        -t 88 \\
        -v sample1.g.vcf.gz -v sample2.g.vcf.gz ...
    ```

    执行方式 | Execution Methods:

    \b
    # 串行执行
    bash run_gtx_joint.sh

    \b
    # GNU Parallel并行执行（4个任务）
    parallel -j 4 < run_gtx_joint.sh

    \b
    # 使用批量提交脚本
    your_batch_submit_script.sh run_gtx_joint.sh

    区间模式后续处理 | Post-processing for Window Mode:

    如果使用-w参数按区间拆分，运行完成后需要合并同一染色体的
    所有区间VCF文件：
    ```
    bcftools concat -o Chr01.merged.vcf.gz Chr01_*.joint.vcf.gz
    ```

    资源建议 | Resource Recommendations:

    - 💾 内存：建议≥32GB/任务，大基因组项目需≥64GB
    - 💽 磁盘：预留输入数据的2-3倍空间
    - ⚡ CPU：根据内存情况调整并行任务数

    故障排除 | Troubleshooting:

    1. "GTX可执行文件不存在"
       - 使用-g参数指定正确的GTX路径
       - 确保GTX文件具有可执行权限

    2. "未找到任何*.g.vcf.gz文件"
       - 检查-i指定的目录是否正确
       - 确认GVCF文件以.g.vcf.gz结尾

    3. "参考基因组索引不存在"
       - 运行：samtools faidx genome.fa

    4. "faketime未找到"
       - 这是警告，不影响功能
       - 程序将不使用faketime
    """

    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    gtx_joint_main = _lazy_import_gtx_joint_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['gtx_joint.py']

    # 必需参数 | Required parameters
    if gtx != '/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx':
        args.extend(['-g', gtx])
    args.extend(['-r', ref])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加，减少命令行长度）| Optional parameters (add only when non-default)
    if threads != 88:
        args.extend(['-t', str(threads)])

    if tmp_dir != './tmp':
        args.extend(['-T', tmp_dir])

    if script != 'run_gtx_joint.sh':
        args.extend(['-s', script])

    if faketime != '2020-10-20 00:00:00':
        args.extend(['-f', faketime])

    if pattern:
        args.extend(['-p', pattern])

    if window:
        args.extend(['-w', str(window)])

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        gtx_joint_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n⚠️ 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
