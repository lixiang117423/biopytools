"""
🧬 VCF文件样品名称重命名工具 | VCF Sample Name Renamer
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_vcf_renamer_main():
    """懒加载vcf-renamer main函数 | Lazy load vcf-renamer main function"""
    try:
        from ...vcf_renamer.main import main as vcf_renamer_main
        return vcf_renamer_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_vcf(file_path):
    """验证输入VCF文件是否存在（仅在非帮助模式下）| Validate input VCF file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("输入VCF文件路径不能为空 | Input VCF path cannot be empty")
    if not os.path.exists(file_path):
        raise click.BadParameter(f"输入VCF文件不存在 | Input VCF not found: {file_path}")
    if not file_path.endswith(('.vcf.gz', '.vcf')):
        raise click.BadParameter(f"输入文件必须是VCF格式 | Input must be VCF format: {file_path}")
    return file_path


@click.command(short_help="VCF样品名称重命名工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_vcf(value) if value else None,
              help='📂 输入VCF文件路径 | Input VCF file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📤 输出VCF文件路径 | Output VCF file path')
@click.option('--prefix', '-p',
              default='S',
              help='🏷️ 新样品名前缀 (默认: S) | New sample name prefix (default: S)')
@click.option('--mapping', '-m',
              type=click.Path(),
              help='📋 映射文件路径 | Mapping file path')
@click.option('--no-mapping',
              is_flag=True,
              help='🗑️ 不保留映射文件 | Don\'t keep mapping file')
def vcf_renamer(input, output, prefix, mapping, no_mapping):
    """
    🧬 VCF文件样品名称重命名工具

    重命名VCF文件中的样品名称（如 S1, S2, S3...），防止样品名称被
    软件截断的问题。特别适用于GATK、GTX等对样品名称长度有限制的
    软件。

    核心功能 | Core Features:
    - 📋 自动提取VCF文件中的所有样品名称
    - 🏷️ 按顺序重命名样品（S1, S2, S3...或自定义前缀）
    - 📝 生成新旧样品名称映射表
    - 🔄 使用bcftools reheader进行重命名
    - 📇 自动创建VCF索引文件
    - 🧹 可选清理映射文件

    使用场景 | Use Cases:
    - ✅ 样品名称过长被软件截断
    - ✅ 需要简化样品名称便于处理
    - ✅ 批量重命名多个样品
    - ✅ 为后续分析准备标准化的样品名

    示例 | Examples:

    \b
    # 🚀 基本用法（使用默认前缀S）
    biopytools vcf-renamer \\
        -i variation.filtered.snp.vcf.gz \\
        -o variation.renamed.vcf.gz

    \b
    # 🏷️ 自定义样品名前缀为Sample
    biopytools vcf-renamer \\
        -i input.vcf.gz \\
        -o output.vcf.gz \\
        -p Sample

    \b
    # 📋 指定映射文件路径
    biopytools vcf-renamer \\
        -i input.vcf.gz \\
        -o output.vcf.gz \\
        -m sample_mapping.txt

    \b
    # 🗑️ 不保留映射文件
    biopytools vcf-renamer \\
        -i input.vcf.gz \\
        -o output.vcf.gz \\
        --no-mapping

    输出说明 | Output Description:

    1. **重命名后的VCF文件** - 样品名称已更新为S1, S2, S3...
    2. **VCF索引文件** - 自动生成.tbi索引文件
    3. **映射文件** - 保存新旧样品名称对应关系（可选）

    映射文件格式 | Mapping File Format:
    ```
    Original_Sample_Name_1    S1
    Original_Sample_Name_2    S2
    Original_Sample_Name_3    S3
    ...
    ```

    依赖软件 | Dependencies:
    - bcftools - 必需，用于重命名和索引操作

    安装依赖 | Install Dependencies:
    ```bash
    conda install -c bioconda bcftools
    ```

    故障排除 | Troubleshooting:

    1. "未找到bcftools命令"
       - 安装bcftools: conda install -c bioconda bcftools
       - 或使用: sudo apt install bcftools

    2. "输入VCF文件不存在"
       - 检查-i参数指定的文件路径
       - 确认文件格式为.vcf.gz或.vcf

    3. "提取样品名称失败"
       - 确认VCF文件格式正确
       - 检查VCF文件是否损坏
       - 尝试运行: bcftools view -h input.vcf.gz

    4. "重命名失败"
       - 检查输出目录是否有写权限
       - 确认磁盘空间充足

    最佳实践 | Best Practices:

    1. 📋 保留映射文件：建议保留映射文件，方便后续追溯原始样品名
    2. 🏷️ 前缀选择：根据项目选择有意义的前缀，如"Sample"、"Pop"等
    3. 📦 批量处理：可结合shell脚本批量处理多个VCF文件
    4. ✅ 验证结果：重命名后使用bcftools query -l验证新样品名

    验证结果 | Verify Results:
    ```bash
    # 查看重命名后的样品名称
    bcftools query -l output.vcf.gz

    # 应该显示:
    # S1
    # S2
    # S3
    # ...
    ```

    后续处理 | Post-processing:

    如果需要恢复原始样品名称，可以反向使用映射文件：
    ```bash
    # 反向映射文件（手动创建）
    S1    Original_Sample_Name_1
    S2    Original_Sample_Name_2

    # 重新使用vcf-renamer或bcftools reheader恢复
    ```
    """

    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    vcf_renamer_main = _lazy_import_vcf_renamer_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['vcf_renamer.py']

    # 必需参数 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加，减少命令行长度）| Optional parameters (add only when non-default)
    if prefix != 'S':
        args.extend(['-p', prefix])

    if mapping:
        args.extend(['-m', mapping])

    if no_mapping:
        args.append('--no-mapping')

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        vcf_renamer_main()
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
