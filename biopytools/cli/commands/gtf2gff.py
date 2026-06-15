"""
GTF到GFF转换命令包装器|GTF to GFF Converter Command Wrapper
"""

import click


@click.command()
@click.option('-i', '--input', required=True, help='输入GTF文件|Input GTF file')
@click.option('-o', '--output', required=True, help='输出GFF文件|Output GFF file')
@click.option('--remove-introns', is_flag=True, help='移除intron特征|Remove intron features')
@click.option('--keep-all-attributes', is_flag=True, help='保留所有原始属性|Keep all original attributes')
@click.option('--no-clean', is_flag=True, help='不清理属性字段|Do not clean attributes')
@click.option('-p', '--prefix', help='ID前缀|ID prefix (e.g., CDRT, AGIS)')
@click.option('-s', '--species', help='物种缩写|Species abbreviation (e.g., Ov, Os)')
@click.option('-t', '--threads', default=12, help='线程数|Number of threads [default: 12]')
def gtf2gff(input, output, remove_introns, keep_all_attributes, no_clean, prefix, species, threads):
    """GTF到GFF文件转换工具|GTF to GFF File Converter

    示例|Examples: biopytools gtf2gff -i input.gtf -o output.gff
    """
    # 延迟导入模块|Lazy import module
    from biopytools.gtf_to_gff.main import main as gtf_to_gff_main

    import sys
    # 构建参数列表|Build arguments list
    sys.argv = [
        'gtf_to_gff',
        '--input', input,
        '--output', output,
        '--threads', str(threads)
    ]

    if remove_introns:
        sys.argv.append('--remove-introns')

    if keep_all_attributes:
        sys.argv.append('--keep-all-attributes')

    if no_clean:
        sys.argv.append('--no-clean')

    if prefix:
        sys.argv.extend(['--prefix', prefix])

    if species:
        sys.argv.extend(['--species', species])

    # 调用主函数|Call main function
    gtf_to_gff_main()
