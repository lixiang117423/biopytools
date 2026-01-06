"""
ANNOVAR变异注释|ANNOVAR Variant Annotation Command

"""

import click
import sys
import os


def _lazy_import_annovar_main():
    """延迟加载annovar主函数|Lazy load annovar main function"""
    try:
        from ...annovar.main import main as annovar_main
        return annovar_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='ANNOVAR变异注释工具|ANNOVAR variant annotation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF变异文件|VCF variant file path')
@click.option('--gff3', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF3注释文件|GFF3 annotation file path')
@click.option('--genome', '-G',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组序列文件|Genome sequence file path')
@click.option('--build-ver', '-b',
              required=True,
              help='基因组版本(例如: OV, KY131)|Genome build version identifier (e.g., OV, KY131)')
@click.option('--annovar-path', '-a',
              default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
              show_default=True,
              help='ANNOVAR软件路径|ANNOVAR software installation path')
@click.option('--database-path', '-d',
              default='./database',
              show_default=True,
              type=click.Path(),
              help='ANNOVAR数据库路径|ANNOVAR database path')
@click.option('--output-dir', '-o',
              default='./annovar_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--qual-threshold', '-q',
              type=int,
              default=20,
              show_default=True,
              help='VCF质量阈值|VCF quality filtering threshold')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4']),
              help='运行指定步骤|Run only specified step:\n'
                   '1: GFF3转换|GFF3 conversion\n'
                   '2: 提取序列|Extract sequences\n'
                   '3: VCF处理|VCF processing\n'
                   '4: 变异注释|Variant annotation')
@click.option('--skip-gff-cleaning',
              is_flag=True,
              help='跳过GFF3清理|Skip GFF3 file format cleaning')
@click.option('--skip-gff-fix',
              is_flag=True,
              help='跳过GFF3修复|Skip automatic GFF3 file fixes')
@click.option('--enable-vcf-filter',
              is_flag=True,
              help='启用VCF过滤(默认跳过)|Enable VCF filtering step (skipped by default)')
def annovar(input, gff3, genome, build_ver, annovar_path, database_path,
           output_dir, qual_threshold, step, skip_gff_cleaning,
           skip_gff_fix, enable_vcf_filter):
    """
    ANNOVAR变异注释工具|ANNOVAR Variant Annotation Tool

    使用ANNOVAR进行VCF文件变异注释和功能预测|Use ANNOVAR for VCF variant annotation and functional prediction

    示例|Example: biopytools annovar -i variants.vcf -g annotation.gff3 -G genome.fa -b OV
    """

    # 延迟加载|Lazy loading
    annovar_main = _lazy_import_annovar_main()

    # 构建参数列表|Build argument list
    args = ['annovar.py']

    # 必需参数|Required parameters (按新顺序)
    args.extend(['-i', input])
    args.extend(['-g', gff3])
    args.extend(['-G', genome])
    args.extend(['-b', build_ver])

    # 可选参数|Optional parameters
    if annovar_path != '/share/org/YZWL/yzwl_lixg/software/annovar/annovar':
        args.extend(['-a', annovar_path])

    if database_path != './database':
        args.extend(['-d', database_path])

    if output_dir != './annovar_output':
        args.extend(['-o', output_dir])

    if qual_threshold != 20:
        args.extend(['-q', str(qual_threshold)])

    # 步骤控制|Step control
    if step:
        args.extend(['-s', step])

    # 布尔选项|Boolean options
    if skip_gff_cleaning:
        args.append('--skip-gff-cleaning')

    if skip_gff_fix:
        args.append('--skip-gff-fix')

    if enable_vcf_filter:
        args.append('--enable-vcf-filter')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        annovar_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
