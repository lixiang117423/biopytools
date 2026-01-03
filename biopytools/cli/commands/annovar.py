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


@click.command(short_help="ANNOVAR",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--gff3', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF3注释文件|GFF3 annotation file path')
@click.option('--genome', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组序列文件|Genome sequence file path')
@click.option('--vcf', '-v',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF变异文件|VCF variant file path')
@click.option('--build-ver', '-b',
              required=True,
              help='基因组版本(例如: OV, KY131)|Genome build version identifier (e.g., OV, KY131)')
@click.option('--annovar-path', '-a',
              default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
              help='ANNOVAR软件路径|ANNOVAR software installation path')
@click.option('--database-path', '-d',
              default='./database',
              type=click.Path(),
              help='ANNOVAR数据库路径|ANNOVAR database path')
@click.option('--output-dir', '-o',
              default='./annovar_output',
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--qual-threshold', '-q',
              type=int,
              default=20,
              help='VCF质量阈值(默认: 20)|VCF quality filtering threshold (default: 20)')
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
def annovar(gff3, genome, vcf, build_ver, annovar_path, database_path, 
           output_dir, qual_threshold, step, skip_gff_cleaning, 
           skip_gff_fix, enable_vcf_filter):
    """
    ANNOVAR变异注释流程

    ANNOVAR功能强大的变异注释工具
    支持VCF文件变异注释和功能预测

    使用示例|Examples:

    \b
    # 基本用法
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        -o annotation_results

    \b
    # 只运行步骤1
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --step 1

    \b
    # 启用VCF过滤
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        -a /path/to/annovar \\
        -d /path/to/database \\
        --enable-vcf-filter \\
        --qual-threshold 30

    \b
    # 跳过GFF3处理步骤
    biopytools annovar \\
        -g clean.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --skip-gff-cleaning \\
        --skip-gff-fix
    """

    # 延迟加载: 仅在实际调用时导入
    annovar_main = _lazy_import_annovar_main()

    # 构建主函数的参数列表
    args = ['annovar.py']

    # 必需参数
    args.extend(['-g', gff3])
    args.extend(['-f', genome])
    args.extend(['-v', vcf])
    args.extend(['-b', build_ver])

    # 可选参数(仅在非默认值时添加)
    if annovar_path != '/share/org/YZWL/yzwl_lixg/software/annovar/annovar':
        args.extend(['-a', annovar_path])

    if database_path != './database':
        args.extend(['-d', database_path])

    if output_dir != './annovar_output':
        args.extend(['-o', output_dir])

    if qual_threshold != 20:
        args.extend(['-q', str(qual_threshold)])

    # 步骤控制
    if step:
        args.extend(['-s', step])

    # 处理选项(布尔标志)
    if skip_gff_cleaning:
        args.append('--skip-gff-cleaning')

    if skip_gff_fix:
        args.append('--skip-gff-fix')

    # VCF过滤逻辑(重要!)
    # 默认跳过VCF过滤,除非明确启用VCF过滤
    if enable_vcf_filter:
        args.append('--enable-vcf-filter')

    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数
        annovar_main()
    except SystemExit as e:
        # 处理正常程序退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv