"""
🧬 ANNOVAR变异注释命令 | ANNOVAR Variant Annotation Command
"""

import click
import sys
from ...annovar.main import main as annovar_main


@click.command(short_help = "ANNOVAR变异注释工具")
@click.option('--gff3', '-g',
              required=True,
              type=click.Path(exists=True),
              help='📂 GFF3注释文件路径 | GFF3 annotation file path')
@click.option('--genome', '-f',
              required=True,
              type=click.Path(exists=True),
              help='🧬 基因组序列文件路径 | Genome sequence file path')
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='📄 VCF变异文件路径 | VCF variant file path')
@click.option('--build-ver', '-b',
              required=True,
              help='🏗️ 基因组构建版本标识符 (如: OV, KY131) | '
                   'Genome build version identifier (e.g., OV, KY131)')
@click.option('--annovar-path', '-a',
              default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
              help='🛠️ ANNOVAR软件安装路径 | ANNOVAR software installation path')
@click.option('--database-path', '-d',
              default='./database',
              type=click.Path(),
              help='💾 ANNOVAR数据库路径 | ANNOVAR database path')
@click.option('--output-dir', '-o',
              default='./annovar_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory')
@click.option('--qual-threshold', '-q',
              type=int,
              default=20,
              help='🎯 VCF质量过滤阈值 (默认: 20) | VCF quality filtering threshold (default: 20)')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4']),
              help='🎯 只运行指定步骤 | Run only specified step:\n'
                   '1: 🔄 GFF3转换 | GFF3 conversion\n'
                   '2: 🧬 提取序列 | Extract sequences\n' 
                   '3: 🔍 VCF处理 | VCF processing\n'
                   '4: 📝 变异注释 | Variant annotation')
@click.option('--skip-gff-cleaning',
              is_flag=True,
              help='⏭️ 跳过GFF3文件的格式清理 | Skip GFF3 file format cleaning')
@click.option('--skip-gff-fix',
              is_flag=True,
              help='⏭️ 跳过GFF3文件的自动修复 | Skip automatic GFF3 file fixes')
@click.option('--enable-vcf-filter',
              is_flag=True,
              help='🔍 启用VCF过滤步骤 (默认跳过) | Enable VCF filtering step (skipped by default)')
def annovar(gff3, genome, vcf, build_ver, annovar_path, database_path, 
           output_dir, qual_threshold, step, skip_gff_cleaning, 
           skip_gff_fix, enable_vcf_filter):
    """
    ANNOVAR变异注释工具.
    
    对VCF文件进行基因变异注释，包括GFF3格式转换、序列提取、
    VCF预处理和最终的变异功能注释。
    
    示例 | Examples:
    
    \b
    # 🎯 完整注释流程
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        -o annotation_results
    
    \b
    # 🔧 只运行特定步骤
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --step 1
    
    \b
    # 🔍 启用VCF过滤和自定义路径
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
    # ⏩ 跳过清理步骤
    biopytools annovar \\
        -g clean.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --skip-gff-cleaning \\
        --skip-gff-fix
    """
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'annovar']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', gff3])
    args.extend(['-f', genome])
    args.extend(['-v', vcf])
    args.extend(['-b', build_ver])
    
    # 可选参数（只在非默认值时添加，减少命令行长度）⚙️ | Optional parameters (add only when non-default)
    if annovar_path != '/share/org/YZWL/yzwl_lixg/software/annovar/annovar':
        args.extend(['-a', annovar_path])
    
    if database_path != './database':
        args.extend(['-d', database_path])
    
    if output_dir != './annovar_output':
        args.extend(['-o', output_dir])
    
    if qual_threshold != 20:
        args.extend(['-q', str(qual_threshold)])
    
    # 步骤控制 🎯 | Step control
    if step:
        args.extend(['-s', step])
    
    # 处理选项（布尔标志）🚩 | Processing options (boolean flags)
    if skip_gff_cleaning:
        args.append('--skip-gff-cleaning')
    
    if skip_gff_fix:
        args.append('--skip-gff-fix')
    
    # VCF过滤逻辑 🔍 | VCF filtering logic (important!)
    # 默认情况下 skip_vcf_filter=True，除非明确启用VCF过滤
    # By default skip_vcf_filter=True, unless explicitly enable VCF filtering
    if enable_vcf_filter:
        args.append('--enable-vcf-filter')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        annovar_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv