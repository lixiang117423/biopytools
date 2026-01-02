"""
🧬 ANNOVAR | ANNOVAR Variant Annotation Command

"""

import click
import sys
import os


def _lazy_import_annovar_main():
    """annovar main | Lazy load annovar main function"""
    try:
        from ...annovar.main import main as annovar_main
        return annovar_main
    except ImportError as e:
        click.echo(f"  | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """ | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f" | File does not exist: {file_path}")
    return file_path


@click.command(short_help="ANNOVAR",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--gff3', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help=' GFF3 | GFF3 annotation file path')
@click.option('--genome', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧬  | Genome sequence file path')
@click.option('--vcf', '-v',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help=' VCF | VCF variant file path')
@click.option('--build-ver', '-b',
              required=True,
              help='  (: OV, KY131) | '
                   'Genome build version identifier (e.g., OV, KY131)')
@click.option('--annovar-path', '-a',
              default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
              help=' ANNOVAR | ANNOVAR software installation path')
@click.option('--database-path', '-d',
              default='./database',
              type=click.Path(),
              help=' ANNOVAR | ANNOVAR database path')
@click.option('--output-dir', '-o',
              default='./annovar_output',
              type=click.Path(),
              help='  | Output directory')
@click.option('--qual-threshold', '-q',
              type=int,
              default=20,
              help=' VCF (: 20) | VCF quality filtering threshold (default: 20)')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4']),
              help='  | Run only specified step:\n'
                   '1:  GFF3 | GFF3 conversion\n'
                   '2: 🧬  | Extract sequences\n' 
                   '3:  VCF | VCF processing\n'
                   '4:   | Variant annotation')
@click.option('--skip-gff-cleaning',
              is_flag=True,
              help='⏭ GFF3 | Skip GFF3 file format cleaning')
@click.option('--skip-gff-fix',
              is_flag=True,
              help='⏭ GFF3 | Skip automatic GFF3 file fixes')
@click.option('--enable-vcf-filter',
              is_flag=True,
              help=' VCF () | Enable VCF filtering step (skipped by default)')
def annovar(gff3, genome, vcf, build_ver, annovar_path, database_path, 
           output_dir, qual_threshold, step, skip_gff_cleaning, 
           skip_gff_fix, enable_vcf_filter):
    """
    ANNOVAR.
    
    VCFGFF3
    VCF
    
     | Examples:
    
    \b
    #  
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        -o annotation_results
    
    \b
    #  
    biopytools annovar \\
        -g annotation.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --step 1
    
    \b
    #  VCF
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
    # ⏩ 
    biopytools annovar \\
        -g clean.gff3 \\
        -f genome.fa \\
        -v variants.vcf \\
        -b OV \\
        --skip-gff-cleaning \\
        --skip-gff-fix
    """
    
    #   | Lazy loading: import only when actually called
    annovar_main = _lazy_import_annovar_main()
    
    # main  | Build argument list for original main function
    args = ['annovar.py']
    
    #   | Required parameters
    args.extend(['-g', gff3])
    args.extend(['-f', genome])
    args.extend(['-v', vcf])
    args.extend(['-b', build_ver])
    
    #  | Optional parameters (add only when non-default)
    if annovar_path != '/share/org/YZWL/yzwl_lixg/software/annovar/annovar':
        args.extend(['-a', annovar_path])
    
    if database_path != './database':
        args.extend(['-d', database_path])
    
    if output_dir != './annovar_output':
        args.extend(['-o', output_dir])
    
    if qual_threshold != 20:
        args.extend(['-q', str(qual_threshold)])
    
    #   | Step control
    if step:
        args.extend(['-s', step])
    
    #  | Processing options (boolean flags)
    if skip_gff_cleaning:
        args.append('--skip-gff-cleaning')
    
    if skip_gff_fix:
        args.append('--skip-gff-fix')
    
    # VCF  | VCF filtering logic (important!)
    #  skip_vcf_filter=TrueVCF
    # By default skip_vcf_filter=True, unless explicitly enable VCF filtering
    if enable_vcf_filter:
        args.append('--enable-vcf-filter')
    
    # sys.argv  | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # main  | Call original main function
        annovar_main()
    except SystemExit as e:
        #   | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"  | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv