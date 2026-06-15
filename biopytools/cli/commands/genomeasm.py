"""
 CLI|Genome Assembly Pipeline CLI Wrapper
Click--help
"""

import click
import sys
import os
from pathlib import Path

def _lazy_import_assembler_main():
    """|Lazy load the genome assembler main function"""
    try:
        # main.py
        #  genomeasm/main.py  genomeasm/cli.py,  from .main import main
        from ...genomeasm import main as assembler_main
        return assembler_main
    except ImportError as e:
        click.echo(f"Import Error: ", err=True)
        click.echo(f"  Details: {e}", err=True)
        sys.exit(1)

def _is_help_request():
    """|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)

def _validate_input_dir(ctx, param, value):
    """| Validate input directory (only in non-help mode)"""
    if not _is_help_request() and not Path(value).is_dir():
        raise click.BadParameter(f"|Input directory does not exist or is not a directory: {value}")
    return value

@click.command(
    name="assemble",
    short_help=" ",
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
    help="""
      ()|Genome Assembly Tool (Multi-data Integration Version)

    Examples:
    
    \b
    # HiFi
    biopytools assemble -i raw_data/ -o assembly_results/
    
    \b
    # HiFi + Hi-C
    biopytools assemble -i data/ -o results/ --hic-strategy complete_juicer
    
    \b
    # 
    biopytools assemble -i input/ -o output/ -n my_genome --genome-size 3g -t 64
    
    \b
    # Hi-C
    biopytools assemble -i data/ -o results/ --hic-strategy simplified_salsa2

    Supported Data Types:
      - HiFi:  ()|High-accuracy long reads (required)
      - Hi-C:Chromosome conformation capture data
      - ONT: Oxford Nanopore|Oxford Nanopore long reads
      - NGS: Illumina|Illumina short reads

     Hi-C|Hi-C Processing Strategies:
      - complete_juicer: Juicer + 3D-DNA ()
      - standard_3ddna: 3D-DNA ()
      - simplified_salsa2: SALSA2 ()
    """
)
# ---Required Arguments ---
@click.option('-i', '--input-dir', 
              required=True, 
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              callback=_validate_input_dir,
              help=' ()')

# ---Basic Arguments ---
@click.option('-o', '--output-dir', default='./assembly_output', type=click.Path(), help='')
@click.option('-n', '--project-name', default='genome_assembly', help='')
@click.option('-t', '--threads', type=int, default=12, help='')

# --- Hi-C|Hi-C Arguments ---
@click.option('--hic-strategy', 
              type=click.Choice(['complete_juicer', 'standard_3ddna', 'simplified_salsa2']),
              default='complete_juicer', help='Hi-C')
@click.option('--restriction-enzyme', 
              type=click.Choice(['MboI', 'DpnII', 'HindIII', 'EcoRI']),
              default='MboI', help='')
@click.option('--min-contig-size', type=int, default=15000, help='contig')
@click.option('--edit-rounds', type=int, default=2, help='3D-DNA')

# ---Assembly Arguments ---
@click.option('--genome-size', default='3g', help=' (e.g., 3g, 500m)')
@click.option('--species-type', 
              type=click.Choice(['diploid', 'haploid', 'polyploid']),
              default='diploid', help='')
@click.option('--telomere-motif', default='CCCTAA', help='motif')
@click.option('--purge-level', type=click.Choice(['0', '1', '2', '3']), default='1', help='Purging')
@click.option('--purge-max', type=int, default=80, help='Purging')
@click.option('--similarity-threshold', type=float, default=0.75, help='')
@click.option('--n-haplotypes', type=int, default=2, help='')

# ---Quality Control Arguments ---
# @click.option('--skip-fastqc', default=True, help='FastQC ()|Skip FastQC quality check (default: skip to save time)')
@click.option('--skip-fastqc', default=True, help='FastQC ()')
@click.option('--min-hifi-coverage', type=int, default=30, help='HiFi')
@click.option('--min-hic-coverage', type=int, default=50, help='Hi-C')
@click.option('--min-mapping-rate', type=float, default=0.7, help='')
@click.option('--busco-lineage', default='auto', help='BUSCO')

# ---Tool Paths ---
@click.option('--hifiasm-path', default='hifiasm', help='Hifiasm')
@click.option('--bwa-path', default='bwa', help='BWA')
@click.option('--samtools-path', default='samtools', help='Samtools')
@click.option('--juicer-path', default='juicer.sh', help='Juicer')
@click.option('--pipeline-3ddna', default='3d-dna/run-asm-pipeline.sh', help='3D-DNA pipeline')
@click.option('--juicer-tools', default='juicer_tools.jar', help='Juicer tools JAR')
@click.option('--salsa2-path', default='run_pipeline.py', help='SALSA2')
def genomeasm(**kwargs):
    """
     Click
    
    
    `argparse`  `sys.argv` 
    """
    #  
    assembler_main = _lazy_import_assembler_main()
    
    #  main
    # argparse
    args = ['genomeasm.py']
    
    # 
    defaults = {param.name: param.default for param in genomeasm.params}
    
    for key, value in kwargs.items():
        # click'-''_'
        cli_option = '--' + key.replace('_', '-')
        
        # 
        # 'input_dir'  (None)
        if value is not None and value != defaults.get(key):
            args.extend([cli_option, str(value)])

    #  sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    # click.echo(" ...")
    # click.echo(f"   : {' '.join(sys.argv)}")
    
    try:
        #  main
        assembler_main()
    except SystemExit as e:
        #   (e.g., sys.exit(0) or sys.exit(1))
        # click
        if e.code != 0:
            click.secho(f" : {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        #  
        click.secho(f" : {e}", fg='red', bold=True, err=True)
        # 
        # import traceback
        # click.echo(traceback.format_exc(), err=True)
        sys.exit(1)
    finally:
        # sys.argv
        sys.argv = original_argv

if __name__ == '__main__':
    genomeasm()