"""
GFF|GFF Format Conversion Command
--help
"""

import click
import sys
import os


def _lazy_import_gff_main():
    """GFF main|Lazy load GFF main function"""
    try:
        from ...gffconverter.main import main as gff_main
        return gff_main
    except ImportError as e:
        click.echo(f"|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"|File does not exist: {file_path}")
    return file_path


@click.command(short_help='GFFID',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='GFF|Input GFF file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='GFF|Output GFF file path')
@click.option('--species-name', '-s',
              required=True,
              type=str,
              help=' (: OV53)|Species name (e.g., OV53)')
@click.option('--species-prefix', '-p',
              required=True,
              type=str,
              help=' (: Ov)|Species prefix (e.g., Ov)')
@click.option('--start-num',
              type=int,
              default=10,
              help='Starting number for gene numbering (default: 10)|Starting number for gene numbering (default: 10)')
@click.option('--step',
              type=int,
              default=10,
              help='Step size for gene numbering (default: 10)|Step size for gene numbering (default: 10)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              help='Number of threads (default: 88)|Number of threads (default: 88)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='Verbose output mode|Verbose output mode')
@click.option('--keep-intermediate',
              is_flag=True,
              help='Keep intermediate files|Keep intermediate files')
@click.option('--show-sample',
              type=int,
              metavar='N',
              help='N|Show N conversion samples and exit')
def gffconverter(input, output, species_name, species_prefix, start_num, step, 
                 threads, verbose, keep_intermediate, show_sample):
    """
    GFF|GFF Format Conversion Tool
    
    GFF/GFF3ID
    
    
    
   Features:
    - GFF/GFF3
    - ID
    - 
    - 
    - 
    - 
    - 
    
   Conversion Pipeline:
    1. GFF
    2. 
    3. ID
    4. 
    5. GFF
    
    ID|ID Naming Convention:
    ID: {species_prefix}g{padded_number}
    : "Ov" + 10 + 10 = Ovg010, Ovg020, Ovg030...
    
   Use Cases:
    - 
    - 
    - 
    - ID
    - 
    - 
    
   Examples:
    
    \b
    # 
    biopytools renamegff -i input.gff -o output.gff \\
        -s OV53 -p Ov
    
    \b
    # 
    biopytools renamegff -i annotation.gff -o standardized.gff \\
        -s "Ecoli_K12" -p Ec --start-num 100 --step 5
    
    \b
    # 
    biopytools renamegff -i data.gff3 -o output.gff \\
        -s Species01 -p Sp --show-sample 10
    
    \b
    # 
    biopytools renamegff -i complex.gff -o processed.gff \\
        -s OV53 -p Ov --verbose --keep-intermediate
    
    \b
    # 
    biopytools renamegff -i large_genome.gff -o converted.gff \\
        -s LargeGenome -p Lg -t 64 --start-num 1 --step 1
    
    \b
    # 
    biopytools renamegff -i species_A.gff -o std_species_A.gff \\
        -s SpeciesA -p SpA --start-num 1000 --step 10
    biopytools renamegff -i species_B.gff -o std_species_B.gff \\
        -s SpeciesB -p SpB --start-num 2000 --step 10
    
   Input File Requirements:
    
    GFF/GFF3:
    - GFF3 () GFF2
    - 
    - 
    -  (9)
    -  (.gz)
    
    :
    - gene: 
    - mRNA/transcript:  ()
    - CDS/exon:  ()
    
    GFF3:
    ##gff-version 3
    chr1    source  gene    1000    5000    .    +    .    ID=gene-001;Name=hypothetical_protein
    chr1    source  mRNA    1000    5000    .    +    .    ID=mRNA-001;Parent=gene-001
    chr1    source  CDS     1200    4800    .    +    0    ID=CDS-001;Parent=mRNA-001
    
   Parameter Details:
    
    :
    --input: GFF.gff.gff3.gtf
    --output: GFF3
    --species-name: 
    --species-prefix: ID
    
    :
    --start-num: 10
    --step: 10
        - 10
        - (1-5)
        - (50-100)
    
    :
    --threads: 
    
    :
    --verbose: 
    --keep-intermediate: 
    --show-sample: 
    
   Output Files:
    
    :
    - {output_file}: GFF3
    - {output_file}.log: 
    - conversion_stats.txt:  (verbose)
    
     (--keep-intermediate):
    - parsed_features.tmp: 
    - id_mappings.txt: IDID
    - validation_report.txt: 
    
    :
    - GFF3
    - 
    - ID
    - 
    - 
    
    ID|ID Generation Rules:
    
    ID: {prefix}g{number}
    - prefix: 
    - g: 
    - number: 
    
    :
    prefix="Ov", start_num=10, step=10:
    - 1: Ovg010
    - 2: Ovg020  
    - 3: Ovg030
    
    prefix="Ec", start_num=1, step=1:
    - 1: Ecg001
    - 2: Ecg002
    - 3: Ecg003
    
    ID:
    - mRNA: {gene_id}.t1, {gene_id}.t2 ()
    - CDS: {mrna_id}.cds
    - protein: {gene_id}.p1, {gene_id}.p2
    
   Performance & System Requirements:
    
    :
    - Python 3.7+
    - pathlib, argparse, logging
    - Python
    
    :
    - RAM: 2GB(>1GB)8GB+
    - : 2
    - CPU: 
    
    :
    - (<10MB): 
    - (10MB-100MB): 1-10
    - (100MB-1GB): 10-60
    - (>1GB): 
    
    :
    - : ~100MB
    - : ~2-3
    - : ~50MB
    
   Troubleshooting:
    
    :
    1. "": 
    2. "": GFF
    3. "": 
    4. "": 
    5. "": UTF-8
    
    :
    - 
    - 9ID
    - 
    - ID
    
    :
    - --show-sample
    - --verbose
    - 
    - 
    
    Best Practices:
    
    1⃣  :
       -  GFF
       -  
       -  
       -  
    
    2⃣  :
       -  
       -  ID
       -  
       - 🧵 
    
    3⃣  :
       -  --show-sample
       -  
       -  ID
       -  GFF3
    
    4⃣  :
       -  
       -  
       -  
       -  
    
    Output Validation:
    
    🤖 :
    -  GFF3
    -  ID
    -  
    -  
    -  
    
     :
    -  ID
    -  
    -  -
    -  
    -  
    
    Citation & References:
    
     :
    -  GFF3 Format Specification: Sequence Ontology Project
    -  Gene Ontology: http://geneontology.org/
    -  FAIR Data Principles: Wilkinson et al. (2016)
    """
    
    #Lazy loading: import only when actually called
    gff_main = _lazy_import_gff_main()
    
    # main|Build argument list for original main function
    args = ['gffconverter.py']
    
    #Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-s', species_name])
    args.extend(['-p', species_prefix])
    
    #|Optional parameters (add only when non-default)
    if start_num != 10:
        args.extend(['--start-num', str(start_num)])
    
    if step != 10:
        args.extend(['--step', str(step)])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    #Boolean options
    if verbose:
        args.append('--verbose')
    
    if keep_intermediate:
        args.append('--keep-intermediate')
    
    #Special parameters
    if show_sample is not None:
        args.extend(['--show-sample', str(show_sample)])
    
    # sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # main|Call original main function
        gff_main()
    except SystemExit as e:
        #Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n GFF|GFF format conversion interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f" GFF|GFF format conversion failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv