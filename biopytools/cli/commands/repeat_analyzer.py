"""
 RepeatAnalyzerCLI|RepeatAnalyzer CLI Command

"""

import click
import sys
import os


def _lazy_import_repeat_main():
    """repeat analyzer main|Lazy load repeat analyzer main function"""
    try:
        from ...repeat_analyzer.main import main as repeat_main
        return repeat_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_file(file_path):
    """| Validate input genome file existence (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(file_path):
            raise click.BadParameter(f"|Input genome file does not exist: {file_path}")
        if not os.path.isfile(file_path):
            raise click.BadParameter(f"|Input path is not a file: {file_path}")
    return file_path


def _validate_output_dir(dir_path):
    """| Validate output directory path (only in non-help mode)"""
    if not _is_help_request():
        # 
        parent_dir = os.path.dirname(os.path.abspath(dir_path))
        if parent_dir and not os.path.exists(parent_dir):
            raise click.BadParameter(f"|Parent directory of output does not exist: {parent_dir}")
    return dir_path


@click.command(short_help='',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_file(value) if value else None,
              help=' FASTA|Input genome FASTA file path')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_output_dir(value) if value else None,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              help='  (: 88)|Number of threads (default: 88)')
@click.option('--skip-modeler',
              is_flag=True,
              help=' RepeatModeler|Skip RepeatModeler step')
@click.option('--skip-ltr',
              is_flag=True,
              help=' LTR|Skip LTR analysis step')
@click.option('--repeatmodeler-path',
              default='RepeatModeler',
              help=' RepeatModeler (: RepeatModeler)|RepeatModeler program path (default: RepeatModeler)')
@click.option('--ltr-finder-path',
              default='ltr_finder',
              help=' LTR_FINDER (: ltr_finder)|LTR_FINDER program path (default: ltr_finder)')
@click.option('--ltrharvest-path',
              default='gt ltrharvest',
              help=' LTRharvest (: gt ltrharvest)|LTRharvest program path (default: gt ltrharvest)')
@click.option('--ltr-retriever-path',
              default='LTR_retriever',
              help=' LTR_retriever (: LTR_retriever)|LTR_retriever program path (default: LTR_retriever)')
@click.option('--repeatmasker-path',
              default='RepeatMasker',
              help=' RepeatMasker (: RepeatMasker)|RepeatMasker program path (default: RepeatMasker)')
@click.option('--tesorter-path',
              default='TEsorter',
              help=' TEsorter (: TEsorter)|TEsorter program path (default: TEsorter)')
def repeat_analyzer(input, output, threads, skip_modeler, skip_ltr,
                   repeatmodeler_path, ltr_finder_path, ltrharvest_path,
                   ltr_retriever_path, repeatmasker_path, tesorter_path):
    """
     
    
    RepeatModelerLTRRepeatMaskerTEsorter
    EDTA
    
    Examples:
    
    \b
    #  
    biopytools repeat-analyzer -i genome.fasta -o repeat_results
    
    \b
    #  
    biopytools repeat-analyzer -i large_genome.fa -o results -t 64
    
    \b
    #  RepeatModeler
    biopytools repeat-analyzer -i genome.fasta -o results --skip-modeler
    
    \b
    #  RepeatMaskerTEsorter
    biopytools repeat-analyzer -i genome.fa -o results \\
        --skip-modeler --skip-ltr
    
    \b
    #  
    biopytools repeat-analyzer -i genome.fasta -o results \\
        --repeatmasker-path ~/.local/bin/RepeatMasker \\
        --tesorter-path ~/.local/bin/TEsorter -t 32
    """
    
    # Lazy loading: import only when actually called
    repeat_main = _lazy_import_repeat_main()
    
    # main|Build argument list for original main function
    args = ['repeat_analyzer.py']  # biopytools
    
    # Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    
    #Optional parameters (add only when non-default)
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    # Boolean options
    if skip_modeler:
        args.append('--skip-modeler')
    
    if skip_ltr:
        args.append('--skip-ltr')
    
    #Tool path parameters (add only when non-default)
    if repeatmodeler_path != 'RepeatModeler':
        args.extend(['--repeatmodeler-path', repeatmodeler_path])
    
    if ltr_finder_path != 'ltr_finder':
        args.extend(['--ltr-finder-path', ltr_finder_path])
        
    if ltrharvest_path != 'gt ltrharvest':
        args.extend(['--ltrharvest-path', ltrharvest_path])
        
    if ltr_retriever_path != 'LTR_retriever':
        args.extend(['--ltr-retriever-path', ltr_retriever_path])
    
    if repeatmasker_path != 'RepeatMasker':
        args.extend(['--repeatmasker-path', repeatmasker_path])
    
    if tesorter_path != 'TEsorter':
        args.extend(['--tesorter-path', tesorter_path])
    
    # sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # main|Call original main function
        repeat_main()
    except SystemExit as e:
        # Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv