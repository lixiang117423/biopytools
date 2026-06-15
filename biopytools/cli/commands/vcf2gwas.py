"""vcf2gwas GWAS分析命令|vcf2gwas GWAS Analysis Command"""

import subprocess
import sys
import click


def _print_help():
    """打印帮助信息|Print help"""
    click.echo("""
Usage: biopytools vcf2gwas [OPTIONS] vcf2gwas_args...

  vcf2gwas GWAS分析工具 - 所有参数直接传递给vcf2gwas
  vcf2gwas GWAS Analysis Tool - all arguments passed through

Options:
  --vcf2gwas-env TEXT  conda环境名|conda env name (default: vcf2gwas_v.0.8.9)
  -h, --help           显示本帮助|Show this help

示例|Examples:
  biopytools vcf2gwas -v input.vcf.gz -pf phenotype.csv -p 1 -lmm
  biopytools vcf2gwas -v input.vcf.gz -pf pheno.csv -cf cov.csv -c 1 -p 1 -lmm -P 3
""")


@click.command(
    short_help='vcf2gwas GWAS分析工具|vcf2gwas GWAS Analysis Tool',
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
        allow_interspersed_args=False,
        help_option_names=[],
        max_content_width=120,
    ),
)
@click.option('--vcf2gwas-env', default='vcf2gwas_v.0.8.9',
              help='[STR] conda环境名|conda env name (default: vcf2gwas_v.0.8.9)')
@click.pass_context
def vcf2gwas(ctx, vcf2gwas_env):
    # 无参数时显示包装器帮助|Show wrapper help when no args
    if not ctx.args:
        _print_help()
        return

    full_cmd = ['conda', 'run', '-n', vcf2gwas_env, '--no-capture-output',
                'vcf2gwas'] + list(ctx.args)

    try:
        result = subprocess.call(full_cmd)
        sys.exit(result)
    except Exception as e:
        click.echo(f"执行错误|Execution error: {e}", err=True)
        sys.exit(1)


__all__ = ['vcf2gwas']
