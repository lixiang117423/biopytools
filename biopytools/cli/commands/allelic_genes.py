"""兼容已废弃命令|Compatibility command for deprecated allelic-genes."""

import click


@click.command(
    short_help='[已废弃,请使用jcvi allelic]|[DEPRECATED, use jcvi allelic]',
    context_settings=dict(help_option_names=['-h', '--help']),
)
def allelic_genes():
    """DEPRECATED: use `biopytools jcvi allelic`|请改用 jcvi allelic."""
    click.echo(
        "已废弃|DEPRECATED: use `biopytools jcvi allelic` instead.",
        err=True,
    )
    raise click.exceptions.Exit(0)


__all__ = ["allelic_genes"]
