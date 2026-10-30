"""Selected CLI help contract tests."""

from click.testing import CliRunner

from biopytools.cli.commands.allelic_genes import allelic_genes
from biopytools.cli.commands.phylo_selector import phylo_selector
from biopytools.cli.commands.vcf2gwas import vcf2gwas


def test_deprecated_allelic_genes_command_loads():
    result = CliRunner().invoke(allelic_genes, ["-h"])
    assert result.exit_code == 0
    assert "DEPRECATED" in result.output


def test_phylo_selector_accepts_short_help():
    result = CliRunner().invoke(phylo_selector, ["-h"])
    assert result.exit_code == 0
    assert "Phylogenetic Tree Sample Selector" in result.output


def test_vcf2gwas_help_does_not_probe_conda(monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("conda should not be resolved for help")

    monkeypatch.setattr(
        "biopytools.cli.commands.vcf2gwas.conda_env_run_prefix",
        fail_if_called,
    )
    result = CliRunner().invoke(vcf2gwas, ["-h"])
    assert result.exit_code == 0
    assert "vcf2gwas GWAS Analysis Tool" in result.output
