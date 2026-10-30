"""JOSS submission packaging checks."""

import json
import re
import tomllib
from pathlib import Path

import biopytools
import yaml
from click.testing import CliRunner
from biopytools.cli.main import cli
from biopytools.cli.main import COMMAND_REGISTRY


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "paper.md"
BIB = ROOT / "paper" / "paper.bib"
README = ROOT / "README.md"


def test_joss_paper_has_required_sections():
    text = PAPER.read_text(encoding="utf-8")
    headings = set(re.findall(r"^# (.+)$", text, flags=re.MULTILINE))
    required = {
        "Summary",
        "Statement of need",
        "State of the field",
        "Software design",
        "Research impact statement",
        "AI usage disclosure",
        "Acknowledgements",
        "References",
    }
    assert required <= headings


def test_all_paper_citations_resolve_in_bib():
    text = PAPER.read_text(encoding="utf-8")
    bib = BIB.read_text(encoding="utf-8")
    cited = set()
    for group in re.findall(r"\[@([^\]]+)\]", text):
        for key in group.split(";"):
            cited.add(key.strip().lstrip("@"))
    entries = set(re.findall(r"^@\w+\{([^,]+),", bib, flags=re.MULTILINE))
    assert cited
    assert cited <= entries


def test_paper_length_is_in_joss_range():
    text = PAPER.read_text(encoding="utf-8")
    body = text.split("---", 2)[2]
    words = len(re.findall(r"\b[\w'-]+\b", body))
    assert 250 <= words <= 1000


def test_repository_has_joss_community_files():
    required = [
        "LICENSE",
        "CONTRIBUTING.md",
        "SUPPORT.md",
        "CODE_OF_CONDUCT.md",
        "CITATION.cff",
        ".zenodo.json",
        "paper/paper.md",
        "paper/paper.bib",
        ".github/workflows/test.yml",
        ".github/workflows/joss-paper.yml",
    ]
    missing = [name for name in required if not (ROOT / name).is_file()]
    assert not missing, f"missing JOSS files: {missing}"


def test_readme_lists_every_registered_command():
    text = README.read_text(encoding="utf-8")
    mentioned = set(re.findall(r"`biopytools\s+([a-zA-Z0-9_-]+)", text))
    registered = {command for _, command, _ in COMMAND_REGISTRY}
    assert registered == mentioned


def test_every_registered_command_has_working_help():
    failures = []
    runner = CliRunner()
    for _, command, _ in COMMAND_REGISTRY:
        result = runner.invoke(cli, [command, "-h"])
        if result.exit_code != 0:
            failures.append(
                (command, result.exit_code, repr(result.exception))
            )
    assert not failures, f"CLI help failures: {failures}"


def test_local_markdown_links_resolve():
    broken = []
    for md in [README, *sorted((ROOT / "docs").rglob("*.md"))]:
        text = md.read_text(encoding="utf-8")
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", text):
            target = target.strip().split("#", 1)[0].strip("<>")
            if (
                not target
                or "://" in target
                or target.startswith("mailto:")
                or target.startswith("{")
                or target.startswith("[^")
            ):
                continue
            if not (md.parent / target).resolve().exists():
                broken.append((str(md.relative_to(ROOT)), target))
    assert not broken, f"broken local Markdown links: {broken}"


def test_release_versions_are_consistent():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    versions = {
        pyproject["project"]["version"],
        str(citation["version"]),
        str(zenodo["version"]),
    }
    assert len(versions) == 1, f"inconsistent release versions: {versions}"


def test_package_version_is_not_placeholder():
    assert biopytools.__version__ != "0.0.0"
