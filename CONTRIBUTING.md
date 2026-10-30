# Contributing to BioPyTools

Contributions are welcome through GitHub issues and pull requests.

## Before opening an issue

1. Search existing issues.
2. Confirm that the problem reproduces with the latest release.
3. Include the command, exit code, relevant log lines, operating system,
   Python version, and BioPyTools version.
4. Remove private data paths, credentials, and unpublished biological data.

Use the issue templates in `.github/ISSUE_TEMPLATE/` when possible.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev,test]"
python -m pytest
```

The test suite uses mocks for external executables. End-to-end tests that invoke
real bioinformatics programs should be run on a compute node, not a login node.

## Pull requests

1. Fork the repository and create a focused branch.
2. Keep one logical change per pull request.
3. Add or update tests for changed behaviour.
4. Update the relevant `docs/<module>.md` page and `CHANGELOG.md`.
5. Run `python -m pytest` before opening the pull request.
6. Describe the problem, the approach, and the validation performed.

Please do not include credentials, private data, generated build products, or
large sequencing files.

## Code conventions

- Use snake_case for functions and variables, and CamelCase for classes.
- Validate paths and parameters before invoking external programs.
- Log the complete command before executing it.
- Preserve checkpointing and resumability in multi-step workflows.
- Add English and Chinese messages using the `中文|English` convention.
- Do not add user-specific absolute paths to distributable code.

## AI-assisted contributions

Contributors may use generative AI tools, but must review and validate all
outputs and disclose material AI assistance in the pull request description.

## Project maintainers

The project is currently maintained by Xiang Li. By participating, you agree to
follow `CODE_OF_CONDUCT.md`.
