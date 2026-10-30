# JOSS submission checklist

Status checked against the official JOSS submission instructions and review
criteria on 2026-10-30.

## Submission-form values

- Title: `BioPyTools: A modular command-line toolkit for reproducible
  bioinformatics workflows`
- Repository: `https://github.com/lixiang117423/biopytools`
- Git branch: leave blank (`paper/paper.md` is on the default `main` branch)
- Paper: `paper/paper.md`
- Version: `v1.75.0`
- Submitting author GitHub handle: `@lixiang117423`
- ORCID: `0000-0002-8400-7116`
- Type of submission: `New submission`
- Main subject: `Bioinformatics and computational biology`
- Suggested tags: bioinformatics, Python, command-line interface, reproducible
  research, high-performance computing, plant genomics

### Message to editors

Confirm these facts before pasting:

> This is a new JOSS submission and not a resubmission or a second JOSS paper.
> To the best of my knowledge, no portion of the submitted software,
> documentation, or manuscript has been published or submitted to another
> peer-reviewed venue, and no related peer-reviewed manuscript is currently
> planned. BioPyTools is a modular Python command-line toolkit that provides a
> shared, reproducible execution layer for bioinformatics analyses on local and
> HPC systems; its scientific algorithms are delegated to established upstream
> tools. No dedicated funding supported this work. The corresponding author
> declares no financial or personal competing interests.

## Completed in the repository

- [x] OSI-approved MIT `LICENSE` with an identified copyright holder.
- [x] Public GitHub repository, issue tracker, and source browsing without login.
- [x] More than six months of public development history (since 2026-01-02).
- [x] More than 500 commits and more than 150 version tags.
- [x] `v1.75.0` tagged, published as a GitHub Release, and published on PyPI.
- [x] PyPI package installation and a lazy unified command-line interface.
- [x] Public automated test suite with continuous integration.
- [x] All 220 registered commands return `-h` successfully from a clean wheel
  installation.
- [x] README module list matches all 220 registered commands; all local Markdown
  links resolve.
- [x] MkDocs documentation builds successfully with `--strict`.
- [x] `CONTRIBUTING.md`, `SUPPORT.md`, `CODE_OF_CONDUCT.md`, issue templates, and
  pull-request template.
- [x] `paper/paper.md` with all required JOSS sections.
- [x] `paper/paper.bib` with DOI-verified references.
- [x] Official Open Journals draft-PDF GitHub workflow.
- [x] Citation metadata in `CITATION.cff`.
- [x] Archive metadata in `.zenodo.json`.

## Must confirm before pressing submit

- [x] Confirm that the author list and order are complete. Add research-group
  collaborators who made substantial scientific or software contributions.
- [ ] Confirm the corresponding author's current affiliation and ORCID.
  Repository metadata currently uses Yazhouwan National Laboratory and ORCID
  0000-0002-8400-7116.
- [ ] Confirm every AI tool and model version in the disclosure. In particular,
  record the exact Anthropic model/version used through Claude Code if available.
- [x] Confirm the funding statement, grant numbers, and sponsor involvement, or
  state explicitly that no specific funding supported the work.
- [ ] Respond to open user issues before submission so the public support record
  reflects the current maintenance practice.
- [ ] Re-run the complete test suite and package build after the final commit:
  `python -m pytest`, `python -m build`, and `python -m twine check dist/*`.
- [ ] Let at least one colleague who did not develop the package install it from
  a clean environment and run the documented quick-start command.
- [ ] Archive the existing `v1.75.0` release on Zenodo or another accepted
  repository, obtain the software DOI, and add the DOI to the paper or repository
  metadata as required by JOSS.
- [x] Download the `paper.pdf` artifact from the latest JOSS workflow run and
  inspect the title, author block, sections, references, and page layout.
- [ ] Confirm that no data, credentials, private paths, or unpublished sensitive
  information appears in the paper, logs, or repository.

## Editorial risk assessment

- Scope is broad and integration-oriented. The paper therefore emphasizes the
  shared execution contract and explicitly delegates scientific algorithms to
  upstream tools.
- Research impact is strongest in PyPI downloads, sustained development, and
  documented internal research use. Additional publications or external user
  testimonials would materially strengthen the submission.
- Claude Code and GLM are already disclosed in the project README; the manuscript
  disclosure must remain consistent with the final development history.
