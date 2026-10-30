---
title: 'BioPyTools: A modular command-line toolkit for reproducible bioinformatics workflows'
tags:
  - bioinformatics
  - Python
  - command-line interface
  - reproducible research
  - high-performance computing
  - plant genomics
authors:
  - name: Xiang Li
    orcid: 0000-0002-8400-7116
    corresponding: true
    affiliation: 1
affiliations:
  - name: Yazhouwan National Laboratory, Sanya, Hainan, China
    index: 1
date: 30 October 2026
bibliography: paper.bib
---

# Summary

`BioPyTools` is a Python command-line toolkit for computational biology. It
provides 220 modular commands spanning sequencing-data quality control, genome
assembly and annotation, variant analysis, pangenomics, phylogenetics,
population genetics, microbiome analysis, and effector prediction. Each command
provides a stable interface for one analysis task while delegating the underlying
algorithm to an established scientific program. `BioPyTools` adds the operational
layer needed to use those programs consistently: input validation, structured
output directories, checkpointing, bilingual logging, software-version manifests,
and portable invocation of tools installed in isolated Conda environments. The
toolkit runs as a conventional `pip`-installable package and is designed for
both local workstations and scheduler-managed high-performance computing (HPC)
systems.

Source code, documentation, and the issue tracker are available at
<https://github.com/lixiang117423/biopytools>.

# Statement of need

Bioinformatics analyses commonly combine dozens of command-line programs with
different argument conventions, environment requirements, output layouts, and
failure modes. Even when individual algorithms are mature, analysts must write
custom shell or Python glue to validate inputs, record versions, resume partial
runs, and keep intermediate files separate from deliverables. This glue is costly
to maintain and is a frequent source of irreproducible analyses.

Existing scientific libraries and workflow engines provide important components,
but they address different levels of the problem. General-purpose libraries such
as Biopython [@cock2009], SeqKit [@shen2016], SAMtools [@li2009], and BCFtools
[@danecek2021] expose reusable primitives, while Snakemake [@koster2012],
Nextflow [@ditommaso2017], and nf-core [@ewels2020] provide workflow
orchestration. Users still need to encode domain-specific task semantics and a
consistent runtime contract for heterogeneous compute environments. `BioPyTools`
targets this gap for bioinformatics analysts, small research groups, and HPC
support teams that need reproducible, task-oriented commands without maintaining
a separate workflow implementation for every analysis.

# State of the field

The software landscape around `BioPyTools` includes three complementary classes
of tools. First, domain libraries and utilities such as Biopython, SeqKit,
SAMtools, and BCFtools provide well-tested algorithms and file-format support.
Second, workflow systems such as Snakemake, Nextflow, and nf-core provide
execution graphs, container integration, and portability across compute
platforms. Third, Galaxy [@afgan2018] provides a graphical environment for
accessible and reproducible analyses. Conda-based distribution through Bioconda
[@gruning2018] further improves installation reproducibility for the underlying
tools.

`BioPyTools` does not replace these projects. It reuses their algorithms and
contributes a task-oriented command layer with uniform validation, logging,
checkpointing, output organization, and tool-discovery behavior. Building this
layer rather than adding every requirement to one upstream project was a
deliberate choice: the required behaviours span tools maintained by many
independent communities, while the shared execution contract is specific to
analysts working across multiple organisms and HPC environments. The design
therefore keeps scientific algorithms upstream and concentrates maintenance on
the integration layer.

# Software design

Each analysis module follows the same internal contract: a typed configuration
object validates parameters and expands paths, a logger separates routine output
from warnings and errors, a calculator implements the analysis, and a command
entry point translates command-line arguments into the core API. This separation
allows modules to be tested without invoking a shell and allows the same logic to
be used from scripts.

A shared runtime layer handles concerns that are otherwise repeated across
modules. External executables are resolved explicitly, and Conda environments are
invoked through the Conda installation that owns the environment prefix rather
than through whichever `conda` happens to be on `PATH`. Commands use
`conda run -p <prefix> --no-capture-output` to avoid output buffering and
memory-related failures, and pipeline construction avoids nested `conda run`
calls. Modules write a common `00_pipeline_info/` directory containing software
versions and run metadata, place intermediates under `output_dir/tmp`, and use
numbered stage directories followed by `99_logs/`. Long-running analyses check
for completed outputs before recomputing them, so interrupted jobs can resume
without manual bookkeeping.

The command registry is loaded lazily: listing help does not import every
scientific dependency, while invoking a command imports only its module. The
repository currently contains more than 800 automated unit tests covering
configuration validation, command construction, parsers, checkpointing, and
module wiring.

# Research impact statement

`BioPyTools` has been developed publicly since January 2026 and has more than 500
commits and 150 tagged releases. It is used in the author's ongoing research on
soybean and wild soybean genome assembly and population analyses, *Phytophthora
sojae* effector and NLR annotation, and plant microbiome workflows. These
applications have driven iterative changes to checkpointing, environment
resolution, output validation, and failure reporting.

Adoption extends beyond the author's local environment. In the month ending
30 October 2026, the package recorded 2,912 non-mirror downloads from the Python
Package Index (<https://pypistats.org/packages/biopytools>). The public issue
tracker has also received usage questions from an independent user, and
the project publishes installation instructions, module-level documentation,
and a public changelog. These signals demonstrate active use and provide
concrete pathways for external users to report problems and request features.

# AI usage disclosure

Generative AI was used during development. Anthropic Claude Code (exact model
version not recorded in the repository), GLM 5.2, and OpenAI Codex based on GPT-5
assisted with code generation and refactoring, test scaffolding, documentation,
release housekeeping, and drafting and editing this manuscript. AI tools were not
used for interactions with JOSS editors or reviewers. The human author reviewed,
edited, executed, and validated all AI-assisted outputs; selected the scientific
scope, architecture, algorithms, and external tools; and remains responsible for
the software and manuscript.

# Acknowledgements

The author thanks the maintainers of the open-source programs integrated by
`BioPyTools`, the Bioconda community, and users who have reported issues or
suggested improvements.

# References
