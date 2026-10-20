"""genomeforge_cli — the `genomeforge` command-line front door.

A thin, stdlib-only wrapper around the pipeline's entry points (Nextflow `main.nf`,
`bin/make_samplesheet.py`, `bin/fetch_resources.sh`, `bin/summarize_qc.py`). It launches,
checks and inspects runs; the tools themselves run inside the CallForge containers that
Nextflow starts per task.
"""
__version__ = "0.1.0"
