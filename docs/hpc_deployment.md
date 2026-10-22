# GenomeForge — HPC deployment (SLURM + Apptainer/Singularity)

GenomeForge is a cluster pipeline: a 30× whole genome needs tens of CPU-hours for
alignment and tens more for calling, and a cohort runs for days. This page covers a central
install, the container images, the reference bundle, site profiles, running, monitoring and
troubleshooting. All `/hpc/...` and `/opt/...` paths are **examples**; substitute your layout.

## 0. One-time: central install

Recommended layout (the same as CallForge): `<root>/repo` (the checkout), `<root>/cli-env`
(a virtual environment holding only the CLI), `<root>/bin/genomeforge` (a wrapper, so the
venv's own `python`/`pip` never shadow the system ones), and `<root>/genomeforge-env.sh`.

```bash
ROOT=/opt/genomeforge
git clone https://github.com/aduton1000/genomeforge.git $ROOT/repo
python3 -m venv $ROOT/cli-env && $ROOT/cli-env/bin/pip install -q -e $ROOT/repo      # PEP 668-safe; no dependencies
mkdir -p $ROOT/bin
printf '#!/bin/bash\nexport GENOMEFORGE_HOME="${GENOMEFORGE_HOME:-%s/repo}"\nexec %s/cli-env/bin/genomeforge "$@"\n' "$ROOT" "$ROOT" > $ROOT/bin/genomeforge
chmod 755 $ROOT/bin/genomeforge
cp $ROOT/repo/share/genomeforge-env.sh.example $ROOT/genomeforge-env.sh      # edit: home, profile, params, run dir, PATH=$ROOT/bin
echo "source $ROOT/genomeforge-env.sh" >> ~/.bashrc && source $ROOT/genomeforge-env.sh
genomeforge --version
genomeforge test                      # compiles and dry-runs the DAG; needs Nextflow ≥ 24.10 and Java 17+
```

Updating later is `git -C $ROOT/repo pull --ff-only`; the editable install follows. Without
a venv, `$ROOT/repo/share/bin` on `PATH` gives the same command straight from the checkout.
No host conda environment is needed; Nextflow itself must be on `PATH` (or come from a site
module).

## 1. Container images

GenomeForge runs every task inside the **CallForge** images:

| image | tools used here |
|:------|:----------------|
| `callforge.sif` | bwa-mem2, samtools, bcftools, htslib, GATK 4.6, picard, fastp, mosdepth, multiqc, python (QC summary) |
| `callforge-annotate.sif` | ensembl-vep 112 (`--vep_container` can point at the official `ensemblorg/ensembl-vep` image instead) |
| `callforge-cohortqc.sif` | somalier, verifybamid2 |
| `callforge-burden.sif` | plink2 |

Build them once per site from the CallForge repository (`apptainer build callforge.sif
env/callforge.def`; `bin/build_stage_images.sh` for the rest) and point the site profile at
them. Verify with a non-login shell:

```bash
apptainer exec /opt/callforge/images/callforge.sif bwa-mem2 version
apptainer exec /opt/callforge/images/callforge.sif gatk --version
apptainer exec /opt/callforge/images/callforge-cohortqc.sif somalier --version
apptainer exec /opt/callforge/images/callforge-burden.sif plink2 --version
```

## 2. Reference bundle (stage once, shared read-only)

All VCFs must use the **same contig naming as the FASTA** (ClinVar is the one exception
and is renamed by the pipeline). `CHECK_RESOURCES` verifies this at the start of every run.

| item | example | notes |
|:-----|:--------|:------|
| genome | `GRCh38_no_alt.chr.fa` + `.fai` + `.dict` + bwa-mem2 index (`.0123 .amb .ann .bwt.2bit.64 .pac`) | chr-style, no alt/decoy contigs |
| calling regions | `wgs_calling_regions.hg38.interval_list` | GATK bundle; `bin/fetch_resources.sh` |
| known sites | `dbsnp.vcf.gz`, `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`, `Homo_sapiens_assembly38.known_indels.vcf.gz` | BQSR |
| VQSR | `hapmap_3.3.hg38.vcf.gz`, `1000G_omni2.5.hg38.vcf.gz`, `1000G_phase1.snps.high_confidence.hg38.vcf.gz` | `bin/fetch_resources.sh` |
| annotation | `clinvar.vcf.gz`, gnomAD per-chromosome sites (`gnomad.exomes.v4.1.sites.chrN.vcf.bgz` or genomes), VEP cache (`<cache>/homo_sapiens/112_GRCh38`) | |
| cohort QC | somalier `sites.hg38.vcf.gz` | verifyBamID2's SVD resource ships inside the cohortqc image |

```bash
bash bin/fetch_resources.sh /refs/GRCh38          # idempotent, resumable, ~2.5 GB
```

Building a bwa-mem2 index for GRCh38 needs ~72 GB RAM and about an hour:

```bash
apptainer exec /opt/callforge/images/callforge.sif bwa-mem2 index GRCh38_no_alt.chr.fa
apptainer exec /opt/callforge/images/callforge.sif samtools faidx GRCh38_no_alt.chr.fa
apptainer exec /opt/callforge/images/callforge.sif gatk CreateSequenceDictionary -R GRCh38_no_alt.chr.fa
```

## 3. Site profiles

- **`-profile noguchi`** — the Noguchi HPC: partition `global`, bind `/hpc`, CallForge images
  and references under `/hpc/opt/callforge`. Nothing else to configure.
- **`-profile slurm`** — generic SLURM + Singularity/Apptainer. Copy
  [`params.example.yaml`](../params.example.yaml), fill in the paths, partition, account and
  bind roots, and pass it with `-params-file`. Add a `conf/<site>.config` modelled on
  `conf/noguchi.config` if you prefer a named profile.
- **`-profile stub`** — local dry run of the DAG with `-stub-run`; `test/make_stub_env.sh`
  builds a fake data layout for it.

Container binds: Apptainer auto-mounts only `$HOME`, `/tmp` and the working directory.
Set `bind_paths` to the root(s) holding the references and the data (comma-separated).

## 4. Run on SLURM

```bash
genomeforge check --profile slurm --params-file /data/projects/cohort/site.params.yaml     # everything present?
genomeforge samplesheet --roots /data/cohort/delivery --metadata phenotypes.csv --out samplesheet.csv

export GENOMEFORGE_PROFILE=slurm GENOMEFORGE_PARAMS=/data/projects/cohort/site.params.yaml   # or GENOMEFORGE_PROFILE=noguchi
export GENOMEFORGE_RUNDIR=/data/projects/cohort/genomeforge_run
genomeforge run --input samplesheet.csv --stub        # seconds: paths, containers, scheduler
genomeforge run --input samplesheet.csv               # the real run
```

`genomeforge run` writes an sbatch head script under `<rundir>/logs/`, exports `NXF_WORK`,
`TMPDIR`, `NXF_TEMP` and `NXF_OFFLINE` into the run directory, raises the open-file limit for
GenomicsDB, logs to `<rundir>/nextflow.log`, and passes `-resume`, so repeating the command
after a failure continues from the last completed task. `--foreground` runs Nextflow in the
current shell instead (for `screen`/`tmux` sessions), `--print-cmd` shows the command without
submitting.

**Keeping nodes free.** `--exclude-nodes node4` (or `GENOMEFORGE_EXCLUDE_NODES`) adds
`--exclude` to the head job and to every task through `--slurm_exclude`, so a cohort run uses
the remaining nodes only.

**Self-resumption.** The head job is submitted with `--requeue` (SLURM restarts it after a node
failure; the queue itself survives a controller restart), refuses to start while another head
job for the same run directory is running, and on any exit without a completed run it sleeps
`--restart-delay` seconds (600) and resubmits itself, up to `--max-restarts` times (20), always
with `-resume`. A `logs/head_<stamp>.complete` marker ends the chain. Failed tasks retry with
scaled memory and time. On a cluster that loses power often, add `--whole_node_align true`
so that no task runs longer than about two hours. Nothing in a run needs the internet, so a network outage only matters
if the data or reference filesystem goes with it, in which case the same mechanism resumes
once it is back.

**Sizing.** Per 30× sample: alignment ≈ 60–100 CPU-h, HaplotypeCaller ≈ 30–60 CPU-h,
peak RAM 64 GB (bwa-mem2 at 32 threads) or 160 GB with `--whole_node_align true`, 72 GB for an
index build. Task durations: lane alignment 1–4 h (half with whole-node), MarkDuplicates 2–3 h,
BQSR shards 10–20 min, HaplotypeCaller shards ≈ 30–60 min, joint-genotyping shards ≈ 1 h. Storage: final BAMs
≈ 1.1× the compressed FASTQ, `work/` peak ≈ 4×, so plan 5× the FASTQ size. A 64-sample 30×
cohort on four 72-thread nodes takes 5–7 days of wall-clock.

## 5. Monitoring and resuming

```bash
genomeforge status --rundir $GENOMEFORGE_RUNDIR        # head job, task counts per process, recent log
tail -f $GENOMEFORGE_RUNDIR/nextflow.log
squeue -u $USER
open $GENOMEFORGE_RUNDIR/results/pipeline_info/report.html      # per-task resources, after completion
```

Nextflow retries failed tasks up to `--max_retries` (default 2) with memory and time scaled
by attempt. To resume after the head job itself was killed, repeat `genomeforge run`. To clean
intermediates after a successful run, `nextflow clean -f -k <run name>` or remove `work/`.

## 6. Validate-at-deployment checklist

1. `-stub-run` on the cluster: paths resolve, containers are found, tasks reach the queue.
2. A real run on **two or three samples** with `--scatter_count 8`: exercises every tool once
   in about a day. Check in particular the gnomAD tag rename in `ANNOTATE_CHROM`
   (`bcftools annotate -c INFO/gnomAD_AF:=AF`) and that `VERIFYBAMID` found its SVD resource
   (`03_qc/verifybamid/<sample>.verifybamid.log`).
3. Review `07_cohort/cohort_qc_summary.md` and `05_joint/vqsr/*.tranches`.
4. The full cohort with `-resume`.

## 7. Troubleshooting

| symptom | cause / fix |
|:--------|:------------|
| `Missing site parameters: …` at start | no site profile or params file; `genomeforge run --profile noguchi` or `--profile slurm --params-file site.params.yaml`; `genomeforge check` lists what is missing |
| `MISMATCH: … contig style` in `resource_check.txt` | a VCF is Ensembl-named while the FASTA is chr-named (or vice versa); rename with `bcftools annotate --rename-chrs` |
| `MISSING bwa-mem2 index part` | build the index next to the FASTA (section 2) |
| tasks fail with exit 137 / OOM | retries scale memory automatically; raise the label in `conf/base.config` if a third attempt is needed |
| `No such file` for a FASTQ inside a task | the data root is not bound into the container; add it to `bind_paths` |
| VEP: cache version mismatch | `--vep_release` must equal the cache directory (`homo_sapiens/<release>_GRCh38`) and the image's VEP version |
| GenomicsDBImport "too many open files" | the launcher raises `ulimit -n`; if the site caps it, lower `--batch-size` in `modules/joint.nf` |
| head job dies at logout | use `genomeforge run` (sbatch head job) rather than `nextflow run` interactively, or `--foreground` inside tmux |

## 8. Noguchi HPC notes

- Login node `node1` is also a compute node; `/tmp` is RAM-backed tmpfs — never point
  `NXF_WORK`, `TMPDIR` or `process.scratch` at it (the launcher and `conf/base.config` do
  not).
- CallForge lives at `/hpc/opt/callforge/repo` with images in `/hpc/opt/callforge/images`
  and references in `/hpc/opt/callforge/refs`; use `GRCh38_chr/GRCh38.full.chr.fa` (chr-style
  no-alt), not `GRCh38_noalt.clean.fa` (Ensembl-named).
- Keeping the CallForge checkout current: its `main` history was rewritten on 2026-10-20
  (commit-message cleanup, identical file contents), so a checkout from before that date
  must be reset rather than pulled; afterwards a plain `git pull --ff-only` works again.

  ```bash
  cd /hpc/opt/callforge/repo
  git fetch origin && git reset --hard origin/main   # untracked .tools/ survives
  rm -f gnomad.exomes.v4.1.sites.chr*.vcf.bgz.tbi   # stray, git-ignored
  callforge --version
  ```
