# GenomeForge

**Whole-genome germline joint calling, FASTQ to an annotated and QC-ed cohort callset.**
A Nextflow DSL2 pipeline for SLURM clusters that runs on the
[CallForge](https://github.com/aduton1000/callforge) container images and reference bundle.
CallForge is the targeted-panel pipeline; GenomeForge is its whole-genome sibling.

```
lanes ──fastp|bwa-mem2|sort──▶ lane BAM ──MarkDuplicates──▶ BQSR ──▶ sample BAM (+CRAM)
sample BAM ──HaplotypeCaller × N shards──▶ gVCF shards ──GatherVcfs──▶ sample gVCF
shards × samples ──GenomicsDBImport──▶ GenotypeGVCFs ──concat──▶ raw cohort VCF
raw ──VQSR (SNP+INDEL) or hard filters──▶ norm/split ──▶ PASS VCF
PASS ──per chromosome: dbSNP IDs, gnomAD AF, ClinVar, VEP──▶ annotated VCF + flat TSV
BAMs/VCF ──samtools, mosdepth, CollectWgsMetrics, verifyBamID2, somalier, plink2 PCA──▶ QC summary + MultiQC
```

Any Illumina lane layout is accepted: samples may span several lanes, flowcells, libraries
or deliveries; each lane-pair is aligned with its own read group and merged per sample.

## Layout

| path | purpose |
|---|---|
| `main.nf`, `modules/*.nf` | workflow and processes (every process has a `stub:` for dry runs) |
| `nextflow.config`, `conf/base.config` | parameters and per-label resources (site-neutral) |
| `conf/noguchi.config` | site profile: Noguchi HPC paths, SLURM partition, CallForge images |
| `conf/slurm.config` | generic SLURM + Singularity/Apptainer profile; pair with `params.example.yaml` |
| `bin/make_samplesheet.py` | lane-level sample sheet from FASTQ directories (Novogene, bcl2fastq or generic names) |
| `bin/fetch_resources.sh` | downloads hapmap/omni/1000G (VQSR) and the GATK WGS calling regions |
| `bin/summarize_qc.py` | per-sample QC table with flags (depth, duplicates, contamination, sex, relatedness, Ti/Tv) |
| `run_genomeforge.sh` | submits the Nextflow head job under sbatch with work/tmp on disk |
| `test/make_stub_env.sh` | fake data layout for `-stub-run` DAG tests |

## Requirements

Nextflow ≥ 24.10 (Java 17+), SLURM, Singularity or Apptainer, and the four CallForge images
(`callforge.sif`, `callforge-annotate.sif`, `callforge-cohortqc.sif`, `callforge-burden.sif`;
built from the CallForge repo's `env/` definitions). Reference bundle: a **chr-style no-alt
GRCh38** with bwa-mem2 index, `.fai` and `.dict`; dbSNP, Mills and known-indels VCFs; hapmap,
omni and 1000G phase-1 SNPs for VQSR; ClinVar; per-chromosome gnomAD sites files; a somalier
sites VCF; a VEP cache. All VCFs must share the reference's contig naming (ClinVar is renamed
automatically).

## Quick start (Noguchi HPC)

```bash
# one-off: stage the VQSR training sets + calling regions next to the other annotation files
bash /hpc/opt/genomeforge/bin/fetch_resources.sh /hpc/opt/callforge/refs/annotation_db

# sample sheet from the delivery directories (read-only; fails loudly on unpaired/truncated files)
python3 /hpc/opt/genomeforge/bin/make_samplesheet.py \
    --roots /path/to/delivery-part1 /path/to/delivery-part2 \
    --metadata phenotypes.csv --out samplesheet.csv
#   --exclude S1,S2            leave samples out      --sample-from dir|filename   override id source
#   --fastq-csv pairs.csv      start from an existing sample_id,fastq_1,fastq_2 list instead of scanning

# dry run of the DAG, then the real run (resumable: re-run the same command after a failure)
export GENOMEFORGE_RUNDIR=/hpc/data/<lab>/<project>/genomeforge_run
bash /hpc/opt/genomeforge/run_genomeforge.sh --input samplesheet.csv -stub-run
bash /hpc/opt/genomeforge/run_genomeforge.sh --input samplesheet.csv
```

`phenotypes.csv` is optional: `sample_id` plus `sex` (M/F) and `phenotype` if known; sex drives
the sex-concordance check, nothing else uses it.

On another cluster: copy `params.example.yaml`, fill in the paths, and run with
`GENOMEFORGE_PROFILE=slurm GENOMEFORGE_PARAMS=site.params.yaml`, or call `nextflow run` directly.

## Outputs (`results/`)

| directory | content |
|---|---|
| `01_fastp/` | per-lane trimming reports |
| `02_bams/` | `<sample>.bam` + `.bai` (BQSR-recalibrated; CRAM with `--archive_cram true`) |
| `03_qc/` | markdup metrics, BQSR tables, samtools, mosdepth, CollectWgsMetrics, verifyBamID2 |
| `04_gvcf/` | per-sample whole-genome gVCF (reusable for future joint calls) |
| `05_joint/` | `cohort.raw.vcf.gz`, `cohort.filtered.vcf.gz`, `cohort.filtered.norm.vcf.gz`, `cohort.pass.vcf.gz`, `vqsr/` tranches |
| `06_annotated/` | `cohort.annotated.vcf.gz` (dbSNP, gnomAD, ClinVar, VEP), `variants.flat.tsv.gz` |
| `07_cohort/` | `cohort_qc_summary.tsv/.md`, somalier relatedness + sex, plink2 PCA, bcftools stats, `cohort_multiqc.html` |
| `pipeline_info/` | resource check, scatter intervals, timeline/report/trace/dag, `versions.txt`, `params.json` |

## Key parameters

| parameter | default | notes |
|---|---|---|
| `--scatter_count` | 48 | HaplotypeCaller / joint-genotyping shards |
| `--calling_regions` | site profile | GATK `wgs_calling_regions.hg38.interval_list`; falls back to whole primary chromosomes |
| `--filtering` | `vqsr` | `hard` for GATK hard filters; automatic fallback when training sets are missing |
| `--snp_tranche` / `--indel_tranche` | 99.7 / 99.7 | VQSR truth-sensitivity levels |
| `--trim` | true | fastp streamed into bwa-mem2 (no trimmed FASTQs written) |
| `--sequencing_center` | null | read-group `CN` tag |
| `--gnomad_pattern` / `--gnomad_fields` | exomes v4.1 / `AF,AF_afr,AF_grpmax,nhomalt` | per-chromosome gnomAD file name template and INFO fields copied |
| `--run_vep`, `--run_pca`, `--run_contamination` | true | switch stages off |
| `--archive_cram` | false | CRAM next to each final BAM |
| `--qc_min_depth`, `--qc_max_dup`, `--qc_max_contam`, `--qc_min_mapped` | 20 / 0.25 / 0.03 / 0.95 | QC flags (nothing is dropped) |

## Design notes

- **Read groups**: `ID=<flowcell>.L<lane>`, `SM=<sample>`, `LB=<library>`, `PU=<flowcell>.<lane>`,
  `PL=ILLUMINA`; lanes are aligned separately and merged inside MarkDuplicates.
- **Scatter**: `SplitIntervals` over the calling regions; GenomicsDBImport merges input intervals.
- **Filtering**: VQSR with hapmap/omni/1000G/dbsnp (SNP) and Mills/dbsnp (INDEL); hard filters
  otherwise. Multi-allelics are split and left-aligned before annotation.
- **Scratch**: every task uses `./tmp` inside its work directory; node `/tmp` is never used.
- **Resources**: `conf/base.config`; memory and time scale with the retry attempt.

## First-run checklist on a new site

1. `-stub-run` to confirm paths, containers and scheduler submission.
2. A real run on two or three samples with `--scatter_count 8` (about a day): exercises every
   tool once, including the gnomAD tag rename (`bcftools annotate -c INFO/gnomAD_AF:=AF`) and
   the verifyBamID2 resource auto-detection inside the cohortqc image.
3. The full cohort with `-resume`.

## Updating the CallForge checkout on the Noguchi server (one-off housekeeping)

The CallForge `main` history on GitHub was rewritten on 2026-10-20 (commit-message cleanup;
file contents unchanged, new head `3f27e63`). A server checkout at `1781846` cannot
fast-forward: reset it instead of pulling. It also carries an untracked `.tools/` directory
and stray gnomAD index files in the repo root.

```bash
cd /hpc/opt/callforge/repo
git status --short                                  # expect only "?? .tools/"
git fetch origin
git reset --hard origin/main                        # -> 3f27e63; .tools/ is untracked and survives
rm -f gnomad.exomes.v4.1.sites.chr*.vcf.bgz.tbi     # stray, already git-ignored
grep -qx '.tools/' .gitignore || printf '\n# local tool binaries (nextflow etc.)\n.tools/\n' >> .gitignore
git add .gitignore
git commit -m "gitignore: ignore local .tools/ binaries directory"
git push origin main
git status --short                                  # should print nothing
callforge --version                                 # 0.1.0
```

## License

MIT (same as CallForge).
