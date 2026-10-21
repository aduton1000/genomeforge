# What GenomeForge produces

This page describes the **shape** of a run's outputs, not any specific cohort. Every run
writes to `<rundir>/results/` (override with `--outdir`). Four deliverables matter to most
users; the rest are intermediates kept for audit and reuse.

## 1. Analysis-ready alignments — `02_alignments/`

- **`<sample>.bam` + `.bai`** — coordinate-sorted, duplicate-marked, base-quality-recalibrated,
  one file per sample regardless of how many lanes it was sequenced on. Read groups carry
  `ID=<flowcell>.L<lane>`, `SM`, `LB`, `PU`, `PL` (and `CN` when `--sequencing_center` is set).
- **`<sample>.cram` + `.crai`** — written alongside when `--archive_cram true`.

## 2. The cohort callset — `05_joint_calling/`

| file | content |
|:-----|:--------|
| `cohort.raw.vcf.gz` | GenotypeGVCFs output over all shards, unfiltered, dbSNP IDs attached |
| `cohort.filtered.vcf.gz` | VQSR-annotated (or hard-filtered) multi-sample VCF; `FILTER` carries the tranche or filter name |
| `cohort.filtered.norm.vcf.gz` | the same after `bcftools norm`: multi-allelics split, indels left-aligned |
| `cohort.pass.vcf.gz` | PASS-only subset of the normalised VCF — the input to annotation and cohort QC |
| `vqsr/SNP.tranches`, `vqsr/INDEL.tranches`, `*.recal` | VQSR model outputs for inspection |

## 3. The annotated callset — `06_annotation/`

- **`cohort.annotated.vcf.gz`** — the PASS VCF with dbSNP `ID`, gnomAD INFO fields
  (`gnomAD_AF`, `gnomAD_AF_afr`, `gnomAD_AF_grpmax`, `gnomAD_nhomalt` by default), ClinVar
  (`CLNSIG`, `CLNDN`, `CLNREVSTAT`, `CLNSIGCONF`, `CLNVC`) and the VEP `CSQ` string
  (symbol, consequence, impact, HGVS, MANE, SIFT/PolyPhen, allele frequencies, pick flag).
- **`variants.flat.tsv.gz`** — one row per variant with the worst-consequence transcript,
  the gnomAD and ClinVar columns, and one genotype column per sample. Loads directly into
  R/pandas for filtering and burden work.

## 4. Reports — `08_reports/`

- **`qc_report.pdf`** — the review document: title page with cohort KPIs and thresholds, the
  per-sample QC table with flagged rows highlighted, then every figure on its own page with a
  caption. Print it or attach it to a results summary.
- **`qc_dashboard.html`** — the same content as a self-contained web page (figures embedded),
  with a navigation bar and the sample table.
- **`figures/`** — every figure as PNG and SVG: coverage per sample, cumulative coverage,
  coverage by chromosome (sex and aneuploidy at a glance), mapping/pairing/duplication, insert
  size, GC content, lane yield and Q30, base quality by cycle, contamination, sex check,
  relatedness histogram and matrix, PCA with scree, per-sample variant metrics, site statistics
  (indel length, substitution types, QUAL, depth), VQSR tranches, and the annotation summary
  (site frequency spectrum, consequences, cohort vs gnomAD AF, ClinVar).
- **`cohort_qc_summary.tsv` / `.md`** — the per-sample QC table (see section 5).
- **`multiqc/cohort_multiqc.html`** — the aggregated tool-level report.

## 5. Cohort QC — `07_cohort_qc/`

- **`08_reports/cohort_qc_summary.tsv`** — one row per sample: reads, mapping and proper-pair rates,
  duplicate rate, mean/median coverage and fraction ≥ 10/20/30×, X and Y coverage ratios,
  sex from the sheet / from coverage / from somalier, verifyBamID2 `FREEMIX`, non-reference
  SNV count, het/hom and Ti/Tv, related samples, `qc_flags` and `qc_status` (PASS or FLAG).
  **Nothing is dropped**: flags are decisions for the analyst. `cohort_qc_summary.md` is the
  human-readable version with cohort medians and the flagged list.
- **`somalier/`** — `cohort.samples.tsv` (inferred sex, X/Y metrics), `cohort.pairs.tsv`
  (pairwise relatedness; ≥ 0.2 is reported in the summary) and the interactive `cohort.html`.
- **`pca/cohort_pca.eigenvec` / `.eigenval`** — plink2 PCs on LD-pruned common autosomal SNPs.
- **`cohort.bcftools_stats.txt`** — per-sample and per-site variant statistics (input to the
  figures and MultiQC).

## Intermediates kept for reuse

- **`04_gvcfs/<sample>.g.vcf.gz`** — whole-genome gVCF per sample; add future samples to a
  joint call without realigning.
- **`03_sample_qc/`** — MarkDuplicates metrics, BQSR tables, samtools stats/flagstat/idxstats,
  mosdepth summaries and 1-kb windows, CollectWgsMetrics, verifyBamID2 `selfSM`.
- **`01_trim_reports/`** — per-lane trimming reports and bwa-mem2 logs.
- **`pipeline_info/`** — `resource_check.txt` (reference/annotation contig-style check),
  `intervals/` (the calling regions and scatter shards), `clinvar.chr.vcf.gz` (renamed
  ClinVar), `versions.txt`, `params.json`, and Nextflow's `report`, `timeline`, `trace` and `dag`.

## Reading the QC summary

| flag | meaning | typical action |
|:-----|:--------|:---------------|
| `low_depth(x)` | mean coverage below `--qc_min_depth` | check lane completeness; consider top-up sequencing |
| `high_dup(x)` | duplicate rate above `--qc_max_dup` | library complexity problem; effective depth is lower than nominal |
| `contamination(x)` | `FREEMIX` above `--qc_max_contam` | exclude or treat genotypes with caution |
| `low_mapping(x)` | mapped fraction below `--qc_min_mapped` | wrong species/reference or heavy adapter content |
| `sex_mismatch(...)` | sheet sex ≠ inferred sex | sample swap or metadata error; verify before association |
| `sex_ambiguous` | X/Y coverage ratios between the male and female bands | possible contamination, aneuploidy, or mixed sample |
| `related:S(r)` | relatedness ≥ 0.2 with sample S | duplicates or relatives; account for in association models |
| `titv(x)` | Ti/Tv outside 1.9–2.2 for whole-genome PASS SNVs | excess false positives; revisit filtering |
