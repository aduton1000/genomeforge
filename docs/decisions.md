# GenomeForge — decisions, caveats, and revisit-on-real-data items

Running log of design decisions and the reasoning behind them. Items marked *revisit* are
to be checked after the first full cohort run.

## Confirmed decisions

- **Reference**: a **chr-style no-alt GRCh38** (primary assembly, 194 contigs, UCSC names).
  Every public resource used downstream (GATK bundle known sites and VQSR sets, dbSNP,
  gnomAD, ClinVar after renaming, somalier sites, VEP cache) is chr-named, so a chr-named
  reference avoids per-file contig renaming. No alt contigs because bwa-mem2 is not run
  alt-aware and the GATK tools downstream assume a primary assembly; `CHECK_RESOURCES`
  fails loud on alt/decoy contigs or a contig-style mismatch.
- **Scope is germline SNV/indel.** No CNV/SV, STR or somatic calling in this release; the
  final BAMs and gVCFs are the inputs those modules would consume later.
- **Per-lane alignment, merge in MarkDuplicates.** Lanes are aligned independently (32
  threads each) so a 2–6-lane sample scales across nodes; MarkDuplicates merges them and
  marks optical duplicates with the NovaSeq-appropriate 2500-pixel distance. `samtools
  merge` was dropped as an extra pass.
- **fastp streamed into bwa-mem2.** Trimming is on by default (`--trim`) but no trimmed
  FASTQs are written; fastp's interleaved stdout feeds `bwa-mem2 mem -p`. `-K 100000000`
  and `-Y` give deterministic output and keep supplementary alignments soft-clipped.
- **BQSR unscattered.** BaseRecalibrator/ApplyBQSR on a 30× genome take a few hours each
  and parallelise across samples; scattering would add a gather step for little wall-clock
  gain on a four-node cluster.
- **HaplotypeCaller scattered over the GATK WGS calling regions** (48 shards by default,
  `BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW`); whole-chromosome fallback uses
  `INTERVAL_SUBDIVISION` so shards stay balanced. GenomicsDBImport and GenotypeGVCFs run on
  the same shards with `--merge-input-intervals`; `--only-output-calls-starting-in-intervals`
  prevents duplicate records at shard boundaries.
- **VQSR by default, hard filters as fallback.** Cohorts of several dozen whole genomes are
  large enough for VQSR; the SNP model uses hapmap/omni/1000G/dbsnp with
  `QD MQRankSum ReadPosRankSum FS MQ SOR DP`, the INDEL model Mills/dbsnp with the same
  minus `MQ`; both apply at 99.7 % truth sensitivity. If the training sets are absent the
  pipeline logs a warning and applies the GATK hard-filter expressions instead.
- **Normalise before annotation.** `bcftools norm -m -any` splits multi-allelics so gnomAD
  and ClinVar annotation (which are per-allele) match exactly; left-alignment makes indels
  comparable across callers.
- **Annotation scattered by chromosome** with `bcftools annotate` for dbSNP/gnomAD/ClinVar
  and VEP (`--fork`) per chromosome; gnomAD is consumed from its native per-chromosome
  files, so no 180-GB concatenation is needed.
- **gnomAD exomes v4.1 by default** (what the CallForge bundle ships). Fine for coding
  variants; non-coding allele frequencies need gnomAD genomes — switch with
  `--gnomad_pattern`.
- **Cohort QC: somalier + verifyBamID2 + plink2.** somalier gives relatedness and sex from
  a sites VCF with no pedigree; verifyBamID2 estimates contamination from its bundled 1000G
  SVD resource; plink2 PCA on LD-pruned common autosomal SNPs gives ancestry/batch axes.
  peddy is not used (same signal, heavier dependency).
- **QC flags, never drops.** Thresholds are parameters; the summary marks samples FLAG
  with explicit reasons and leaves the decision to the analyst.
- **Scratch on disk.** Each task uses `./tmp` under its work directory and sets
  `-Djava.io.tmpdir`; `process.scratch` is off. Clusters whose `/tmp` is RAM-backed would
  otherwise sort BAMs in memory.
- **Containers, not conda.** All tools come from the CallForge images; nothing is
  installed on the host.

## Revisit on real data

- **VQSR tranche levels** (99.7/99.7): inspect `vqsr/*.tranches` and Ti/Tv per tranche on
  the first cohort; 99.5 for indels is common.
- **`DP` in the VQSR annotations** assumes uniform whole-genome depth; drop it if the
  cohort mixes very different depths.
- **gnomAD tag rename** (`bcftools annotate -c INFO/gnomAD_AF:=AF`) and **verifyBamID2
  resource auto-detection** are validated on the first real run (see the first-run
  checklist in `docs/hpc_deployment.md`).
- **MarkDuplicates memory** on 60×+ samples: raise the `markdup` label if retries occur.
- **Sex inference thresholds** (Y/autosome > 0.15 → M, < 0.05 → F) are conservative;
  somalier's call is reported alongside and disagreements are flagged.
- **CollectWgsMetrics** runs with `--USE_FAST_ALGORITHM`; compare its mean coverage with
  mosdepth once and keep only one if they agree.
- **Future modules**: SV/CNV (Manta/smoove, CNVkit in `wgs` mode), STR (ExpansionHunter),
  DeepVariant + GLnexus as an alternative caller, GIAB benchmarking when a control is
  sequenced.
