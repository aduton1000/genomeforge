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
- **BQSR scattered over whole-contig groups.** BaseRecalibrator and ApplyBQSR each run per
  contiguous group of whole contigs (24 groups by default, balanced by length over the entire
  reference in `.dict` order), then `GatherBQSRReports` and `GatherBamFiles` reassemble the
  table and the BAM. Whole contigs mean no read is split between shards; the last group adds
  `-L unmapped` so unmapped reads survive into the final BAM. Originally unscattered (simpler,
  a few hours per sample), changed after the first cohort cluster showed daily power cuts:
  a 10–20 minute shard is lost to an outage, a 4-hour task is lost repeatedly.
- **Alignment task size is a parameter.** Default two 32-thread tasks per 72-thread node;
  `--whole_node_align` gives one full-node task (72 threads; same throughput, half the wall-clock
  per task) for clusters where the work at risk per outage matters more than fine packing.
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
- **MarkDuplicates memory** on 60×+ samples: raise the `markdup` label if retries occur. It is
  the longest remaining unscattered task (2–3 h on a 30× sample, single-threaded); switching
  to `MarkDuplicatesSpark` or `samtools markdup` would bring it under an hour if outages demand it.
- **Sex inference thresholds** (Y/autosome > 0.15 → M, < 0.05 → F) are conservative;
  somalier's call is reported alongside and disagreements are flagged.
- **CollectWgsMetrics** runs with `--USE_FAST_ALGORITHM`; compare its mean coverage with
  mosdepth once and keep only one if they agree.
- **Future modules**: SV/CNV (Manta/smoove, CNVkit in `wgs` mode), STR (ExpansionHunter),
  DeepVariant + GLnexus as an alternative caller, GIAB benchmarking when a control is
  sequenced.

## Advisory steps never abort a run (October 2026 pilot)

The 3-sample pilot exposed two ways a run could die for reasons unrelated to the data:
plink2 refuses LD pruning below 50 samples, and bcftools 1.21 segfaults when asked to
annotate an empty chromosome slice (chrM had no passing variants) with a contig-less
ClinVar file. PCA now adapts to the cohort size and degrades to empty outputs on failure,
empty slices bypass annotation, and ClinVar is reheadered with the reference contigs.
The head job also stops resubmitting when the same process fails in two consecutive
attempts: that is a deterministic error to fix, not an outage to wait out.

## Stage gathers emit per group, not per stage (October 2026 pilot)

Every per-sample or per-interval gather (`MARKDUP`, `GATHER_BQSR`, `GATHER_BAM`, `GATHER_GVCF`,
`GENOMICSDB_IMPORT`) used a plain `groupTuple()`. Without a known group size Nextflow holds all
groups until the upstream channel closes, which turned the workflow into a sequence of cohort-wide
barriers: no sample could start BQSR until the last sample's MarkDuplicates finished, no sample could
reach ApplyBQSR until every BQSR shard of every sample was done, and so on. In the pilot that left
three reserved nodes idle for over an hour while one MarkDuplicates task finished; for 67 samples the
idle tails would have added many hours.

Each gather now wraps its key with `groupKey(key, size)`, where the size is known up front: lanes per
sample from the sample sheet, BQSR groups and HaplotypeCaller intervals from counting the
`PREP_INTERVALS` outputs, and the number of samples for the per-interval joint-genotyping groups.
A group is emitted the moment its last member arrives, so fast samples flow through the whole
per-sample chain while slow ones are still aligning. Verified locally: the change does not alter
task hashes (a `-resume` after the rewiring reported every task cached), and the arrival order of
files within a group replays identically on resume, so no `sort:` option is needed.
