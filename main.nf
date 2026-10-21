// GenomeForge — whole-genome germline joint calling (FASTQ -> annotated, QC-ed cohort callset).
// Runs on the CallForge container images and reference bundle. See README.md.

include { CHECK_RESOURCES; PREP_INTERVALS; PREP_CLINVAR }                 from './modules/prep'
include { ALIGN_LANE; MARKDUP; BASERECAL; APPLYBQSR }                     from './modules/align'
include { SAMTOOLS_STATS; MOSDEPTH; WGS_METRICS; VERIFYBAMID; SOMALIER_EXTRACT } from './modules/qc'
include { HAPLOTYPECALLER; GATHER_GVCF }                                  from './modules/calling'
include { GENOMICSDB_IMPORT; GENOTYPE_GVCFS; CONCAT_VCFS; VQSR_RECAL; APPLY_VQSR; HARD_FILTER; NORMALIZE } from './modules/joint'
include { ANNOTATE_CHROM; VEP_CHROM; CONCAT_ANNOTATED; FLATTEN }          from './modules/annotate'
include { SOMALIER_RELATE; BCFTOOLS_STATS; PLINK_PCA; SUMMARIZE_QC; PLOT_QC; MULTIQC; VERSIONS } from './modules/cohort'

def vcfOf(p) { return file(p, checkIfExists: true) }
def tbiOf(p) { return file("${p}.tbi", checkIfExists: true) }

def helpMessage() {
    return """
    GenomeForge ${workflow.manifest.version} — whole-genome germline joint calling (CallForge images + refs)

    Usage:
      nextflow run main.nf -profile noguchi --input samplesheet.csv --outdir results
      nextflow run main.nf -profile slurm -params-file site.params.yaml --input samplesheet.csv

    Required:
      --input          sample sheet CSV (one row per lane-pair): sample_id,lane_id,library,flowcell,lane,fastq_1,fastq_2[,sex,phenotype]
                       build it with: bin/make_samplesheet.py --roots <delivery dirs> --out samplesheet.csv
    Key options (defaults in nextflow.config; site paths from a profile or params.example.yaml):
      --fasta          chr-style no-alt GRCh38 with bwa-mem2 index, .fai, .dict
      --calling_regions  GATK wgs_calling_regions.hg38.interval_list (bin/fetch_resources.sh) or BED
      --scatter_count  48      --filtering vqsr|hard      --trim true|false
      --run_vep / --run_pca / --run_contamination  true|false      --archive_cram false
    """.stripIndent()
}

workflow {

    if (params.help) {
        log.info(helpMessage())
        System.exit(0)
    }
    if (!params.input) {
        error("Missing --input sample sheet. Build one with bin/make_samplesheet.py (see README.md).")
    }
    def required = ['fasta', 'dbsnp', 'mills', 'known_indels', 'clinvar', 'gnomad_dir', 'somalier_sites', 'callforge_sif', 'cohortqc_sif']
    def missing  = required.findAll { !params[it] }
    if (missing) {
        error("Missing site parameters: ${missing.join(', ')}. Use -profile noguchi, or -profile slurm with a -params-file (see params.example.yaml).")
    }
    if (params.run_vep && !params.vep_cache) {
        error("--run_vep needs --vep_cache (or set run_vep false)")
    }
    if (params.run_pca && !params.burden_sif) {
        error("--run_pca needs --burden_sif (plink2) or set run_pca false")
    }

    // ---- reference bundle -------------------------------------------------------
    def fasta     = file(params.fasta, checkIfExists: true)
    def fai       = file("${params.fasta}.fai", checkIfExists: true)
    def dict      = file(params.fasta.replaceFirst(/\.(fa|fasta|fna)$/, '.dict'), checkIfExists: true)
    def ref       = channel.value([fasta, fai, dict])
    def bwa_index = channel.fromPath("${params.fasta}.{0123,amb,ann,bwt.2bit.64,pac}", checkIfExists: true).collect()

    def dbsnp     = channel.value([vcfOf(params.dbsnp), tbiOf(params.dbsnp)])
    def known_vcf = channel.value([params.dbsnp, params.mills, params.known_indels].collect { vcfOf(it) })
    def known_tbi = channel.value([params.dbsnp, params.mills, params.known_indels].collect { tbiOf(it) })
    def somalier_sites = channel.value(vcfOf(params.somalier_sites))      // somalier reads it sequentially; no index needed
    def chroms    = params.chromosomes.tokenize(',')

    def vqsr_paths = [params.hapmap, params.omni, params.g1000, params.mills, params.dbsnp]
    def use_vqsr   = params.filtering == 'vqsr' && vqsr_paths.every { it && file(it).exists() && file("${it}.tbi").exists() }
    if (params.filtering == 'vqsr' && !use_vqsr) {
        log.warn("VQSR training sets (hapmap/omni/1000G) not found: falling back to GATK hard filters. Run bin/fetch_resources.sh to stage them.")
    }
    def regions = params.calling_regions ? file(params.calling_regions, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    if (!params.calling_regions) {
        log.warn("No --calling_regions given: calling on whole primary chromosomes (centromeres/gaps included). bin/fetch_resources.sh stages the GATK wgs_calling_regions list.")
    }

    def gnomad_chr22 = "${params.gnomad_dir}/${params.gnomad_pattern.replace('{chrom}', 'chr22')}"
    def check_list   = ([params.dbsnp, params.mills, params.known_indels, params.somalier_sites, gnomad_chr22] +
                        (use_vqsr ? [params.hapmap, params.omni, params.g1000] : [])).findAll { file(it).exists() }
    def check_vcfs   = channel.value(check_list.collect { file(it) })
    def check_tbis   = channel.value(check_list.findAll { file("${it}.tbi").exists() }.collect { file("${it}.tbi") })

    // ---- sample sheet: one row per lane-pair ------------------------------------
    def samplesheet = file(params.input, checkIfExists: true)
    def lanes = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample_id, lane: row.lane_id, library: row.library ?: 'NA',
                        flowcell: row.flowcell ?: 'NA', lane_no: row.lane ?: 'NA']
            tuple(meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }

    // ---- stage 0 ------------------------------------------------------------------
    CHECK_RESOURCES(ref, bwa_index, check_vcfs, check_tbis)
    PREP_INTERVALS(ref, regions)
    def intervals = PREP_INTERVALS.out.intervals.flatten().map { f -> tuple(f.name.replace('-scattered.interval_list', ''), f) }
    PREP_CLINVAR(file(params.clinvar, checkIfExists: true), ref)

    // ---- alignment -> analysis-ready BAM ---------------------------------------
    ALIGN_LANE(lanes, ref, bwa_index)
    def per_sample = ALIGN_LANE.out.bam.map { meta, bam, bai -> tuple(meta.id, bam, bai) }.groupTuple()
    MARKDUP(per_sample)
    BASERECAL(MARKDUP.out.bam, ref, known_vcf, known_tbi)
    APPLYBQSR(MARKDUP.out.bam.join(BASERECAL.out), ref)
    def bams = APPLYBQSR.out.bam

    // ---- per-sample QC -----------------------------------------------------------
    SAMTOOLS_STATS(bams)
    MOSDEPTH(bams)
    WGS_METRICS(bams, ref)
    SOMALIER_EXTRACT(bams, ref, somalier_sites)
    def selfsm = channel.empty()
    if (params.run_contamination) {
        VERIFYBAMID(bams, ref)
        selfsm = VERIFYBAMID.out.selfsm
    }

    // ---- variant calling ---------------------------------------------------------
    HAPLOTYPECALLER(bams.combine(intervals), ref)
    GATHER_GVCF(HAPLOTYPECALLER.out.map { s, iv, g, t -> tuple(s, g, t) }.groupTuple())
    def by_interval = HAPLOTYPECALLER.out.map { s, iv, g, t -> tuple(iv, s, g, t) }.groupTuple().join(intervals)
    GENOMICSDB_IMPORT(by_interval)
    GENOTYPE_GVCFS(GENOMICSDB_IMPORT.out, ref, dbsnp)
    CONCAT_VCFS(GENOTYPE_GVCFS.out.map { iv, v, t -> v }.collect(), GENOTYPE_GVCFS.out.map { iv, v, t -> t }.collect())

    def filtered = channel.empty()
    if (use_vqsr) {
        VQSR_RECAL(CONCAT_VCFS.out, ['SNP', 'INDEL'], ref,
                   vcfOf(params.hapmap), vcfOf(params.omni), vcfOf(params.g1000), vcfOf(params.mills), vcfOf(params.dbsnp),
                   channel.value(vqsr_paths.collect { tbiOf(it) }))
        APPLY_VQSR(CONCAT_VCFS.out, VQSR_RECAL.out.map { m, r, i, t -> [r, i, t] }.flatten().collect(), ref)
        filtered = APPLY_VQSR.out
    } else {
        HARD_FILTER(CONCAT_VCFS.out, ref)
        filtered = HARD_FILTER.out
    }
    NORMALIZE(filtered, ref)

    // ---- annotation (scattered by chromosome) ------------------------------------
    ANNOTATE_CHROM(NORMALIZE.out.pass.combine(channel.fromList(chroms)), dbsnp, PREP_CLINVAR.out, file(params.gnomad_dir, checkIfExists: true))
    def annotated = ANNOTATE_CHROM.out
    if (params.run_vep) {
        VEP_CHROM(ANNOTATE_CHROM.out, ref, file(params.vep_cache, checkIfExists: true))
        annotated = VEP_CHROM.out
    }
    CONCAT_ANNOTATED(annotated.map { c, v, t -> v }.collect(), annotated.map { c, v, t -> t }.collect(), chroms)
    FLATTEN(CONCAT_ANNOTATED.out)

    // ---- cohort QC ---------------------------------------------------------------
    SOMALIER_RELATE(SOMALIER_EXTRACT.out.collect())
    BCFTOOLS_STATS(NORMALIZE.out.pass, ref)
    def pca_files = channel.value([])
    if (params.run_pca) {
        PLINK_PCA(NORMALIZE.out.pass)
        pca_files = PLINK_PCA.out.eigenvec.mix(PLINK_PCA.out.eigenval).collect()
    }
    def tranche_files = use_vqsr ? VQSR_RECAL.out.map { m, r, i, t -> t }.collect() : channel.value([])
    SUMMARIZE_QC(SAMTOOLS_STATS.out.flagstat.collect(), MARKDUP.out.metrics.collect(), MOSDEPTH.out.summary.collect(),
                 WGS_METRICS.out.collect(), selfsm.collect().ifEmpty([]),
                 SOMALIER_RELATE.out.samples, SOMALIER_RELATE.out.pairs, BCFTOOLS_STATS.out, samplesheet)
    PLOT_QC(SUMMARIZE_QC.out.tsv, samplesheet,
            MOSDEPTH.out.summary.mix(MOSDEPTH.out.dist).collect(), SAMTOOLS_STATS.out.stats.collect(),
            ALIGN_LANE.out.fastp.collect().ifEmpty([]), SOMALIER_RELATE.out.pairs, pca_files, BCFTOOLS_STATS.out,
            tranche_files, FLATTEN.out)
    def mqc = ALIGN_LANE.out.fastp.mix(SAMTOOLS_STATS.out.flagstat, SAMTOOLS_STATS.out.stats, SAMTOOLS_STATS.out.idxstats,
                                       MARKDUP.out.metrics, MOSDEPTH.out.summary, MOSDEPTH.out.dist, WGS_METRICS.out,
                                       selfsm, SOMALIER_RELATE.out.samples, SOMALIER_RELATE.out.pairs, BCFTOOLS_STATS.out)
    MULTIQC(mqc.collect())
    VERSIONS(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params)))
}
