// Stage 10: cohort-level QC (relatedness/sex, variant stats, PCA), QC summary, MultiQC, provenance.

process SOMALIER_RELATE {
    label 'small'
    publishDir "${params.outdir}/07_cohort_qc/somalier", mode: params.publish_mode

    input:
    path(somaliers)

    output:
    path("cohort.samples.tsv"), emit: samples
    path("cohort.pairs.tsv"),   emit: pairs
    path("cohort.html"),        emit: html

    script:
    """
    somalier relate --infer -o cohort ${somaliers}
    """

    stub:
    """
    touch cohort.samples.tsv cohort.pairs.tsv cohort.html
    """
}

process BCFTOOLS_STATS {
    label 'medium'
    publishDir "${params.outdir}/07_cohort_qc", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    path("cohort.bcftools_stats.txt")

    script:
    """
    bcftools stats --threads ${task.cpus} -s - --fasta-ref ${fasta} ${vcf} > cohort.bcftools_stats.txt
    """

    stub:
    """
    touch cohort.bcftools_stats.txt
    """
}

process PLINK_PCA {
    label 'medium'
    publishDir "${params.outdir}/07_cohort_qc/pca", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)

    output:
    path("cohort_pca.eigenvec"), emit: eigenvec
    path("cohort_pca.eigenval"), emit: eigenval
    path("cohort_pca.log"),      emit: log

    script:
    def mem = Math.max(2000, ((task.memory.toMega() * 0.8) as int))
    """
    # PCA is advisory: it must never abort a run. Small cohorts (plink2 refuses LD pruning
    # below 50 samples and cannot return more PCs than samples-1) are handled, and any
    # plink2 failure leaves empty outputs plus an explanation in cohort_pca.log.
    n=\$(gzip -dc ${vcf} | grep -m1 '^#CHROM' | awk '{print NF - 9}')
    pcs=10; [ "\$n" -le 10 ] && pcs=\$((n - 1))
    badld=""; [ "\$n" -lt 50 ] && badld="--bad-ld"
    : > cohort_pca.eigenvec; : > cohort_pca.eigenval
    if [ "\$n" -lt 3 ]; then
        echo "PCA skipped: \$n sample(s) in the cohort (need >= 3)" | tee cohort_pca.log >&2
        exit 0
    fi
    set +e
    plink2 --vcf ${vcf} --vcf-half-call m --snps-only just-acgt --max-alleles 2 --autosome \\
        --set-missing-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 50 missing \\
        --maf 0.05 --geno 0.05 --hwe 1e-6 --indep-pairwise 200 50 0.2 \$badld \\
        --threads ${task.cpus} --memory ${mem} --out prune \\
    && plink2 --vcf ${vcf} --vcf-half-call m --snps-only just-acgt --max-alleles 2 --autosome \\
        --set-missing-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 50 missing \\
        --extract prune.prune.in --pca \$pcs --threads ${task.cpus} --memory ${mem} --out cohort_pca
    rc=\$?
    set -e
    if [ \$rc -ne 0 ] || [ ! -s cohort_pca.eigenvec ]; then
        echo "plink2 PCA failed (exit \$rc) for n=\$n samples; PCA figures will be skipped" | tee -a cohort_pca.log >&2
        : > cohort_pca.eigenvec; : > cohort_pca.eigenval
    fi
    touch cohort_pca.log
    """

    stub:
    """
    touch cohort_pca.eigenvec cohort_pca.eigenval cohort_pca.log
    """
}

process SUMMARIZE_QC {
    label 'small'
    publishDir "${params.outdir}/08_reports", mode: params.publish_mode

    input:
    path(flagstats,  stageAs: 'flagstat/*')
    path(markdups,   stageAs: 'markdup/*')
    path(mosdepths,  stageAs: 'mosdepth/*')
    path(wgsmetrics, stageAs: 'wgs_metrics/*')
    path(selfsms,    stageAs: 'verifybamid/*')
    path(somalier_samples)
    path(somalier_pairs)
    path(bcfstats)
    path(samplesheet)

    output:
    path("cohort_qc_summary.tsv"), emit: tsv
    path("cohort_qc_summary.md"),  emit: md

    script:
    """
    summarize_qc.py --flagstat flagstat --markdup markdup --mosdepth mosdepth --wgs-metrics wgs_metrics \\
        --verifybamid verifybamid --somalier-samples ${somalier_samples} --somalier-pairs ${somalier_pairs} \\
        --bcftools-stats ${bcfstats} --samplesheet ${samplesheet} \\
        --min-depth ${params.qc_min_depth} --max-dup ${params.qc_max_dup} --max-contam ${params.qc_max_contam} \\
        --min-mapped ${params.qc_min_mapped} --out-tsv cohort_qc_summary.tsv --out-md cohort_qc_summary.md
    """

    stub:
    """
    touch cohort_qc_summary.tsv cohort_qc_summary.md
    """
}

process MULTIQC {
    label 'medium'
    publishDir "${params.outdir}/08_reports/multiqc", mode: params.publish_mode

    input:
    path(files, stageAs: 'inputs/*')

    output:
    path("cohort_multiqc.html"),    emit: html
    path("cohort_multiqc_data"),    emit: data

    script:
    """
    multiqc -f -n cohort_multiqc.html inputs
    """

    stub:
    """
    touch cohort_multiqc.html; mkdir -p cohort_multiqc_data
    """
}

process PLOT_QC {
    label 'medium'
    publishDir "${params.outdir}/08_reports", mode: params.publish_mode

    input:
    path(qc_summary)
    path(samplesheet)
    path(mosdepth,  stageAs: 'mosdepth/*')
    path(samtools,  stageAs: 'samtools/*')
    path(fastp,     stageAs: 'fastp/*')
    path(somalier_pairs)
    path(pca,       stageAs: 'pca/*')
    path(bcfstats)
    path(tranches,  stageAs: 'vqsr/*')
    path(flat)

    output:
    path("figures"),           emit: figures
    path("qc_dashboard.html"), emit: dashboard
    path("qc_report.pdf"),     emit: pdf

    script:
    """
    mkdir -p mosdepth samtools fastp pca vqsr
    plot_qc.py --qc-summary ${qc_summary} --samplesheet ${samplesheet} \
        --mosdepth mosdepth --samtools samtools --fastp fastp --somalier-pairs ${somalier_pairs} \
        --pca-eigenvec pca/cohort_pca.eigenvec --pca-eigenval pca/cohort_pca.eigenval \
        --bcftools-stats ${bcfstats} --tranches vqsr --flat ${flat} \
        --min-depth ${params.qc_min_depth} --max-dup ${params.qc_max_dup} --max-contam ${params.qc_max_contam} --min-mapped ${params.qc_min_mapped} \
        --out-dir figures --dashboard qc_dashboard.html --pdf qc_report.pdf --label "run ${workflow.runName} · ${workflow.start}"
    """

    stub:
    """
    mkdir -p figures; touch figures/coverage_per_sample.png qc_dashboard.html qc_report.pdf
    """
}

process VERSIONS {
    label 'small'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    input:
    val(params_json)

    output:
    path("versions.txt")
    path("params.json")

    script:
    """
    cat > params.json <<'EOF'
${params_json}
EOF
    {
      echo "GenomeForge ${workflow.manifest.version}  nextflow ${workflow.nextflow.version}  run ${workflow.runName}  started ${workflow.start}"
      echo "container: ${params.callforge_sif}"
      echo "bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -1)"
      echo "samtools: \$(samtools --version | head -1)"
      echo "bcftools: \$(bcftools --version | head -1)"
      echo "gatk: \$(gatk --version 2>/dev/null | grep -i 'GATK v' | head -1 || gatk --version 2>&1 | tail -1)"
      echo "fastp: \$(fastp --version 2>&1)"
      echo "mosdepth: \$(mosdepth --version)"
      echo "multiqc: \$(multiqc --version 2>&1 | head -1)"
    } > versions.txt
    """

    stub:
    """
    echo stub > versions.txt; echo '{}' > params.json
    """
}
