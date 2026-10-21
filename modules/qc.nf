// Stage 4: per-sample alignment QC.

process SAMTOOLS_STATS {
    tag "${sample}"
    label 'qc'
    publishDir "${params.outdir}/03_sample_qc/samtools", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path("${sample}.flagstat"), emit: flagstat
    path("${sample}.stats"),    emit: stats
    path("${sample}.idxstats"), emit: idxstats

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${sample}.flagstat
    samtools stats    -@ ${task.cpus} ${bam} > ${sample}.stats
    samtools idxstats ${bam} > ${sample}.idxstats
    """

    stub:
    """
    touch ${sample}.flagstat ${sample}.stats ${sample}.idxstats
    """
}

process MOSDEPTH {
    tag "${sample}"
    label 'qc'
    publishDir "${params.outdir}/03_sample_qc/mosdepth", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path("${sample}.mosdepth.summary.txt"),     emit: summary
    path("${sample}.mosdepth.global.dist.txt"), emit: dist
    path("${sample}.mosdepth.region.dist.txt"), emit: region_dist
    path("${sample}.regions.bed.gz*"),          emit: regions

    script:
    """
    mosdepth -t ${task.cpus} -n --fast-mode --by 1000 --mapq 20 ${sample} ${bam}
    """

    stub:
    """
    touch ${sample}.mosdepth.summary.txt ${sample}.mosdepth.global.dist.txt ${sample}.mosdepth.region.dist.txt ${sample}.regions.bed.gz ${sample}.regions.bed.gz.csi
    """
}

process WGS_METRICS {
    tag "${sample}"
    label 'wgsmetrics'
    publishDir "${params.outdir}/03_sample_qc/wgs_metrics", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    path("${sample}.wgs_metrics.txt")

    script:
    def mem = Math.max(4, (task.memory.toGiga() * 0.8) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" CollectWgsMetrics \\
        -I ${bam} -O ${sample}.wgs_metrics.txt -R ${fasta} \\
        --USE_FAST_ALGORITHM true --INCLUDE_BQ_HISTOGRAM false --VALIDATION_STRINGENCY SILENT
    """

    stub:
    """
    touch ${sample}.wgs_metrics.txt
    """
}

process VERIFYBAMID {
    tag "${sample}"
    label 'qc'
    publishDir "${params.outdir}/03_sample_qc/verifybamid", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    path("${sample}.selfSM"),          emit: selfsm
    path("${sample}.verifybamid.log"), emit: log

    script:
    def svd = params.verifybamid_svd ?: ''
    """
    svd="${svd}"
    if [ -z "\$svd" ]; then
        f=\$(ls /opt/conda/share/verifybamid2*/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD 2>/dev/null | head -1 || true)
        [ -n "\$f" ] && svd="\${f%.UD}"
    fi
    if [ -n "\$svd" ] && [ -s "\$svd.UD" ]; then
        verifybamid2 --SVDPrefix "\$svd" --Reference ${fasta} --BamFile ${bam} --NumPC 4 --NumThread ${task.cpus} \\
            --Output ${sample} > ${sample}.verifybamid.log 2>&1
    else
        echo "verifybamid2 SVD resource not found (set --verifybamid_svd); contamination not estimated" > ${sample}.verifybamid.log
        printf '#SEQ_ID\\tRG\\tCHIP_ID\\t#SNPS\\t#READS\\tAVG_DP\\tFREEMIX\\n${sample}\\tALL\\tNA\\t0\\t0\\t0\\tNA\\n' > ${sample}.selfSM
    fi
    """

    stub:
    """
    printf '#SEQ_ID\\tRG\\tCHIP_ID\\t#SNPS\\t#READS\\tAVG_DP\\tFREEMIX\\n${sample}\\tALL\\tNA\\t0\\t0\\t0\\t0.001\\n' > ${sample}.selfSM
    touch ${sample}.verifybamid.log
    """
}

process SOMALIER_EXTRACT {
    tag "${sample}"
    label 'small'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path(sites)

    output:
    path("*.somalier")

    script:
    """
    somalier extract -d . --sites ${sites} -f ${fasta} ${bam}
    """

    stub:
    """
    touch ${sample}.somalier
    """
}
