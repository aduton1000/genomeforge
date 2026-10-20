// Stage 5: scattered HaplotypeCaller (gVCF) and per-sample gVCF gather.

process HAPLOTYPECALLER {
    tag "${sample}:${iv}"
    label 'hc'

    input:
    tuple val(sample), path(bam), path(bai), val(iv), path(interval)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(sample), val(iv), path("${sample}.${iv}.g.vcf.gz"), path("${sample}.${iv}.g.vcf.gz.tbi")

    script:
    def mem = Math.max(3, (task.memory.toGiga() * 0.75) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" HaplotypeCaller \\
        -R ${fasta} -I ${bam} -L ${interval} -O ${sample}.${iv}.g.vcf.gz \\
        -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation \\
        --native-pair-hmm-threads ${task.cpus}
    """

    stub:
    """
    touch ${sample}.${iv}.g.vcf.gz ${sample}.${iv}.g.vcf.gz.tbi
    """
}

process GATHER_GVCF {
    tag "${sample}"
    label 'small'
    publishDir "${params.outdir}/04_gvcf", mode: params.publish_mode

    input:
    tuple val(sample), path(gvcfs), path(tbis)

    output:
    tuple val(sample), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi")

    script:
    def list   = (gvcfs instanceof List ? gvcfs : [gvcfs]).sort(false) { it.name }
    def inputs = list.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx3g -Djava.io.tmpdir=./tmp" GatherVcfs ${inputs} -O ${sample}.g.vcf.gz
    tabix -p vcf ${sample}.g.vcf.gz
    """

    stub:
    """
    touch ${sample}.g.vcf.gz ${sample}.g.vcf.gz.tbi
    """
}
