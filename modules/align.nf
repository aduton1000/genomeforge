// Stages 1-3: per-lane alignment, per-sample duplicate marking, BQSR scattered over contig groups.

process ALIGN_LANE {
    tag "${meta.id}:${meta.lane}"
    label 'align'
    publishDir "${params.outdir}/01_trim_reports", mode: params.publish_mode, pattern: "*.fastp.{json,html}"
    publishDir "${params.outdir}/01_trim_reports/logs", mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(meta), path(r1), path(r2)
    tuple path(fasta), path(fai), path(dict)
    path(bwa_index)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.lane.bam"), path("${meta.id}.${meta.lane}.lane.bam.bai"), emit: bam
    path("*.fastp.json"), emit: fastp, optional: true
    path("*.{fastp.html,log}"), emit: logs, optional: true

    script:
    def prefix      = "${meta.id}.${meta.lane}"
    def cn          = params.sequencing_center ? "\\tCN:${params.sequencing_center}" : ''
    def rg          = "@RG\\tID:${meta.lane}\\tSM:${meta.id}\\tLB:${meta.library}\\tPU:${meta.flowcell}.${meta.lane_no}\\tPL:ILLUMINA${cn}"
    def sort_cpus   = Math.max(4, Math.min(8, (task.cpus / 8) as int))      // 4 at 32 cpus, 8 at 64
    def trim        = params.trim.toString().toLowerCase() in ['true', '1', 'yes']
    def fastp_cpus  = trim ? sort_cpus : 0
    def bwa_cpus    = Math.max(4, task.cpus - sort_cpus - fastp_cpus)
    if (trim)
        """
        fastp --in1 ${r1} --in2 ${r2} --stdout --thread ${fastp_cpus} \\
              --qualified_quality_phred ${params.fastp_qual} --length_required ${params.fastp_min_len} \\
              --detect_adapter_for_pe --json ${prefix}.fastp.json --html ${prefix}.fastp.html 2> ${prefix}.fastp.log \\
          | bwa-mem2 mem -t ${bwa_cpus} -K 100000000 -Y -p -R '${rg}' ${fasta} - 2> ${prefix}.bwa.log \\
          | samtools sort -@ ${sort_cpus} -m 2G -T ./tmp/${prefix} -o ${prefix}.lane.bam -
        samtools index -@ ${sort_cpus} ${prefix}.lane.bam
        """
    else
        """
        bwa-mem2 mem -t ${bwa_cpus} -K 100000000 -Y -R '${rg}' ${fasta} ${r1} ${r2} 2> ${prefix}.bwa.log \\
          | samtools sort -@ ${sort_cpus} -m 2G -T ./tmp/${prefix} -o ${prefix}.lane.bam -
        samtools index -@ ${sort_cpus} ${prefix}.lane.bam
        """

    stub:
    """
    touch ${meta.id}.${meta.lane}.lane.bam ${meta.id}.${meta.lane}.lane.bam.bai ${meta.id}.${meta.lane}.fastp.json ${meta.id}.${meta.lane}.fastp.html ${meta.id}.${meta.lane}.bwa.log
    """
}

process MARKDUP {
    tag "${sample}"
    label 'markdup'
    publishDir "${params.outdir}/03_sample_qc/markdup", mode: params.publish_mode, pattern: "*.metrics.txt"

    input:
    tuple val(sample), path(bams), path(bais)

    output:
    tuple val(sample), path("${sample}.md.bam"), path("${sample}.md.bam.bai"), emit: bam
    path("${sample}.markdup.metrics.txt"), emit: metrics

    script:
    def list   = bams instanceof List ? bams : [bams]
    def inputs = list.collect { "-I ${it}" }.join(' ')
    def mem    = Math.max(4, (task.memory.toGiga() * 0.75) as int)
    """
    gatk --java-options "-Xmx${mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=./tmp" MarkDuplicates \\
        ${inputs} -O ${sample}.md.bam -M ${sample}.markdup.metrics.txt \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX false --VALIDATION_STRINGENCY SILENT --MAX_RECORDS_IN_RAM 4000000 --TMP_DIR ./tmp
    samtools index -@ ${task.cpus} ${sample}.md.bam
    """

    stub:
    """
    touch ${sample}.md.bam ${sample}.md.bam.bai ${sample}.markdup.metrics.txt
    """
}

process BASERECAL_SHARD {
    tag "${sample}:${grp}"
    label 'bqsr_shard'

    input:
    tuple val(sample), path(bam), path(bai), val(grp), path(intervals)
    tuple path(fasta), path(fai), path(dict)
    path(known)
    path(known_tbi)

    output:
    tuple val(sample), val(grp), path("${sample}.${grp}.recal.table")

    script:
    def list = known instanceof List ? known : [known]
    def ks   = list.collect { "--known-sites ${it}" }.join(' ')
    def mem  = Math.max(3, (task.memory.toGiga() * 0.75) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" BaseRecalibrator \\
        -R ${fasta} -I ${bam} -L ${intervals} ${ks} -O ${sample}.${grp}.recal.table
    """

    stub:
    """
    touch ${sample}.${grp}.recal.table
    """
}

process GATHER_BQSR {
    tag "${sample}"
    label 'small'
    publishDir "${params.outdir}/03_sample_qc/bqsr", mode: params.publish_mode

    input:
    tuple val(sample), path(tables)

    output:
    tuple val(sample), path("${sample}.recal.table")

    script:
    def list   = (tables instanceof List ? tables : [tables]).sort(false) { it.name }
    def inputs = list.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx3g -Djava.io.tmpdir=./tmp" GatherBQSRReports ${inputs} -O ${sample}.recal.table
    """

    stub:
    """
    touch ${sample}.recal.table
    """
}

process APPLYBQSR_SHARD {
    tag "${sample}:${grp}"
    label 'bqsr_shard'

    input:
    tuple val(sample), path(bam), path(bai), path(table), val(grp), path(intervals)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(sample), val(grp), path("${sample}.${grp}.bam")

    script:
    def mem      = Math.max(3, (task.memory.toGiga() * 0.75) as int)
    def unmapped = intervals.name.endsWith('.last.intervals') ? '-L unmapped' : ''
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" ApplyBQSR \\
        -R ${fasta} -I ${bam} --bqsr-recal-file ${table} -L ${intervals} ${unmapped} \\
        -O ${sample}.${grp}.bam --create-output-bam-index false
    """

    stub:
    """
    touch ${sample}.${grp}.bam
    """
}

process GATHER_BAM {
    tag "${sample}"
    label 'gather'
    publishDir "${params.outdir}/02_alignments", mode: params.publish_mode, pattern: "*.{bam,bai,cram,crai}"

    input:
    tuple val(sample), path(bams)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    path("${sample}.cram*"), emit: cram, optional: true

    script:
    def list   = (bams instanceof List ? bams : [bams]).sort(false) { it.name }   // 0000 … NNNN.last = reference order, unmapped last
    def inputs = list.collect { "-I ${it}" }.join(' ')
    def mem    = Math.max(4, (task.memory.toGiga() * 0.75) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" GatherBamFiles ${inputs} -O ${sample}.bam --CREATE_INDEX false
    samtools index -@ ${task.cpus} ${sample}.bam
    if [ "${params.archive_cram}" = "true" ]; then
        samtools view -@ ${task.cpus} -C -T ${fasta} -o ${sample}.cram ${sample}.bam
        samtools index ${sample}.cram
    fi
    """

    stub:
    """
    touch ${sample}.bam ${sample}.bam.bai
    """
}
