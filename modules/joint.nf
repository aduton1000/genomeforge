// Stages 6-8: joint genotyping per shard, concatenation, VQSR (or hard filters), normalisation.

process GENOMICSDB_IMPORT {
    tag "${iv}"
    label 'gdb'

    input:
    tuple val(iv), val(samples), path(gvcfs), path(tbis), path(interval)

    output:
    tuple val(iv), path("gdb_${iv}"), path(interval)

    script:
    def slist = samples instanceof List ? samples : [samples]
    def glist = gvcfs   instanceof List ? gvcfs   : [gvcfs]
    def rows  = [slist, glist].transpose().collect { "${it[0]}\t${it[1].name}" }.join('\n')
    def mem   = Math.max(4, (task.memory.toGiga() * 0.7) as int)
    """
    cat > sample_map.txt <<'EOF'
${rows}
EOF
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" GenomicsDBImport \\
        --genomicsdb-workspace-path gdb_${iv} --sample-name-map sample_map.txt \\
        -L ${interval} --merge-input-intervals true --batch-size 16 \\
        --reader-threads ${task.cpus} --genomicsdb-shared-posixfs-optimizations true --tmp-dir ./tmp
    """

    stub:
    """
    mkdir -p gdb_${iv}; touch gdb_${iv}/callset.json
    """
}

process GENOTYPE_GVCFS {
    tag "${iv}"
    label 'genotype'

    input:
    tuple val(iv), path(gdb), path(interval)
    tuple path(fasta), path(fai), path(dict)
    tuple path(dbsnp), path(dbsnp_tbi)

    output:
    tuple val(iv), path("${iv}.vcf.gz"), path("${iv}.vcf.gz.tbi")

    script:
    def mem = Math.max(4, (task.memory.toGiga() * 0.75) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" GenotypeGVCFs \\
        -R ${fasta} -V gendb://${gdb} -L ${interval} -O ${iv}.vcf.gz -D ${dbsnp} \\
        -G StandardAnnotation --only-output-calls-starting-in-intervals \\
        --merge-input-intervals true --genomicsdb-shared-posixfs-optimizations true
    """

    stub:
    """
    touch ${iv}.vcf.gz ${iv}.vcf.gz.tbi
    """
}

process CONCAT_VCFS {
    label 'medium'
    publishDir "${params.outdir}/05_joint", mode: params.publish_mode

    input:
    path(vcfs)
    path(tbis)

    output:
    tuple path("cohort.raw.vcf.gz"), path("cohort.raw.vcf.gz.tbi")

    script:
    def list = (vcfs instanceof List ? vcfs : [vcfs]).sort(false) { it.name }
    """
    bcftools concat --threads ${task.cpus} -Oz -o cohort.raw.vcf.gz ${list.join(' ')}
    tabix -p vcf cohort.raw.vcf.gz
    """

    stub:
    """
    touch cohort.raw.vcf.gz cohort.raw.vcf.gz.tbi
    """
}

process VQSR_RECAL {
    tag "${mode}"
    label 'vqsr'
    publishDir "${params.outdir}/05_joint/vqsr", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)
    each mode
    tuple path(fasta), path(fai), path(dict)
    path(hapmap)
    path(omni)
    path(g1000)
    path(mills)
    path(dbsnp)
    path(tbis)

    output:
    tuple val(mode), path("${mode}.recal"), path("${mode}.recal.idx"), path("${mode}.tranches")

    script:
    def mem = Math.max(8, (task.memory.toGiga() * 0.8) as int)
    def res = mode == 'SNP' ?
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} " +
        "--resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} " +
        "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${g1000} " +
        "--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${dbsnp}" :
        "--resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} " +
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp}"
    def an   = mode == 'SNP' ? '-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP' :
                               '-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP'
    def maxg = mode == 'SNP' ? 6 : 4
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" VariantRecalibrator \\
        -R ${fasta} -V ${vcf} ${res} ${an} -mode ${mode} --max-gaussians ${maxg} --trust-all-polymorphic \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 \\
        -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \\
        -O ${mode}.recal --tranches-file ${mode}.tranches
    """

    stub:
    """
    touch ${mode}.recal ${mode}.recal.idx ${mode}.tranches
    """
}

process APPLY_VQSR {
    label 'vqsr'
    publishDir "${params.outdir}/05_joint", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)
    path(recal_files)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi")

    script:
    def mem = Math.max(4, (task.memory.toGiga() * 0.6) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" ApplyVQSR \\
        -R ${fasta} -V ${vcf} -O indel_recal.vcf.gz -mode INDEL \\
        --recal-file INDEL.recal --tranches-file INDEL.tranches \\
        --truth-sensitivity-filter-level ${params.indel_tranche} --create-output-variant-index true
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" ApplyVQSR \\
        -R ${fasta} -V indel_recal.vcf.gz -O cohort.filtered.vcf.gz -mode SNP \\
        --recal-file SNP.recal --tranches-file SNP.tranches \\
        --truth-sensitivity-filter-level ${params.snp_tranche} --create-output-variant-index true
    [ -s cohort.filtered.vcf.gz.tbi ] || tabix -p vcf cohort.filtered.vcf.gz
    """

    stub:
    """
    touch cohort.filtered.vcf.gz cohort.filtered.vcf.gz.tbi
    """
}

process HARD_FILTER {
    label 'vqsr'
    publishDir "${params.outdir}/05_joint", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi")

    script:
    def mem = Math.max(4, (task.memory.toGiga() * 0.6) as int)
    """
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" SelectVariants -R ${fasta} -V ${vcf} -select-type SNP -O snps.vcf.gz
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" SelectVariants -R ${fasta} -V ${vcf} -select-type INDEL -select-type MIXED -select-type MNP -O indels.vcf.gz
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" VariantFiltration -R ${fasta} -V snps.vcf.gz \\
        --filter-expression "${params.snp_filter_expr}" --filter-name HARD_SNP -O snps.filtered.vcf.gz
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" VariantFiltration -R ${fasta} -V indels.vcf.gz \\
        --filter-expression "${params.indel_filter_expr}" --filter-name HARD_INDEL -O indels.filtered.vcf.gz
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O cohort.filtered.vcf.gz
    [ -s cohort.filtered.vcf.gz.tbi ] || tabix -p vcf cohort.filtered.vcf.gz
    """

    stub:
    """
    touch cohort.filtered.vcf.gz cohort.filtered.vcf.gz.tbi
    """
}

process NORMALIZE {
    label 'medium'
    publishDir "${params.outdir}/05_joint", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple path("cohort.filtered.norm.vcf.gz"), path("cohort.filtered.norm.vcf.gz.tbi"), emit: all
    tuple path("cohort.pass.vcf.gz"),          path("cohort.pass.vcf.gz.tbi"),          emit: pass

    script:
    """
    bcftools norm --threads ${task.cpus} -m -any -f ${fasta} --check-ref w -Oz -o cohort.filtered.norm.vcf.gz ${vcf}
    tabix -p vcf cohort.filtered.norm.vcf.gz
    bcftools view --threads ${task.cpus} -f PASS -Oz -o cohort.pass.vcf.gz cohort.filtered.norm.vcf.gz
    tabix -p vcf cohort.pass.vcf.gz
    """

    stub:
    """
    touch cohort.filtered.norm.vcf.gz cohort.filtered.norm.vcf.gz.tbi cohort.pass.vcf.gz cohort.pass.vcf.gz.tbi
    """
}
