// Stage 0: resource validation, calling-interval scatter, ClinVar contig rename.

process CHECK_RESOURCES {
    label 'small'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    input:
    tuple path(fasta), path(fai), path(dict)
    path(bwa_index)
    path(vcfs)
    path(tbis)

    output:
    path("resource_check.txt")

    script:
    """
    {
      echo "reference: ${fasta}"
      nctg=\$(wc -l < ${fai})
      nalt=\$(cut -f1 ${fai} | grep -cE '_alt|_decoy|hs38d1|EBV|HLA-' || true)
      style=\$(cut -f1 ${fai} | head -1 | grep -q '^chr' && echo chr || echo nochr)
      echo "contigs=\$nctg alt_or_decoy=\$nalt contig_style=\$style"
      [ "\$nalt" -eq 0 ] || echo "WARNING: reference carries alt/decoy contigs"
      for ext in 0123 amb ann bwt.2bit.64 pac; do
        [ -s "${fasta}.\$ext" ] && echo "bwa-mem2 index ok: .\$ext" || echo "MISSING bwa-mem2 index part: ${fasta}.\$ext"
      done
      for v in ${vcfs}; do
        c=\$( (tabix -l "\$v" 2>/dev/null || bcftools view -H "\$v" 2>/dev/null | cut -f1) | head -1 )
        s=\$(echo "\$c" | grep -q '^chr' && echo chr || echo nochr)
        echo "vcf \$v first_contig=\$c style=\$s"
        [ "\$s" = "\$style" ] || echo "MISMATCH: \$v contig style (\$s) differs from the reference (\$style)"
      done
    } > resource_check.txt
    cat resource_check.txt
    if grep -qE '^(MISSING|MISMATCH)' resource_check.txt; then echo "resource check FAILED" >&2; exit 1; fi
    """

    stub:
    """
    echo stub > resource_check.txt
    """
}

process PREP_INTERVALS {
    label 'medium'
    publishDir "${params.outdir}/pipeline_info/intervals", mode: params.publish_mode

    input:
    tuple path(fasta), path(fai), path(dict)
    path(regions)

    output:
    path("calling_regions.bed"),          emit: bed
    path("scattered/*.interval_list"),    emit: intervals

    script:
    def have_regions = regions.name != 'NO_FILE'
    def mode = params.scatter_mode ?: (have_regions ? 'BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW' : 'INTERVAL_SUBDIVISION')
    def mem  = Math.max(2, (task.memory.toGiga() * 0.8) as int)
    """
    if [ "${have_regions}" = "false" ]; then
        awk 'BEGIN{OFS="\\t"} \$1 ~ /^chr([0-9]+|X|Y|M)\$/ {print \$1,0,\$2}' ${fai} > calling_regions.bed
    elif [[ "${regions}" == *.interval_list ]]; then
        grep -v '^@' ${regions} | awk 'BEGIN{OFS="\\t"} \$1 ~ /^chr([0-9]+|X|Y|M)\$/ {print \$1,\$2-1,\$3}' > calling_regions.bed
    else
        awk 'BEGIN{OFS="\\t"} \$1 ~ /^chr([0-9]+|X|Y|M)\$/ {print \$1,\$2,\$3}' ${regions} > calling_regions.bed
    fi
    [ -s calling_regions.bed ] || { echo "no calling regions derived" >&2; exit 1; }
    mkdir -p scattered
    gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" SplitIntervals \\
        -R ${fasta} -L calling_regions.bed --scatter-count ${params.scatter_count} \\
        --subdivision-mode ${mode} -O scattered
    """

    stub:
    """
    touch calling_regions.bed
    mkdir -p scattered
    for i in \$(seq 0 \$((${params.scatter_count} - 1))); do touch scattered/\$(printf '%04d' \$i)-scattered.interval_list; done
    """
}

process PREP_CLINVAR {
    label 'medium'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode, pattern: "clinvar.chr.vcf.gz*"

    input:
    path(clinvar)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple path("clinvar.chr.vcf.gz"), path("clinvar.chr.vcf.gz.tbi")

    script:
    """
    { for i in \$(seq 1 22) X Y; do echo "\$i chr\$i"; done; echo "MT chrM"; } > rename.txt
    bcftools annotate --rename-chrs rename.txt ${clinvar} -Ou \\
      | bcftools sort -T ./tmp -m ${Math.max(1, (task.memory.toGiga() * 0.5) as int)}G -Oz -o clinvar.chr.vcf.gz
    tabix -p vcf clinvar.chr.vcf.gz
    """

    stub:
    """
    touch clinvar.chr.vcf.gz clinvar.chr.vcf.gz.tbi
    """
}
