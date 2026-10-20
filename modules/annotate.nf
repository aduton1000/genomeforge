// Stage 9: per-chromosome annotation (dbSNP IDs, gnomAD exomes AF, ClinVar, VEP) and flattening.

process ANNOTATE_CHROM {
    tag "${chrom}"
    label 'medium'

    input:
    tuple path(vcf), path(tbi), val(chrom)
    tuple path(dbsnp), path(dbsnp_tbi)
    tuple path(clinvar), path(clinvar_tbi)
    path(gnomad_dir)

    output:
    tuple val(chrom), path("${chrom}.anno.vcf.gz"), path("${chrom}.anno.vcf.gz.tbi")

    script:
    def gfile  = params.gnomad_pattern.replace('{chrom}', chrom)
    def fields = params.gnomad_fields.tokenize(',').collect { "INFO/gnomAD_${it}:=${it}" }.join(',')
    """
    bcftools view --threads ${task.cpus} -r ${chrom} -Oz -o ${chrom}.in.vcf.gz ${vcf}
    tabix -p vcf ${chrom}.in.vcf.gz

    # dbSNP rsIDs
    bcftools annotate --threads ${task.cpus} -a ${dbsnp} -c ID -Oz -o ${chrom}.s1.vcf.gz ${chrom}.in.vcf.gz
    tabix -p vcf ${chrom}.s1.vcf.gz

    # gnomAD (per-chromosome sites file; absent for chrM)
    gn="${gnomad_dir}/${gfile}"
    if [ -s "\$gn" ]; then
        bcftools annotate --threads ${task.cpus} -a "\$gn" -c ${fields} -Oz -o ${chrom}.s2.vcf.gz ${chrom}.s1.vcf.gz
    else
        echo "no gnomAD file for ${chrom} (\$gn): gnomAD fields left empty" >&2
        cp ${chrom}.s1.vcf.gz ${chrom}.s2.vcf.gz
    fi
    tabix -p vcf ${chrom}.s2.vcf.gz

    # ClinVar (contigs renamed to chr-style in PREP_CLINVAR)
    bcftools annotate --threads ${task.cpus} -a ${clinvar} -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNSIGCONF,INFO/CLNVC \\
        -Oz -o ${chrom}.anno.vcf.gz ${chrom}.s2.vcf.gz
    tabix -p vcf ${chrom}.anno.vcf.gz
    rm -f ${chrom}.in.vcf.gz* ${chrom}.s1.vcf.gz* ${chrom}.s2.vcf.gz*
    """

    stub:
    """
    touch ${chrom}.anno.vcf.gz ${chrom}.anno.vcf.gz.tbi
    """
}

process VEP_CHROM {
    tag "${chrom}"
    label 'vep'

    input:
    tuple val(chrom), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    path(vep_cache)

    output:
    tuple val(chrom), path("${chrom}.vep.vcf.gz"), path("${chrom}.vep.vcf.gz.tbi")

    script:
    """
    vep --input_file ${vcf} --output_file ${chrom}.vep.vcf.gz --vcf --compress_output bgzip --format vcf \\
        --cache --offline --dir_cache ${vep_cache} --cache_version ${params.vep_release} \\
        --assembly GRCh38 --species homo_sapiens --fasta ${fasta} \\
        --fork ${task.cpus} --buffer_size 10000 --no_stats --force_overwrite \\
        --symbol --canonical --mane --biotype --hgvs --hgvsg --protein --numbers --domains --variant_class \\
        --sift b --polyphen b --af --max_af --af_gnomade --af_gnomadg --pubmed --flag_pick ${params.vep_extra}
    tabix -p vcf ${chrom}.vep.vcf.gz
    """

    stub:
    """
    touch ${chrom}.vep.vcf.gz ${chrom}.vep.vcf.gz.tbi
    """
}

process CONCAT_ANNOTATED {
    label 'medium'
    publishDir "${params.outdir}/06_annotated", mode: params.publish_mode

    input:
    path(vcfs)
    path(tbis)
    val(chroms)

    output:
    tuple path("cohort.annotated.vcf.gz"), path("cohort.annotated.vcf.gz.tbi")

    script:
    def list = (vcfs instanceof List ? vcfs : [vcfs]).sort(false) { chroms.indexOf(it.name.tokenize('.')[0]) }
    """
    bcftools concat --threads ${task.cpus} -Oz -o cohort.annotated.vcf.gz ${list.join(' ')}
    tabix -p vcf cohort.annotated.vcf.gz
    """

    stub:
    """
    touch cohort.annotated.vcf.gz cohort.annotated.vcf.gz.tbi
    """
}

process FLATTEN {
    label 'medium'
    publishDir "${params.outdir}/06_annotated", mode: params.publish_mode

    input:
    tuple path(vcf), path(tbi)

    output:
    path("variants.flat.tsv.gz")

    script:
    def gcols = params.gnomad_fields.tokenize(',').collect { "%gnomAD_${it}" }.join('\\t')
    def vepcols = params.run_vep ? '\\t%SYMBOL\\t%Gene\\t%Consequence\\t%IMPACT\\t%Feature\\t%HGVSc\\t%HGVSp\\t%MANE_SELECT\\t%CANONICAL\\t%SIFT\\t%PolyPhen\\t%MAX_AF' : ''
    def plugin  = params.run_vep ? '+split-vep -s worst' : 'query'
    """
    bcftools ${plugin} -H \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%AC\\t%AN\\t%AF${vepcols}\\t${gcols}\\t%CLNSIG\\t%CLNDN\\t%CLNREVSTAT[\\t%GT]\\n' \\
        ${vcf} | bgzip -@ ${task.cpus} > variants.flat.tsv.gz
    """

    stub:
    """
    touch variants.flat.tsv.gz
    """
}
