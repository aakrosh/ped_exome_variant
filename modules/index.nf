process BGZIP_INDEX_VARIANTS {
    tag "bgzip tabix"
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'
    publishDir 'results', mode: 'copy'

    input:
    path vcf_file

    output:
    tuple path("${vcf_file}.gz"), path("${vcf_file}.gz.tbi")

    script:
    """
    bgzip -@ ${task.cpus} ${vcf_file}
    tabix ${vcf_file}.gz
    """
}