// Description: A Nextflow script to run a workflow that aligns a CRAM file using BWA-MEM.
nextflow.enable.dsl = 2

process EXTRACT_FASTQ {
    tag "$sample"
    container 'aakrosh/bwa-suite:alpine'

    input:
    tuple val(sample), path(cram)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz") 
    
    script:
    """
    samtools sort -n -@ ${task.cpus} ${cram} \
    | samtools fastq -@ ${task.cpus} -1 ${sample}_1.fastq.gz -2 ${sample}_2.fastq.gz - \
    > /dev/null 2>&1
    """
}

process ALIGN_FASTQ {
    tag "$sample"
    container "aakrosh/bwa-suite:alpine"  
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(fq1), path(fq2)
    
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} -Y -R \"@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina\" -K 100000000 ${params.reference} ${fq1} ${fq2} \
    | samblaster --addMateTags \
    | samtools view -b -h - \
    | samtools sort -@ ${task.cpus} -O BAM -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}

process GET_STATS {
    tag "${sample}"
    container "aakrosh/bwa-suite:alpine" 
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai) 

    output:
    tuple val(sample), path("${sample}.stats.txt") 

    script:
    """
    alignstats -i ${bam} -o ${sample}.stats.txt -p -W -P ${task.cpus} -t ${params.target_regions}
    """
}

process CALL_SMALL_VARIANTS {
    tag "joint calling"
    container 'community.wave.seqera.io/library/freebayes:1.3.9--8b4bf7c06d07bf77'

    input:
    path bams_and_bais 

    output:
    path "variants.vcf"

    script:
    def bam_files = bams_and_bais.findAll { it.name.endsWith('.bam') }
    def bam_files_str = bam_files.join(' ')

    """
    freebayes-parallel <(awk '{{print \$1":"\$2"-"\$3}}' ${params.target_regions}) ${task.cpus} --fasta-reference ${params.reference} ${bam_files_str}  \
    | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" -t "PASS" -F "FAIL" \
    > variants.vcf
    """
}

process BGZIP_INDEX_VARIANTS {
    tag "bgzip tabix"
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'
    publishDir 'results', mode: 'copy'

    input:
    path "variants.vcf"

    output:
    tuple path("variants.vcf.gz"), path("variants.vcf.gz.tbi")

    script:
    """
    bgzip -@ ${task.cpus} variants.vcf
    tabix variants.vcf.gz
    """
}

process PRIORITIZE_VARIANTS {
    tag "exomiser_prioritisation"
    container "openjdk:21-jdk"
    publishDir 'results/prioritized_variants', mode: 'copy'
    
    input:
    tuple path("variants.vcf.gz"), path("variants.vcf.gz.tbi") 
    path phenotype_file  // A file containing HPO terms
    path pedigree_file
    
    output:
    path "exomiser_results/*.html"
    path "exomiser_results/*.json"
    path "exomiser_results/*.tsv"
    path "exomiser_results/*.vcf"
    
    script:
    """
    # Create an Exomiser analysis YAML configuration file
    cat > analysis.yml << EOF
    ---
    analysis:
      genomeAssembly: ${params.genome_assembly}
      vcf: variants.vcf.gz
      ped: ${pedigree_file}
      proband: ${params.proband}
      hpoIds: [ \$(cat ${phenotype_file} | paste -sd, -) ]
      analysisMode: PASS_ONLY
      frequencySources: [THOUSAND_GENOMES, TOPMED, UK10K, GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_ASJ, GNOMAD_E_EAS, GNOMAD_E_FIN, GNOMAD_E_NFE, GNOMAD_E_OTH, GNOMAD_E_SAS, GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_ASJ, GNOMAD_G_EAS, GNOMAD_G_FIN, GNOMAD_G_NFE, GNOMAD_G_OTH, GNOMAD_G_SAS]
      pathogenicitySources: [REVEL, MVP, CADD]
      inheritanceModes: {
        AUTOSOMAL_DOMINANT: 0.1,
        AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
        AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
        X_DOMINANT: 0.1,
        X_RECESSIVE_HOM_ALT: 0.1,
        X_RECESSIVE_COMP_HET: 2.0,
        MITOCHONDRIAL: 0.2
      }
    steps: [
        variantEffectFilter: {
            remove: [
                FIVE_PRIME_UTR_EXON_VARIANT,
                FIVE_PRIME_UTR_INTRON_VARIANT,
                THREE_PRIME_UTR_EXON_VARIANT,
                THREE_PRIME_UTR_INTRON_VARIANT,
                NON_CODING_TRANSCRIPT_EXON_VARIANT,
                UPSTREAM_GENE_VARIANT,
                INTERGENIC_VARIANT,
                REGULATORY_REGION_VARIANT,
                CODING_TRANSCRIPT_INTRON_VARIANT,
                NON_CODING_TRANSCRIPT_INTRON_VARIANT,
                DOWNSTREAM_GENE_VARIANT
            ]
        },
        frequencyFilter: {maxFrequency: 2.0},
        pathogenicityFilter: {keepNonPathogenic: true},
        inheritanceFilter: {},
        omimPrioritiser: {},
        priorityScoreFilter: {priorityType: HIPHIVE_PRIORITY, minPriorityScore: 0.501},
        hiPhivePrioritiser: {},
    ]
    outputOptions:
      outputContributingVariantsOnly: true
      numGenes: 20
      outputPrefix: exomiser_results/prioritized
      outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
    EOF
    
    # Create the output directory
    mkdir -p exomiser_results
    
    # Run Exomiser CLI
    java -Xms2g -Xmx4g -jar ${params.exomiser_cli_jar} \
        --analysis analysis.yml \
        --spring.config.location=${params.exomiser_config}
    """
}

workflow {
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.cram_location)) }
        .set { input_channel }

    // extract the fastq files from the CRAM
    fq_channel = EXTRACT_FASTQ(input_channel)

    // align the reads into a BAM file
    bam_channel = ALIGN_FASTQ(fq_channel)

    // calculate the alignment stats
    stats_channel = GET_STATS(bam_channel)

    // call variants using FreeBayes
    bam_channel
        .map { sample, bam, bai -> [bam, bai] }
        .flatten()
        .collect()
        .set { all_bams_for_calling }
    variants_channel = CALL_SMALL_VARIANTS(all_bams_for_calling)

    // bgzip and index the variants
    zipped_variants_channel = BGZIP_INDEX_VARIANTS(variants_channel)

    // prioritize the variants using exomiser
    PRIORITIZE_VARIANTS(zipped_variants_channel,
                        file(params.phenotype_file),
                        file(params.pedigree))

}
