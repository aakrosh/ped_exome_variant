// Description: Run a workflow that aligns a CRAM file using BWA-MEM.
nextflow.enable.dsl = 2
include { bgzip_index_variants as bgzip_index_raw; 
          bgzip_index_variants as bgzip_index_standard;
          bgzip_index_variants as bgzip_index_vep } from "./modules/index"

process extract_fastq {
    tag "$sample"
    container 'aakrosh/bwa-suite:alpine'

    input:
    tuple val(sample), path(cram)

    output:
    tuple val(sample),path("${sample}_1.fastq.gz"),path("${sample}_2.fastq.gz") 
    
    script:
    """
    samtools sort -n -@ ${task.cpus} ${cram} \
    | samtools fastq -@ ${task.cpus} \
      -1 ${sample}_1.fastq.gz -2 ${sample}_2.fastq.gz - \
    > /dev/null 2>&1
    """
}

process align_fastq {
    tag "$sample"
    container "aakrosh/bwa-suite:alpine"  
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(fq1), path(fq2)
    
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} -Y \
      -R \"@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina\" \
      -K 100000000 ${params.reference} ${fq1} ${fq2} \
    | samblaster --addMateTags \
    | samtools view -b -h - \
    | samtools sort -@ ${task.cpus} -O BAM -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}

process get_stats {
    tag "${sample}"
    container "aakrosh/bwa-suite:alpine" 
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai) 

    output:
    tuple val(sample), path("${sample}.stats.txt") 

    script:
    """
    alignstats -i ${bam} -o ${sample}.stats.txt -p -W -P ${task.cpus} \
        -t ${params.target_regions}
    """
}

process call_small_variants {
    tag "joint calling"
    container 'community.wave.seqera.io/library/freebayes:1.3.9--8b4bf7c06d07bf77'

    input:
    path(bams_and_bais)

    output:
    path("variants.vcf")

    script:
    def bam_files = bams_and_bais.findAll { it.name.endsWith('.bam') }
    def bam_files_str = bam_files.join(' ')

    """
    freebayes-parallel <(awk '{{print \$1":"\$2"-"\$3}}' ${params.target_regions}) \
        ${task.cpus} --genotype-qualities \
        --fasta-reference ${params.reference} ${bam_files_str}  \
    > variants.vcf
    """
}

process filter_variants {
    tag "filter"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"

    input:
    path(vcf)
    path(exclude_bed)

    output:
    path("${vcf.simpleName}.flt.vcf")

    script:
    """
    bcftools filter -i 'QUAL>1 && QUAL/INFO/AO>10 && INFO/SAF>0 && INFO/SAR>0 && INFO/RPR>1 && INFO/RPL>1' ${vcf} \
    | bcftools +setGT -- -t q -n . -i 'FMT/GQ<20' \
    | bcftools annotate -a ${exclude_bed} -c CHROM,FROM,TO \
        -h <(echo '##FILTER=<ID=BED_OVERLAP,Description="Variant overlaps BED regions">') \
        -m BED_OVERLAP -Ov -o ${vcf.simpleName}.flt.vcf
    """
}

process calculate_variant_stats {
    tag "variant stats"
    container "openjdk:21-jdk"
    publishDir 'results', mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path(pedigree_file)

    output:
    path("variant_qc.html")

    script:
    """
    java -jar ${params.variantqc_jar} VariantQC -R ${params.reference} \
        -ped ${pedigree_file} -V ${vcf} -O variant_qc.html
    """
}

process leftalign_split {
    tag "standardize"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
    publishDir 'results', mode: 'copy'    

    input:
    tuple path(vcf), path(tbi)

    output:
    path("${vcf.baseName}.norm.vcf")

    script:
    """
    bcftools norm -a -d all -f ${params.reference} -m -both -Ov ${vcf} \
    | bcftools view --trim-alt-alleles -Ov -o ${vcf.baseName}.norm.vcf
    """
}

process prioritize_variants {
    tag "exomiser_prioritisation"
    container "openjdk:21-jdk"
    publishDir 'results', mode: 'copy'
    
    input:
    tuple path("variants.norm.vcf.gz"), path("variants.norm.vcf.gz.tbi") 
    path(phenotype_file)  // A file containing HPO terms
    path(pedigree_file)
    
    output:
    path("exomiser_results/prioritized.html")
    path("exomiser_results/prioritized.json")
    path("exomiser_results/prioritized.genes.tsv")
    path("exomiser_results/prioritized.variants.tsv")
    path("exomiser_results/prioritized.vcf.gz")
    path("exomiser_results/prioritized.vcf.gz.tbi")
    
    script:
    """
    # Create an Exomiser analysis YAML configuration file
    cat > analysis.yml << EOF
    analysis:
        genomeAssembly: ${params.genome_assembly}
        vcf: variants.norm.vcf.gz
        ped: ${pedigree_file}
        proband: ${params.proband}
        hpoIds: [ \$(cat ${phenotype_file} | paste -sd, -) ]
        analysisMode: PASS_ONLY
        frequencySources: [
            THOUSAND_GENOMES, 
            TOPMED, 
            UK10K, 
            GNOMAD_E_AFR, 
            GNOMAD_E_AMR, 
            GNOMAD_E_ASJ, 
            GNOMAD_E_EAS, 
            GNOMAD_E_FIN, 
            GNOMAD_E_NFE, 
            GNOMAD_E_OTH, 
            GNOMAD_E_SAS, 
            GNOMAD_G_AFR, 
            GNOMAD_G_AMR, 
            GNOMAD_G_ASJ, 
            GNOMAD_G_EAS, 
            GNOMAD_G_FIN, 
            GNOMAD_G_NFE, 
            GNOMAD_G_OTH, 
            GNOMAD_G_SAS]
        pathogenicitySources: [REVEL, MVP, SPLICE_AI, ALPHA_MISSENSE]
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
            hiPhivePrioritiser: {runParams: 'human'},
        ]
    outputOptions:
        outputContributingVariantsOnly: true
        numGenes: 30
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

process partition_variants {
    tag 'partition'
    container 'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11'

    input:
    tuple path(vcf), path(tbi)

    output:
    path('*.vcf')

    script:
    """
    for sample in \$(bcftools query -l ${vcf}); do
        bcftools view -c1 -s "\${sample}" -Ou --threads ${task.cpus} ${vcf} \
        | bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ -Ov \
          -o "\${sample}.tmp.vcf" --threads ${task.cpus} 

        awk 'BEGIN {OFS="\\t"}
        /^#/ {print; next}
        {
         split(\$9, fmt, ":");
         split(\$10, dat, ":");
         for (i in fmt) map[fmt[i]] = dat[i];
         new_fmt = "GT:AD:DP:GQ";
         \$9 = new_fmt; 
         \$10 = map["GT"]":"map["AD"]":"map["DP"]":"map["GQ"];
         print
        }' "\${sample}.tmp.vcf" > "\${sample}.vcf"

        rm "\${sample}.tmp.vcf"
    done
    """
}

process annotate_vep {
    tag "VEP"
    container 'ensemblorg/ensembl-vep:release_114.0'

    input:
    path(vcf)
   
    output:
    tuple val("${vcf.baseName}"), val("VEP"), path("${vcf.simpleName}.vep.vcf")

    script:
    """
    vep --offline --cache --fork ${task.cpus} \
    --dir_cache ${params.vep_cachedir} --fasta ${params.reference} \
    --use_given_ref --species homo_sapiens --assembly ${params.genome_assembly} \
    --xref_refseq --hgvs --hgvsg --canonical --symbol --distance 0 \
    --exclude_predicted --flag_pick --lookup_ref --force --input_file ${vcf} \
    --output_file ${vcf.simpleName}.vep.vcf --format vcf --vcf --no_stats \
    --numbers
    """
}

process annotate_annovar {
    tag "Annovar"

    input:
    tuple val(sample), val(type), path(vcf)

    output:
    tuple val(sample), val("ANNOVAR"), path("${sample}.*_multianno.txt")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    perl ${params.annovar_dir}/table_annovar.pl ${vcf} \
        ${params.annovar_database} --buildver ${suffix} --out ${sample} \
        --remove --protocol gnomad211_exome,gnomad211_genome --operation f,f \
        --vcfinput --thread ${task.cpus}
    """
}

process annotate_intervar {
    tag "Intervar"
    container 'aakrosh/intervar:latest'

    input:
    tuple val(sample), val(type), path(vcf)

    output:
    tuple val(sample), val("INTERVAR"), path("${vcf.simpleName}.*_multianno.txt.intervar")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    /usr/bin/python3 /Intervar.py -b ${suffix} -i ${vcf} --input_type=VCF \
        -o ${vcf.simpleName} -t ${params.intervar_database} \
        --table_annovar ${params.annovar_dir}/table_annovar.pl \
        --database_locat ${params.annovar_database} \
        --convert2annovar ${params.annovar_dir}/convert2annovar.pl \
        --annotate_variation ${params.annovar_dir}/annotate_variation.pl
    """
}

process annotate_autopvs1 {
    tag "autopvs1"
    container 'aakrosh/autopvs1:latest'    

    input:
    tuple val(sample), val(type), path(vcf)
    path(autopvs1_data)
    path(autopvs1_config)

    output:
    tuple val(sample), val("AUTOPVS1"), path("${vcf.simpleName}.autopvs1.txt")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    /usr/bin/python3 /autopvs1/autoPVS1_from_VEP_vcf.py \
        --genome_version ${suffix} \
        --vep_vcf ${vcf} > ${vcf.simpleName}.autopvs1.txt
    """
}

process select_clinvar {
    tag "clinvar"
    container 'pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.3'    
    containerOptions = { "--bind ${params.autogvp_dir}:/home/rstudio/AutoGVP" }

    output:
    path("ClinVar-selected-submissions.tsv")

    script:
    """
    Rscript /home/rstudio/AutoGVP/scripts/select-clinVar-submissions.R \
        --variant_summary /home/rstudio/AutoGVP/data/variant_summary.txt.gz \
        --submission_summary /home/rstudio/AutoGVP/data/submission_summary.txt.gz \
        --outdir . \
        --conceptID_list /home/rstudio/AutoGVP/data/clinvar_cpg_concept_ids.txt \
        --conflict_res "latest"
    """
}

process run_autogvp {
    tag "autogvp"
    container 'pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.3'    
    containerOptions = { "--bind ${params.autogvp_dir}:/home/rstudio/AutoGVP" }
    publishDir 'results', mode: 'copy'
   
    input:
    tuple val(sample), path(vep_file), path(annovar_file), 
            path(intervar_file), path(autopvs1_file)
    path(clinvar_file)

    output:
    path("${vep_file.baseName}*")

    script:
    """
    bash /home/rstudio/AutoGVP/run_autogvp.sh --workflow="custom" \
    --vcf=${vep_file} \
    --clinvar=/home/rstudio/AutoGVP/data/clinvar.vcf.gz \
    --intervar=${intervar_file} \
    --multianno=${annovar_file} \
    --autopvs1=${autopvs1_file} \
    --outdir=. \
    --out="${vep_file.baseName}" \
    --selected_clinvar_submissions=${clinvar_file} \
    --variant_summary=/home/rstudio/AutoGVP/data/variant_summary.txt.gz \
    --submission_summary=/home/rstudio/AutoGVP/data/submission_summary.txt.gz \
    --conceptIDs=/home/rstudio/AutoGVP/data/clinvar_cpg_concept_ids.txt \
    --conflict_res="latest"   
    """
}

workflow {
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.cram_location)) }
        .set { input_channel }

    // extract the fastq files from the CRAM
    extract_fastq(input_channel)
        .set { fq_channel }

    // align the reads into a BAM file
    align_fastq(fq_channel)
        .set { bam_channel }

    // calculate the alignment stats
    get_stats(bam_channel)
        .set { stats_channel }

    // call variants using FreeBayes
    bam_channel
        .map { sample, bam, bai -> [bam, bai] }
        .flatten()
        .collect()
        .set { all_bams_for_calling }
    call_small_variants(all_bams_for_calling)
        .set { variants_channel }

    // filter the variants to remove low-quality variants and genotypes
    filter_variants(variants_channel, file(params.exclusion_bed))
        .set { filtered_variants_channel }

    // bgzip and index the variants
    bgzip_index_raw(filtered_variants_channel)
        .set { zipped_variants_channel }

    // collect some stats
    calculate_variant_stats(zipped_variants_channel, file(params.pedigree))

    // left-align and normalize indels, split multiallelic sites
    leftalign_split(zipped_variants_channel)
        .set { standardize_channel }

    // bzgip and index the variants
    bgzip_index_standard(standardize_channel)
        .set { zipped_standard_channel }

    // prioritize the variants using exomiser
    prioritize_variants(zipped_standard_channel,
                        file(params.phenotype_file),
                        file(params.pedigree))

    // here on, we want to calculate variant effect per sample in the VCF file
    // applying ACMG guidelines. First partition the multi-vcf into per-sample
    // VCF. Also, AUTOGVP only expects GT,AD,DP,GQ for sample
    partition_variants(zipped_standard_channel)
        .set { per_sample_vcfs_ch }

    // annotate the file using VEP
    per_sample_vcfs_ch
        .flatten()
        .set { single_vcf_ch }

    annotate_vep(single_vcf_ch)
        .set { vep_annotated_channel }

    // annotate the variants using ANNOVAR
    annotate_annovar(vep_annotated_channel)    
        .set { annovar_output_channel }

    // annotate the variants using INTERVAR
    annotate_intervar(vep_annotated_channel)
        .set { intervar_output_channel }

    // annotate the variants using AUTOPVS1
    annotate_autopvs1(vep_annotated_channel, 
                      file(params.autopvs1_data), 
                      file(params.autopvs1_config))
        .set { autopvs1_output_channel }

    // create the select clinvar file
    select_clinvar()
        .set { clinvar_channel }

    // run the AutoGVP workflow
    all_outputs = vep_annotated_channel
        .mix(annovar_output_channel)
        .mix(intervar_output_channel)
        .mix(autopvs1_output_channel)
        .map { sample, tool, file -> tuple(sample, tuple(tool, file)) }

    grouped = all_outputs
        .groupTuple()
        .map { sample, records ->
            def file_map = records.collectEntries { [(it[0]): it[1]] }
            def vep      = file_map.get('VEP')
            def annovar  = file_map.get('ANNOVAR')
            def intervar = file_map.get('INTERVAR')
            def autopvs1 = file_map.get('AUTOPVS1')

            if ([vep, annovar, intervar, autopvs1].any { it == null }) return null
            return tuple(sample, vep, annovar, intervar, autopvs1)
        }
        .filter { it != null }

    run_autogvp(grouped, clinvar_channel)

}
