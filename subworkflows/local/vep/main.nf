process prepare_gtf {
    label "short_slurm_job"
    module "gcc/12.3:StdEnv/2023:ngstools/1.0.1"
    storeDir "nextflow_results/vep/prepare_gtf"

    input:
    path annotation_gtf

    output:
    path "${annotation_gtf}.gz", emit: gtf_gz
    path "${annotation_gtf}.gz.tbi", emit: gtf_gz_tbi

    script:
    """
    grep -v "#" $annotation_gtf | sort -k1,1 -k4,4n -k5,5n -t\$'\t' | bgzip -c > "${annotation_gtf}.gz"
    tabix -p gff ${annotation_gtf}.gz
    """
}

process vep_utr_annotator {
    label "short_slurm_job"
    container "ensembl-vep.sif"
    storeDir "nextflow_results/vep/${sample_id}"

    input:
    tuple val(sample_id), path(vcf_file)
    path uorf_data_file
    path genome_fasta
    path gtf_gz
    path gtf_gz_tbi

    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf"), emit: annotated_vcf
    tuple val(sample_id), path("${sample_id}.vep.vcf_summary.html"), emit: summary

    script:
    """
    vep \\
        --fasta $genome_fasta \\
        --input_file $vcf_file \\
        --output_file ${sample_id}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --gtf $gtf_gz \\
        --plugin UTRAnnotator,file=$uorf_data_file \\
        --fork ${task.cpus} \\
        --force_overwrite
    """
}

workflow RUN_VEP {
    take:
    vcf_ch        // channel of [sample_id, vcf_file]
    uorf_data_file
    genome_fasta
    annotation_gtf
    main:
    prepare_gtf(annotation_gtf)
    vep_utr_annotator(vcf_ch, uorf_data_file, genome_fasta, prepare_gtf.out.gtf_gz, prepare_gtf.out.gtf_gz_tbi)

    emit:
    annotated_vcf = vep_utr_annotator.out.annotated_vcf
    summary       = vep_utr_annotator.out.summary
}
