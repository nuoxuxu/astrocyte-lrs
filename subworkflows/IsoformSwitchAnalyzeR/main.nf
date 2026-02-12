process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${param_set_name}"

    input:
    tuple val(param_set_name), path(final_expression), path(orfanage_gtf), path(final_classification), path(primer_to_sample), path(corrected_fasta), path(annotation_gtf)

    output:
    path("IsoformSwitchAnalyzeR.rds")

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --corrected_fasta $corrected_fasta \\
        --orfanage_gtf $orfanage_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification
    """
}

workflow ISOFORMSWITCH {
    take:
    final_expression
    primer_to_sample
    corrected_fasta
    orfanage_gtf
    annotation_gtf
    final_classification

    main:
    
    final_expression
        .join(orfanage_gtf)
        .join(final_classification)
        .combine(primer_to_sample)
        .combine(corrected_fasta)
        .combine(annotation_gtf)
    | IsoseqsSwitchList
    // IsoseqsSwitchList(final_expression, primer_to_sample, corrected_fasta, orfanage_gtf, annotation_gtf, final_classification)
}