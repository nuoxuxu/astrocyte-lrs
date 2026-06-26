process FILTER_RIBOTIE_RESULTS {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/filtered"

    input:
    tuple path(input_gtf), path(input_fasta), path(input_expression), path(input_classification)

    output:
    path "filtered_RiboTIE.fasta",                            emit: filtered_RiboTIE_fasta
    path "ribotie_filtered_final_expression.parquet",         emit: filtered_RiboTIE_expression
    path "ribotie_filtered_final_classification.parquet",     emit: filtered_RiboTIE_classification

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --input_fasta $input_fasta \\
        --input_expression $input_expression \\
        --input_classification $input_classification \\
        --output_fasta filtered_RiboTIE.fasta \\
        --output_expression ribotie_filtered_final_expression.parquet \\
        --output_classification ribotie_filtered_final_classification.parquet
    """
}

process translateORFs {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/filtered"

    input:
    path ref_genome_fasta
    path filtered_RiboTIE_gtf

    output:
    path("filtered_RiboTIE_proteins.fasta"), emit: filtered_RiboTIE_proteins

    script:
    """
    gffread -y filtered_RiboTIE_proteins.fasta \\
        -g $ref_genome_fasta \\
        $filtered_RiboTIE_gtf
    """
}

workflow FILTER_RIBOTIE {
    take:
    ribotie_filtered_gtf
    final_fasta
    final_expression
    final_classification
    ref_genome_fasta

    main:
    ribotie_filtered_gtf
        .combine(final_fasta)
        .combine(final_expression)
        .combine(final_classification)
        .set { joined }
    FILTER_RIBOTIE_RESULTS(joined)

    translateORFs(ref_genome_fasta, ribotie_filtered_gtf)

    emit:
    filtered_RiboTIE_fasta          = FILTER_RIBOTIE_RESULTS.out.filtered_RiboTIE_fasta
    filtered_RiboTIE_expression     = FILTER_RIBOTIE_RESULTS.out.filtered_RiboTIE_expression
    filtered_RiboTIE_classification = FILTER_RIBOTIE_RESULTS.out.filtered_RiboTIE_classification
    filtered_RiboTIE_proteins       = translateORFs.out.filtered_RiboTIE_proteins
}
