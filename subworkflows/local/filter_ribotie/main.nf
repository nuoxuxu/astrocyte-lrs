process FILTER_RIBOTIE_RESULTS {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}/filtered"

    input:
    tuple val(param_set_name), path(input_gtf), path(input_fasta), path(input_expression)
    path ribotie_cpm1_3sample

    output:
    tuple val(param_set_name), path("filtered_RiboTIE.gtf"), path("filtered_RiboTIE.fasta"), path("ribotie_filtered_final_expression.parquet")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --ribotie_cpm1_3sample $ribotie_cpm1_3sample \\
        --input_fasta $input_fasta \\
        --input_expression $input_expression \\
        --output_gtf filtered_RiboTIE.gtf \\
        --output_fasta filtered_RiboTIE.fasta \\
        --output_expression ribotie_filtered_final_expression.parquet
    """
}

workflow FILTER_RIBOTIE {
    take:
    ribotie_filtered_gtf  // channel: [param_set_name, gtf]
    ribotie_cpm1_3sample  // value channel: path
    final_fasta           // channel: [param_set_name, fasta]
    final_expression      // channel: [param_set_name, expression]

    main:
    ribotie_filtered_gtf
        .join(final_fasta)
        .join(final_expression)
        .set { joined }
    FILTER_RIBOTIE_RESULTS(joined, ribotie_cpm1_3sample)

    emit:
    filtered = FILTER_RIBOTIE_RESULTS.out
}
