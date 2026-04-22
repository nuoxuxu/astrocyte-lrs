process FILTER_RIBOTIE_RESULTS {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}/filtered"

    input:
    tuple val(param_set_name), path(input_gtf)
    path ribotie_cpm1_3sample

    output:
    tuple val(param_set_name), path("filtered_RiboTIE.gtf")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --ribotie_cpm1_3sample $ribotie_cpm1_3sample \\
        --output_gtf filtered_RiboTIE.gtf
    """
}

workflow FILTER_RIBOTIE {
    take:
    ribotie_filtered_gtf  // channel: [param_set_name, gtf]
    ribotie_cpm1_3sample  // value channel: path

    main:
    FILTER_RIBOTIE_RESULTS(ribotie_filtered_gtf, ribotie_cpm1_3sample)

    emit:
    filtered_gtf = FILTER_RIBOTIE_RESULTS.out
}
