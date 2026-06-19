process RUN_RIBOTIE {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "nextflow_results/ribotie/${input_gtf_name}/${merged_or_separate}"

    input:
    tuple val(input_gtf_name), path(gtf_h5), path(ribotie_h5), val(merged_or_separate), path(ribotie_yml)

    output:
    tuple val(input_gtf_name), path("ribotie_res_${merged_or_separate}.redundant.gtf"),                                              emit: redundant_gtf
    path("ribotie_res_${merged_or_separate}.redundant.csv"),                                                                        emit: redundant_csv
    path("ribotie_res_${merged_or_separate}.novel.gtf"),                                                                            emit: novel_gtf
    path("ribotie_res_${merged_or_separate}.novel.csv"),                                                                            emit: novel_csv
    path("ribotie_res_${merged_or_separate}.gtf"),                                                                                  emit: filtered_gtf
    path("ribotie_res_${merged_or_separate}.csv"),                                                                                  emit: filtered_csv
    tuple path("ribotie_res_*.ckpt"), path("ribotie_res_*.npy"), path("ribotie_res_*.yml"), path("models"), path("multiqc"), emit: other

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotie $ribotie_yml
    """
}

process save_to_parquet {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "export/ribotie/${input_gtf_name}"

    input:
    tuple val(input_gtf_name), path(input_gtf)

    output:
    path("${input_gtf.baseName}.parquet")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    convert_gtf_to_parquet.py --input_gtf $input_gtf --output_gtf "${input_gtf.baseName}.parquet"
    """
}

workflow {
    channel.fromPath(params.ribotie_training_inputs)
        .splitJson()
        .map { entry ->
            tuple(
                entry.input_gtf_name,
                file(entry.gtf_h5),
                file(entry.ribotie_h5),
                entry.merged_or_separate,
                file(entry.ribotie_yml)
            )
        }
        .set { gpu_input_ch }

    RUN_RIBOTIE(gpu_input_ch)
    // save_to_parquet(RUN_RIBOTIE.out.redundant_gtf)
}
