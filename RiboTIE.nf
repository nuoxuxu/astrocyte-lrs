process RUN_RIBOTIE {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "nextflow_results/ribotie/mid_stringency/${mode}"

    input:
    tuple path(gtf_h5), path(ribotie_h5), val(mode), path(ribotie_yml)

    output:
    path("ribotie_res_${mode}.redundant.gtf"), emit: redundant_gtf
    path("ribotie_res_${mode}.redundant.csv"), emit: redundant_csv
    path("ribotie_res_${mode}.novel.gtf"),     emit: novel_gtf
    path("ribotie_res_${mode}.novel.csv"),     emit: novel_csv
    path("ribotie_res_${mode}.gtf"),           emit: filtered_gtf
    path("ribotie_res_${mode}.csv"),           emit: filtered_csv
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
    storeDir "export/local"

    input:
    path(input_gtf)

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
                file(entry.gtf_h5),
                file(entry.ribotie_h5),
                entry.mode,
                file(entry.ribotie_yml)
            )
        }
        .set { gpu_input_ch }

    RUN_RIBOTIE(gpu_input_ch)
    save_to_parquet(RUN_RIBOTIE.out.redundant_gtf)
}
