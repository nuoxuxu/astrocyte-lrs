process RUN_RIBOTIE {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "${base_dir}"
    
    input:
    tuple val(param_set_name), path(gtf_h5), path(ribotie_h5), val(base_dir), path(ribotie_yml)

    output:
    tuple val(param_set_name), val(base_dir), path("ribotie_res_*"), path("models"), path("multiqc")
    
    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotie $ribotie_yml
    """
}

workflow {
    channel.fromPath(params.ribotie_training_inputs)
        .splitJson()
        .map { entry ->
            tuple(
                entry.param_set_name,
                file(entry.gtf_h5),
                file(entry.ribotie_h5),
                file(entry.base_dir),
                file(entry.ribotie_yml)
            )
        }
        .set { gpu_input_ch }
        
    RUN_RIBOTIE(gpu_input_ch)
}