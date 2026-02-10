process RUN_RIBOTIE {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}"
    
    input:
    tuple val(param_set_name), path(gtf_h5), path(ribotie_h5), path(ribotie_yml)

    output:
    path("ribotie_res_*")
    path("models")
    path("multiqc")
    
    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotie $ribotie_yml
    """
}

workflow {
def manifest_data = new groovy.json.JsonSlurper().parse(file(params.manifest))

// Create channel from JSON
channel.fromList(manifest_data)
    .map { entry ->
        tuple(
            entry.param_set_name,
            file(entry.gtf_h5),
            file(entry.ribotie_h5),
            file(entry.ribotie_yml)
        )
    }
    .set { gpu_input_ch }    
    RUN_RIBOTIE(gpu_input_ch)
}