process RUN_RIBOTIE {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}"
    
    input:
    tuple val(param_set_name), path(gtf_h5), path(ribotie_h5), path(ribotie_yml)

    output:
    tuple val(param_set_name), path("ribotie_res_*"), path("models"), path("multiqc")
    
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
                file(entry.ribotie_yml)
            )
        }
        .set { gpu_input_ch }
        RUN_RIBOTIE(gpu_input_ch)
    
    RUN_RIBOTIE
        .out
        .toList()
        .map { results ->
            def outputs = results.collect { param_set_name, _ribotie_res_, _models, _multiqc ->
                [
                    param_set_name: param_set_name,
                    ribotie_unstim_novel_gtf: "nextflow_results/ribotie/${param_set_name}/ribotie_res_Unstim.novel.gtf",
                    ribotie_unstim_novel_csv: "nextflow_results/ribotie/${param_set_name}/ribotie_res_Unstim.novel.csv",
                    ribotie_stim_novel_gtf: "nextflow_results/ribotie/${param_set_name}/ribotie_res_Stim.novel.gtf",
                    ribotie_stim_novel_csv: "nextflow_results/ribotie/${param_set_name}/ribotie_res_Stim.novel.csv",
                ]
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'ribotie_training_outputs.json', storeDir: 'nextflow_results/manifests')
}