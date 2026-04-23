process phylocsfpp {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/mid_stringency"

    input:
    tuple val(condition), path(ribotie_output_gtf), path(phyloCSF_db)

    output:
    tuple val(condition), path("${ribotie_output_gtf.baseName}.PhyloCSF++.gtf")

    script:
    """
    phylocsf++ annotate-with-tracks $phyloCSF_db/PhyloCSF+1.bw $ribotie_output_gtf
    """
}

workflow GET_QUALITY_METRICS {
    take:
    ribotie_training_outputs
    PhyloCSFpp_db

    main:
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple("unstim", file(entry.ribotie_unstim_novel_gtf))
        }
        .set { unstim_input_ch }

    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple("stim", file(entry.ribotie_stim_novel_gtf))
        }
        .set { stim_input_ch }

    unstim_input_ch
        .mix(stim_input_ch)
        .combine(PhyloCSFpp_db)
        | phylocsfpp
}
