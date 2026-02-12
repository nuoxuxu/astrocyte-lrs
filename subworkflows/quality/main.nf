process pfam_scan {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/${param_set_name}"

    input:
    tuple val(param_set_name), path(translation_fasta)
    path(pfamdb)
    
    output:
    tuple val(param_set_name), path("pfam_scan_results.csv")
    
    script:
    """
    pfam_scan.py \\
        -out pfam_scan_results.csv \\
        -outfmt csv \\
        $translation_fasta \\
        $pfamdb
    """
}

process phylocsfpp {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/${param_set_name}"

    input:
    tuple val(param_set_name), val(condition), path(ribotie_output_gtf), path(phyloCSF_db)

    output:
    tuple val(param_set_name), val(condition), path("${ribotie_output_gtf.baseName}.PhyloCSF++.gtf")

    script:
    """
    phylocsf++ annotate-with-tracks $phyloCSF_db/PhyloCSF+1.bw $ribotie_output_gtf
    """
}

workflow GET_QUALITY_METRICS {
    take:
    ribotie_training_outputs
    PhyloCSFpp_db
    translation_fasta
    pfamdb

    main:
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple(entry.param_set_name, "unstim", file(entry.ribotie_unstim_novel_gtf))
        }
        .set { unstim_input_ch }

    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple(entry.param_set_name, "stim", file(entry.ribotie_stim_novel_gtf))
        }
        .set { stim_input_ch }        

    unstim_input_ch
        .mix(stim_input_ch)
        .combine(PhyloCSFpp_db)
        | phylocsfpp

    // pfam_scan(translation_fasta, channel.fromPath(pfamdb))
}