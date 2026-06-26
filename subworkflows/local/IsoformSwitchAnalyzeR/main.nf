process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}"

    input:
    tuple val(version), path(final_expression), path(predicted_cds_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    path("IsoformSwitchAnalyzeR.rds"), emit: rds
    path("isoformFeatures.csv"), emit: isoform_features_csv

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --final_fasta $final_fasta \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification \\
    """
}

process PlotIsoformConsequences {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/figures"

    input:
    val(version)
    path(rds)

    output:
    path("splicing_consequences.pdf")
    path("switch_consequences.pdf")

    script:
    """
    plot_isoform_consequences.R --rds $rds
    """
}

process volcano_plot {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/figures"

    input:
    val(version)
    path(rds)

    output:
    path("volcano_DGE.pdf")
    path("volcano_DTU.pdf")
    path("volcano_iso_DGE.pdf")

    script:
    """
    IsoformSwitchAnalyzeR_volcano.R --switchlist $rds
    """
}

process prepare_shinyApp {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "astrocyte_vis_app/data/${version}"

    input:
    val(version)
    path(rds)

    output:
    path("switchPlotFromTables.RData")

    script:
    """
    extract_rds_data.R --rds $rds --output switchPlotFromTables.RData
    """
}

workflow ISOFORMSWITCH {
    take:
    version
    final_expression
    primer_to_sample
    final_fasta
    predicted_cds_gtf
    annotation_gtf
    final_classification

    main:

    final_expression
        .combine(predicted_cds_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(final_fasta)
        .combine(annotation_gtf)
        .combine(version)
        .map { expr, cds_gtf, classif, primer, fasta, annot, ver ->
            tuple(ver, expr, cds_gtf, classif, primer, fasta, annot)
        }
    | IsoseqsSwitchList

    PlotIsoformConsequences(version, IsoseqsSwitchList.out.rds)
    volcano_plot(version, IsoseqsSwitchList.out.rds)
    prepare_shinyApp(version, IsoseqsSwitchList.out.rds)

    emit:
    isoform_features_csv = IsoseqsSwitchList.out.isoform_features_csv
}
