process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}"

    input:
    tuple val(version), path(final_expression), path(predicted_cds_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    tuple val(version), path("IsoformSwitchAnalyzeR.rds"), emit: rds
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
    conda "/scratch/nxu/astrocytes/env"
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
    conda "/scratch/nxu/astrocytes/env"
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
    conda "/scratch/nxu/astrocytes/env"
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

    IsoseqsSwitchList.out.rds
        .multiMap { ver, rds ->
            versions: ver
            rdss: rds
        }
        .set { split_rds }

    PlotIsoformConsequences(split_rds.versions, split_rds.rdss)
    volcano_plot(split_rds.versions, split_rds.rdss)
    prepare_shinyApp(split_rds.versions, split_rds.rdss)

    emit:
    isoform_features_csv = IsoseqsSwitchList.out.isoform_features_csv
}

workflow ISOFORMSWITCH_MULTI {
    take:
    isoform_ch
    primer_to_sample

    main:

    isoform_ch
        .combine(primer_to_sample)
        .map { ver, expr, cds_gtf, classif, fasta, annot, primer ->
            tuple(ver, expr, cds_gtf, classif, primer, fasta, annot)
        }
    | IsoseqsSwitchList

    IsoseqsSwitchList.out.rds
        .multiMap { ver, rds ->
            versions: ver
            rdss: rds
        }
        .set { split_rds }

    PlotIsoformConsequences(split_rds.versions, split_rds.rdss)
    volcano_plot(split_rds.versions, split_rds.rdss)
    prepare_shinyApp(split_rds.versions, split_rds.rdss)

    emit:
    isoform_features_csv = IsoseqsSwitchList.out.isoform_features_csv
}
