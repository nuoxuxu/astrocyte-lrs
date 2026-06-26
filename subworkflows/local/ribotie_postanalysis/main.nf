process plot_ribotie_figures {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/visualization"

    input:
    path mid_ribotie_merged, stageAs: 'mid_ribotie.csv'
    path mid_orfanage_gtf, stageAs: 'mid_orfanage.gtf'
    path mid_expression, stageAs: 'mid_expression.parquet'

    output:
    path "ribotie_proportion_by_expression.png"
    path "ribotie_expression_density.png"

    script:
    """
    plot_ribotie_figures.py \\
        --isoseq "Mid stringency" mid_ribotie.csv mid_ribotie.csv mid_orfanage.gtf mid_expression.parquet \\
        --output_dir .
    """
}

workflow RIBOTIE_POSTANALYSIS {
    take:
    ribotie_training_outputs
    orfanage_gtf
    final_expression

    main:
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .toList()
        .set { ribotie_entries }

    mid_ribotie_merged = ribotie_entries.map { entries ->
        file(entries[0].ribotie_merged_csv)
    }

    plot_ribotie_figures(mid_ribotie_merged, orfanage_gtf, final_expression)
}
