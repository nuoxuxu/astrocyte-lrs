process plot_ribotie_figures {
    conda "/scratch/nxu/astrocytes/env"
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
    // Parse ribotie JSON manifest to get CSV paths per param_set
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .toList()
        .set { ribotie_entries }

    // Convert upstream channels to value channels (lists) for multiple access
    orfanage_gtf.toList().set { orfanage_list }
    final_expression.toList().set { expression_list }

    // Extract ribotie merged CSV for mid_stringency
    mid_ribotie_merged = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'mid_stringency' }.ribotie_merged_csv)
    }

    // Extract orfanage GTF and expression for mid_stringency
    mid_orfanage = orfanage_list.map { list -> list.find { item -> item[0] == 'mid_stringency' }[1] }
    mid_expression = expression_list.map { list -> list.find { item -> item[0] == 'mid_stringency' }[1] }

    plot_ribotie_figures(mid_ribotie_merged, mid_orfanage, mid_expression)
}
