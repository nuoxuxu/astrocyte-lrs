process plot_ribotie_figures {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/visualization"

    input:
    path low_ribotie_stim, stageAs: 'low_ribotie_stim.csv'
    path low_ribotie_unstim, stageAs: 'low_ribotie_unstim.csv'
    path low_orfanage_gtf, stageAs: 'low_orfanage.gtf'
    path low_expression, stageAs: 'low_expression.parquet'
    path high_ribotie_stim, stageAs: 'high_ribotie_stim.csv'
    path high_ribotie_unstim, stageAs: 'high_ribotie_unstim.csv'
    path high_orfanage_gtf, stageAs: 'high_orfanage.gtf'
    path high_expression, stageAs: 'high_expression.parquet'
    path gencode_ribotie_stim, stageAs: 'gencode_ribotie_stim.csv'
    path gencode_ribotie_unstim, stageAs: 'gencode_ribotie_unstim.csv'
    path gencode_classification, stageAs: 'gencode_classification.parquet'
    path gencode_expression, stageAs: 'gencode_expression.parquet'

    output:
    path "ribotie_proportion_by_expression.png"
    path "ribotie_expression_density.png"

    script:
    """
    plot_ribotie_figures.py \\
        --isoseq "Low stringency" low_ribotie_stim.csv low_ribotie_unstim.csv low_orfanage.gtf low_expression.parquet \\
        --isoseq "High stringency" high_ribotie_stim.csv high_ribotie_unstim.csv high_orfanage.gtf high_expression.parquet \\
        --gencode GENCODE gencode_ribotie_stim.csv gencode_ribotie_unstim.csv gencode_classification.parquet gencode_expression.parquet \\
        --output_dir .
    """
}

workflow RIBOTIE_VISUALIZATION {
    take:
    ribotie_training_outputs
    orfanage_gtf
    final_expression
    final_classification

    main:
    // Parse ribotie JSON manifest to get CSV paths per param_set
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .toList()
        .set { ribotie_entries }

    // Convert upstream channels to value channels (lists) for multiple access
    orfanage_gtf.toList().set { orfanage_list }
    final_expression.toList().set { expression_list }
    final_classification.toList().set { classification_list }

    // Extract ribotie CSVs per param_set
    low_ribotie_stim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'low_stringency' }.ribotie_stim_csv)
    }
    low_ribotie_unstim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'low_stringency' }.ribotie_unstim_csv)
    }
    high_ribotie_stim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'high_stringency' }.ribotie_stim_csv)
    }
    high_ribotie_unstim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'high_stringency' }.ribotie_unstim_csv)
    }
    gencode_ribotie_stim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'gencode' }.ribotie_stim_csv)
    }
    gencode_ribotie_unstim = ribotie_entries.map { entries ->
        file(entries.find { entry -> entry.param_set_name == 'gencode' }.ribotie_unstim_csv)
    }

    // Extract orfanage GTFs per param_set
    low_orfanage = orfanage_list.map { list -> list.find { item -> item[0] == 'low_stringency' }[1] }
    high_orfanage = orfanage_list.map { list -> list.find { item -> item[0] == 'high_stringency' }[1] }

    // Extract expression parquets per param_set
    low_expression = expression_list.map { list -> list.find { item -> item[0] == 'low_stringency' }[1] }
    high_expression = expression_list.map { list -> list.find { item -> item[0] == 'high_stringency' }[1] }

    // Low stringency classification for GENCODE mapping
    low_classification = classification_list.map { list -> list.find { item -> item[0] == 'low_stringency' }[1] }

    plot_ribotie_figures(
        low_ribotie_stim, low_ribotie_unstim, low_orfanage, low_expression,
        high_ribotie_stim, high_ribotie_unstim, high_orfanage, high_expression,
        gencode_ribotie_stim, gencode_ribotie_unstim, low_classification, low_expression
    )
}
