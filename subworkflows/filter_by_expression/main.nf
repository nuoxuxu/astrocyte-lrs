process filter_by_expression {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_set_name}"

    input:
    tuple path(oarfish_quant_files), path(filtered_classification), path(filtered_gtf), path(corrected_fasta), val(param_set_name), val(min_reads), val(min_n_sample)    

    output:
    tuple val(param_set_name), path("final_classification.parquet"), emit: final_classification
    tuple val(param_set_name), path("final_transcripts.gtf"), emit: final_transcripts_gtf
    tuple val(param_set_name), path("final_expression.parquet"), emit: final_expression
    tuple val(param_set_name), path("final_transcripts.fasta"), emit: final_fasta

    script:
    """
    filter_by_expression.py \\
        --min_reads ${min_reads} \\
        --min_n_sample ${min_n_sample} \\
        --fasta $corrected_fasta
    """
}

process GffCompare {
    module "StdEnv/2023:gffcompare/0.12.6"
    label "short_slurm_job"
    storeDir "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_set_name}"

    input:
    path annotation_gtf
    tuple val(param_set_name), path(final_transcripts_gtf)

    output:
    tuple val(param_set_name), path("gffcmp.annotated.gtf"), emit: gffcmp_annotated_sgtf
    tuple val(param_set_name), path("gffcmp.loci"), emit: gffcmp_loci
    tuple val(param_set_name), path("gffcmp.stats"), emit: gffcmp_stats
    tuple val(param_set_name), path("gffcmp.tracking"), emit: gffcmp_tracking
    tuple val(param_set_name), path("gffcmp.final_transcripts.gtf.refmap"), emit: refmap
    tuple val(param_set_name), path("gffcmp.final_transcripts.gtf.tmap"), emit: tmap

    script:
    """
    gffcompare -r $annotation_gtf $final_transcripts_gtf
    """
}

process transcript_visualization {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/figures"

    input:
    tuple val(param_set_name), path(final_classification)
    tuple val(param_set_name), path(final_transcripts_gtf)
    tuple val(param_set_name), path(final_expression)
    tuple val(param_set_name), path(final_fasta)
    
    output:
    path("${param_set_name}_transcripts.pdf")

    script:
    """
    transcript_classification.R \\
        --classification $final_classification \\
        --expression $final_expression \\
        --output "${param_set_name}_transcripts.pdf"
    """    
}

workflow FILTER_BY_EXPRESSION {
    take:
    oarfish_quant
    filtered_classification
    filtered_gtf
    sqanti_corrected_fasta
    filter_configs
    annotation_gtf

    main:
    channel.fromList(filter_configs).set { isoform_exp_filter_params }

    filter_by_expression_input_ch = oarfish_quant
        .collect()
        .toList()
        .combine(filtered_classification)
        .combine(filtered_gtf)
        .combine(sqanti_corrected_fasta)
        .combine(isoform_exp_filter_params)
    filter_by_expression(filter_by_expression_input_ch)
    GffCompare(annotation_gtf, filter_by_expression.out.final_transcripts_gtf)

    emit:
    final_classification = filter_by_expression.out.final_classification
    final_transcripts_gtf = filter_by_expression.out.final_transcripts_gtf
    final_expression = filter_by_expression.out.final_expression
    final_fasta = filter_by_expression.out.final_fasta
    tmap = GffCompare.out.tmap
}
