include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { LABEL_ORF_TYPE_GENCODE } from "./subworkflows/local/quality"
include { ISOFORMSWITCH_MULTI } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"

process filter_ribotie_for_isoformswitch {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/filtered_for_isoformswitch"

    input:
    tuple val(name), path(input_gtf), path(input_fasta), path(input_expression), path(input_classification)

    output:
    tuple val(name), path("${name}_filtered.fasta"),                            emit: fasta
    tuple val(name), path("${name}_filtered_expression.parquet"),               emit: expression
    tuple val(name), path("${name}_filtered_classification.parquet"),           emit: classification

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --input_fasta $input_fasta \\
        --input_expression $input_expression \\
        --input_classification $input_classification \\
        --output_fasta ${name}_filtered.fasta \\
        --output_expression ${name}_filtered_expression.parquet \\
        --output_classification ${name}_filtered_classification.parquet
    """
}

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }

    LABEL_ORF_TYPE_GENCODE(
        params.main_pipeline_outputs,
        params.ribotie_training_outputs,
        annotation_gtf,
        channel.value(params.orbl_alignment_set)
    )

    // Build per-version channel from manifests (exclude gencode)
    channel.fromPath(params.main_pipeline_outputs)
        .map { f ->
            def name = f.baseName.replaceAll(/^main_pipeline_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.final_expression), file(entry.final_fasta), file(entry.final_classification))
        }
        .set { main_outputs_ch }

    channel.fromPath(params.ribotie_training_outputs)
        .filter { f -> f.baseName.replaceAll(/^ribotie_training_outputs_/, '') != 'gencode' }
        .map { f ->
            def name = f.baseName.replaceAll(/^ribotie_training_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.ribotie_merged_gtf))
        }
        .join(main_outputs_ch)
        .map { name, ribotie_gtf, expr, fasta, classif ->
            tuple(name, ribotie_gtf, fasta, expr, classif)
        }
        .set { ribotie_inputs_ch }

    // Filter RiboTIE GTF to match expression data (versioned)
    filter_ribotie_for_isoformswitch(ribotie_inputs_ch)

    // Extract GTF separately (needed for ISOFORMSWITCH)
    channel.fromPath(params.ribotie_training_outputs)
        .filter { f -> f.baseName.replaceAll(/^ribotie_training_outputs_/, '') != 'gencode' }
        .map { f ->
            def name = f.baseName.replaceAll(/^ribotie_training_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.ribotie_merged_gtf))
        }
        .set { gtf_ch }

    // Build per-version isoform_ch for ISOFORMSWITCH using filtered outputs
    gtf_ch
        .join(filter_ribotie_for_isoformswitch.out.expression)
        .join(filter_ribotie_for_isoformswitch.out.fasta)
        .join(filter_ribotie_for_isoformswitch.out.classification)
        .combine(annotation_gtf)
        .map { name, gtf, expr, fasta, classif, annot ->
            tuple(name, expr, gtf, classif, fasta, annot)
        }
        .set { isoform_ch }

    ISOFORMSWITCH_MULTI(isoform_ch, primer_to_sample)
}
