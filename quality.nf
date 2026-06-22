include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { LABEL_ORF_TYPE_GENCODE } from "./subworkflows/local/quality"
include { ISOFORMSWITCH_MULTI } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }

    LABEL_ORF_TYPE_GENCODE(
        params.main_pipeline_outputs,
        params.ribotie_training_outputs,
        annotation_gtf,
        channel.value(params.orbl_alignment_set)
    )

    // Build per-version isoform_ch from manifests (exclude gencode)
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
        .combine(annotation_gtf)
        .map { name, gtf, expr, fasta, classif, annot ->
            tuple(name, expr, gtf, classif, fasta, annot)
        }
        .set { isoform_ch }

    ISOFORMSWITCH_MULTI(isoform_ch, primer_to_sample)
}
