include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { LABEL_ORF_TYPE_GENCODE } from "./subworkflows/local/quality"

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }

    LABEL_ORF_TYPE_GENCODE(
        params.main_pipeline_outputs,
        params.ribotie_training_outputs,
        annotation_gtf,
        channel.value(params.orbl_alignment_set)
    )
}
