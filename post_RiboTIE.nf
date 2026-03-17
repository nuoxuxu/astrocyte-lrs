include { GET_QUALITY_METRICS } from "./subworkflows/quality"
include { ISOFORMSWITCH } from "./subworkflows/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/ribotie_postanalysis/main.nf"

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }
    channel.value(file(params.pfamdb)).set { pfamdb }
    channel.value(file(params.PhyloCSFpp_db)).set { PhyloCSFpp_db }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> tuple(entry.param_set_name, file(entry.final_expression)) }
        .set { final_expression }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> tuple(entry.param_set_name, file(entry.final_fasta)) }
        .set { final_fasta }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> tuple(entry.param_set_name, file(entry.final_classification)) }
        .set { final_classification }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> tuple(entry.param_set_name, file(entry.orfanage_gtf)) }
        .set { orfanage_gtf }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> tuple(entry.param_set_name, file(entry.orfanage_proteins)) }
        .set { orfanage_proteins }

    channel.fromPath(params.main_pipeline_outputs)
        .splitJson()
        .map { entry -> entry.ribotie_output_gtf }
        .set { ribotie_output_gtf }
    
    ribotie_output_gtf.view()

    ISOFORMSWITCH(final_expression, primer_to_sample, final_fasta, orfanage_gtf, annotation_gtf, final_classification, orfanage_proteins, pfamdb, file(params.Human_coding_transcripts_CDS), file(params.Human_noncoding_transcripts_RNA), file(params.Human_logitModel))
    GET_QUALITY_METRICS(params.ribotie_training_outputs, PhyloCSFpp_db)
    RIBOTIE_POSTANALYSIS(params.ribotie_training_outputs, orfanage_gtf, final_expression, final_classification)
}
