include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { ISOFORMSWITCH } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/local/ribotie_postanalysis/main.nf"
include { FILTER_RIBOTIE } from "./subworkflows/local/filter_ribotie/main.nf"
include { CDS_LENGTH_DISTRIBUTION } from "./subworkflows/local/cds_length_distribution/main.nf"

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }
    channel.value(file(params.pfamdb)).set { pfamdb }
    // channel.value(file(params.PhyloCSFpp_db)).set { PhyloCSFpp_db }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_transcripts.fasta")).set { final_fasta }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_expression.parquet")).set { final_expression }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_classification.parquet")).set { final_classification }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }

    channel.fromPath(params.ribotie_training_outputs)
        .map { json_file ->
            def orfanage_mode = (json_file.baseName =~ /ribotie_training_outputs_(.+)/)[0][1]
            def entries = new groovy.json.JsonSlurper().parse(json_file)
            def merged_entry = entries.find { it.containsKey('ribotie_merged_gtf') }
            tuple(orfanage_mode, file(merged_entry.ribotie_merged_gtf))
        }
        .set { ribotie_input_ch }

    CDS_LENGTH_DISTRIBUTION(ribotie_input_ch, ref_genome_fasta)

    // FILTER_RIBOTIE(ribotie_input_ch.map { it[1] }, channel.value(file(params.ribotie_cpm1_3sample)), final_fasta, final_expression, final_classification, ref_genome_fasta)

    // ISOFORMSWITCH(FILTER_RIBOTIE.out.filtered_RiboTIE_expression, primer_to_sample, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, FILTER_RIBOTIE.out.filtered_RiboTIE_gtf, annotation_gtf, FILTER_RIBOTIE.out.filtered_RiboTIE_classification, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, file(params.Human_coding_transcripts_CDS), file(params.Human_noncoding_transcripts_RNA), file(params.Human_logitModel))

    // GET_QUALITY_METRICS(params.ribotie_training_outputs, PhyloCSFpp_db)
    // RIBOTIE_POSTANALYSIS(params.ribotie_training_outputs, FILTER_RIBOTIE.out.filtered_RiboTIE_gtf, FILTER_RIBOTIE.out.filtered_RiboTIE_expression)
}
