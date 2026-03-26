include { PREPROCESSING } from "./subworkflows/local/preprocessing"
include { ISOSEQ } from "./subworkflows/local/isoseq"
include { RUN_OARFISH } from "./subworkflows/local/oarfish"
include { SQANTI } from "./subworkflows/local/sqanti"
include { FILTER_BY_EXPRESSION } from "./subworkflows/local/filter_by_expression"
include { RUN_ORFANAGE } from "./subworkflows/local/orfanage"
include { PREPARE_RIBOTIE } from "./subworkflows/local/prepare_ribotie"
include { RUN_RIBOTISH } from "./subworkflows/local/ribotish"
include { QUANTIFY_PSEUDO_ALIGNMENT } from "./subworkflows/nf-core/quantify_pseudo_alignment"
include { SALMON_INDEX } from "./modules/nf-core/salmon/index"
include { ANOTA2SEQ_ANOTA2SEQRUN } from "./modules/nf-core/anota2seq/anota2seqrun"

workflow {
    channel.value(file(params.kinnex_adapters)).set { kinnex_adapters }
    channel.value(file(params.isoseq_primers)).set{ isoseq_primers }
    channel.value(file(params.biosamples_csv)).set{ biosamples_csv }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file(params.annotation_gtf)).set{ annotation_gtf }
    channel.value(file(params.refTSS)).set { refTSS }
    channel.value(file(params.polyA_motif_list)).set { polyA_motif_list }
    channel.value(params.star_genomeDir_name).set { star_genomeDir_name }
    channel.value(params.chrom_sizes).set { chrom_sizes }
    channel.value(params.chromAlias).set { chromAlias }

    PREPROCESSING(params.hifi_reads_bam, kinnex_adapters, isoseq_primers, biosamples_csv)
    ISOSEQ(PREPROCESSING.out.flnc_bam, ref_genome_fasta)
    RUN_OARFISH(ISOSEQ.out.merged_sorted_collapsed_gtf, ref_genome_fasta, PREPROCESSING.out.flnc_bam)
    SQANTI(params.short_read_fastqs, annotation_gtf, ref_genome_fasta, refTSS, polyA_motif_list, ISOSEQ.out.merged_sorted_collapsed_gtf, star_genomeDir_name)
    FILTER_BY_EXPRESSION(RUN_OARFISH.out.oarfish_quant, SQANTI.out.filtered_classification, SQANTI.out.filtered_gtf, SQANTI.out.sqanti_corrected_fasta, params.filter_configs, annotation_gtf)
    RUN_ORFANAGE(ref_genome_fasta, FILTER_BY_EXPRESSION.out.final_transcripts_gtf, annotation_gtf)
    PREPARE_RIBOTIE(FILTER_BY_EXPRESSION.out.final_transcripts_gtf, FILTER_BY_EXPRESSION.out.final_classification, annotation_gtf, SQANTI.out.star_genomeDir, params.riboseq_unmapped_to_contaminants, ref_genome_fasta, chrom_sizes, chromAlias)

    PREPARE_RIBOTIE.out.ribotie_db
        .toList()
        .map { results ->
            def outputs = results.collect { param_name, gtf_h5, ribotie_h5, mode ->
                [
                    param_set_name: param_name,
                    gtf_h5: "nextflow_results/prepare_ribotie/${param_name}/${mode}/${gtf_h5.name}",
                    ribotie_h5: "nextflow_results/prepare_ribotie/${param_name}/${mode}/${ribotie_h5.name}",
                    mode: "${mode}",
                    ribotie_yml: "nextflow_results/prepare_ribotie/${param_name}/${mode}/RiboTIE.yml"
                ]
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'ribotie_training_inputs.json', storeDir: 'nextflow_results/manifests')

    PREPARE_RIBOTIE.out.ribotie_db
        .toList()
        .map { results ->
            def outputs = results.collect { param_name, _gtf_h5, _ribotie_h5, mode ->
                def base_dir = "nextflow_results/ribotie/${param_name}/${mode}"
                def entry = [
                    param_set_name: param_name,
                ]
                if (mode == "merged") {
                    entry += [
                        ribotie_merged_gtf: "${base_dir}/ribotie_res_merged.gtf",
                        ribotie_merged_csv: "${base_dir}/ribotie_res_merged.csv",
                        ribotie_merged_novel_gtf: "${base_dir}/ribotie_res_merged.novel.gtf",
                        ribotie_merged_novel_csv: "${base_dir}/ribotie_res_merged.novel.csv",
                        ribotie_merged_redundant_gtf: "${base_dir}/ribotie_res_merged.redundant.gtf",
                        ribotie_merged_redundant_csv: "${base_dir}/ribotie_res_merged.redundant.csv",
                    ]
                } else {
                    entry += [
                        ribotie_unstim_gtf: "${base_dir}/ribotie_res_Unstim.gtf",
                        ribotie_unstim_csv: "${base_dir}/ribotie_res_Unstim.csv",
                        ribotie_stim_gtf: "${base_dir}/ribotie_res_Stim.gtf",
                        ribotie_stim_csv: "${base_dir}/ribotie_res_Stim.csv",
                        ribotie_unstim_novel_gtf: "${base_dir}/ribotie_res_Unstim.novel.gtf",
                        ribotie_unstim_novel_csv: "${base_dir}/ribotie_res_Unstim.novel.csv",
                        ribotie_stim_novel_gtf: "${base_dir}/ribotie_res_Stim.novel.gtf",
                        ribotie_stim_novel_csv: "${base_dir}/ribotie_res_Stim.novel.csv",
                    ]
                }
                entry
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'ribotie_training_outputs.json', storeDir: 'nextflow_results/manifests')

    FILTER_BY_EXPRESSION.out.final_expression
        .filter { param_set_name, _final_expression ->
            param_set_name == "mid_stringency"
        }
        .join(FILTER_BY_EXPRESSION.out.final_fasta)
        .join(FILTER_BY_EXPRESSION.out.final_classification)
        .join(RUN_ORFANAGE.out.orfanage_gtf)
        .join(RUN_ORFANAGE.out.orfanage_proteins)
        .toList()
        .map { results ->
            def outputs = results.collect { param_name, _final_expression, _final_fasta, _final_classification, _orfanage_gtf, _orfanage_proteins ->
                [
                    param_set_name: param_name,
                    final_expression: "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_name}/final_expression.parquet",
                    final_fasta: "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_name}/final_transcripts.fasta",
                    final_classification: "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_name}/final_classification.parquet",
                    orfanage_gtf: "nextflow_results/orf_prediction/orfanage/${param_name}/orfanage.gtf",
                    orfanage_proteins: "nextflow_results/orf_prediction/orfanage/${param_name}/orfanage_proteins.fasta",
                ]
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'main_pipeline_outputs.json', storeDir: 'nextflow_results/manifests')
    
    RUN_RIBOTISH(PREPARE_RIBOTIE.out.genome_bam, FILTER_BY_EXPRESSION.out.final_transcripts_gtf, RUN_ORFANAGE.out.orfanage_gtf, ref_genome_fasta)
    
    // reads = PREPARE_RIBOTIE.out.transcriptome_bam
    //     .map{
    //         _param_set_name, bam ->
    //         def meta = [:]
    //         meta.id = bam.name.minus("_Unmapped.Aligned.toTranscriptome.out.bam")
    //         [ meta, bam ]
    //     }
    // ch_samplesheet = channel.value(file("data/samplesheet.csv", checkIfExists: true))

    // transcript_fasta = FILTER_BY_EXPRESSION.out.final_fasta
    //     .filter{ param_set_name, _final_fasta -> param_set_name == "mid_stringency" }
    //     .map{ _param_set_name, final_fasta -> final_fasta }

    // ref_genome_fasta
    //     .combine(transcript_fasta)
    //     | SALMON_INDEX

    // QUANTIFY_PSEUDO_ALIGNMENT(
    //     ch_samplesheet.map { [ [:], it ] },
    //     reads,
    //     SALMON_INDEX.out,
    //     transcript_fasta,
    //     FILTER_BY_EXPRESSION.out.final_transcripts_gtf.filter{ param_set_name, _final_transcripts_gtf -> param_set_name == "mid_stringency" },
    //     "gene_id",
    //     null,
    //     "salmon",
    //     null,
    //     null,
    //     null
    // )

    // channel.value(file(params.ch_contrasts_file)).set{ ch_contrasts_file }
    // if (ch_contrasts_file){

    //     ch_contrasts = ch_contrasts_file
    //         .splitCsv ( header:true, sep:',' )
    //         .map{[it, it.variable, it.reference, it.target]}

    //     ch_samplesheet_matrix = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled
    //         .combine(ch_samplesheet)
    //         .map{[it[0], it[2], it[1]]}
    //         .first()

    //     ANOTA2SEQ_ANOTA2SEQRUN(
    //         ch_contrasts,
    //         ch_samplesheet_matrix
    //     )
    // }
}