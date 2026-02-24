include { PREPROCESSING } from "./subworkflows/preprocessing"
include { ISOSEQ } from "./subworkflows/isoseq"
include { RUN_OARFISH } from "./subworkflows/oarfish"
include { SQANTI } from "./subworkflows/sqanti"
include { FILTER_BY_EXPRESSION } from "./subworkflows/filter_by_expression"
include { RUN_ORFANAGE } from "./subworkflows/orfanage"
include { PREPARE_RIBOTIE } from "./subworkflows/riboseq"

workflow {
    channel.value(file(params.kinnex_adapters)).set { kinnex_adapters }
    channel.value(file(params.isoseq_primers)).set{ isoseq_primers }
    channel.value(file(params.biosamples_csv)).set{ biosamples_csv }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file(params.annotation_gtf)).set{ annotation_gtf }
    channel.value(file(params.refTSS)).set { refTSS }
    channel.value(file(params.polyA_motif_list)).set { polyA_motif_list }
    channel.value(params.star_genomeDir_name).set { star_genomeDir_name }

    PREPROCESSING(params.hifi_reads_bam, kinnex_adapters, isoseq_primers, biosamples_csv)
    ISOSEQ(PREPROCESSING.out.flnc_bam, ref_genome_fasta)
    RUN_OARFISH(ISOSEQ.out.merged_sorted_collapsed_gtf, ref_genome_fasta, PREPROCESSING.out.flnc_bam)
    SQANTI(params.short_read_fastqs, annotation_gtf, ref_genome_fasta, refTSS, polyA_motif_list, ISOSEQ.out.merged_sorted_collapsed_gtf, star_genomeDir_name)
    FILTER_BY_EXPRESSION(RUN_OARFISH.out.oarfish_quant, SQANTI.out.filtered_classification, SQANTI.out.filtered_gtf, SQANTI.out.sqanti_corrected_fasta, params.filter_configs, annotation_gtf)
    RUN_ORFANAGE(ref_genome_fasta, FILTER_BY_EXPRESSION.out.final_transcripts_gtf, annotation_gtf)
    PREPARE_RIBOTIE(RUN_ORFANAGE.out.orfanage_gtf, FILTER_BY_EXPRESSION.out.final_classification, annotation_gtf, SQANTI.out.star_genomeDir, params.riboseq_unmapped_to_contaminants, ref_genome_fasta, FILTER_BY_EXPRESSION.out.tmap)

    PREPARE_RIBOTIE.out.ribotie_db
        .toList()
        .map { results ->
            def outputs = results.collect { param_name, gtf_h5, ribotie_h5, mode ->
                [
                    param_set_name: param_name,
                    gtf_h5: "nextflow_results/ribotie/${param_name}/${mode}/${gtf_h5.name}",
                    ribotie_h5: "nextflow_results/ribotie/${param_name}/${mode}/${ribotie_h5.name}",
                    base_dir: "nextflow_results/ribotie/${param_name}/${mode}",
                    ribotie_yml: "nextflow_results/ribotie/${param_name}/${mode}/RiboTIE.yml"
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
                    orfanage_gtf: "nextflow_results/orfanage/${param_name}/orfanage.gtf",
                    orfanage_proteins: "nextflow_results/orfanage/${param_name}/orfanage_proteins.fasta",
                ]
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'main_pipeline_outputs.json', storeDir: 'nextflow_results/manifests')
}