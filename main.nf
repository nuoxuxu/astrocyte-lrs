include { PREPROCESSING } from "./subworkflows/local/preprocessing"
include { ISOSEQ } from "./subworkflows/local/isoseq"
include { RUN_OARFISH } from "./subworkflows/local/oarfish"
include { SQANTI } from "./subworkflows/local/sqanti"
include { FILTER_BY_EXPRESSION } from "./subworkflows/local/filter_by_expression"
include { RUN_ORFANAGE } from "./subworkflows/local/orfanage"
include { PREPARE_RIBOTIE } from "./subworkflows/local/prepare_ribotie"
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
    FILTER_BY_EXPRESSION(RUN_OARFISH.out.oarfish_quant, SQANTI.out.filtered_classification, SQANTI.out.filtered_gtf, SQANTI.out.sqanti_corrected_fasta, annotation_gtf)
    RUN_ORFANAGE(ref_genome_fasta, FILTER_BY_EXPRESSION.out.final_transcripts_gtf, FILTER_BY_EXPRESSION.out.final_classification, annotation_gtf)
    gtf_for_ribotie = channel.of(["gencode", file(params.annotation_gtf)])
        .concat(RUN_ORFANAGE.out.orfanage_gtf)
    PREPARE_RIBOTIE(gtf_for_ribotie, SQANTI.out.star_genomeDir, params.riboseq_unmapped_to_contaminants, ref_genome_fasta, chrom_sizes, chromAlias, channel.fromPath(params.stimulation_labels))

    PREPARE_RIBOTIE.out.ribotie_db
        .map { input_gtf_name, gtf_h5, ribotie_h5, merged_or_separate ->
            def output = [
                input_gtf_name: "${input_gtf_name}",
                gtf_h5: "nextflow_results/prepare_ribotie/${input_gtf_name}/${merged_or_separate}/${gtf_h5.name}",
                ribotie_h5: "nextflow_results/prepare_ribotie/${input_gtf_name}/${merged_or_separate}/${ribotie_h5.name}",
                merged_or_separate: "${merged_or_separate}",
                ribotie_yml: "nextflow_results/prepare_ribotie/${input_gtf_name}/${merged_or_separate}/RiboTIE.yml"
            ]
            [input_gtf_name, groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([output]))]
        }
        .collectFile(storeDir: 'nextflow_results/manifests') { input_gtf_name, text ->
            ["ribotie_training_inputs_${input_gtf_name}.json", text]
        }

    PREPARE_RIBOTIE.out.ribotie_db
        .map { input_gtf_name, _gtf_h5, _ribotie_h5, merged_or_separate ->
            def base_dir = "nextflow_results/ribotie/${input_gtf_name}/${merged_or_separate}"
            def entry = [:]
            if (merged_or_separate == "merged") {
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
            [input_gtf_name, groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([entry]))]
        }
        .collectFile(storeDir: 'nextflow_results/manifests') { input_gtf_name, text ->
            ["ribotie_training_outputs_${input_gtf_name}.json", text]
        }

    RUN_ORFANAGE.out.orfanage_gtf
        .join(RUN_ORFANAGE.out.orfanage_proteins)
        .map { input_gtf_name, _orf_gtf, _orf_prot ->
            def output = [[
                final_expression: "nextflow_results/sqanti3/isoseq/sqanti3_filter/final_expression.parquet",
                final_fasta: "nextflow_results/sqanti3/isoseq/sqanti3_filter/final_transcripts.fasta",
                final_classification: "nextflow_results/sqanti3/isoseq/sqanti3_filter/final_classification.parquet",
                orfanage_gtf: "nextflow_results/orfanage/${input_gtf_name}/orfanage_numbered_exons.gtf",
                orfanage_proteins: "nextflow_results/orfanage/${input_gtf_name}/orfanage_proteins.fasta",
            ]]
            [input_gtf_name, groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(output))]
        }
        .collectFile(storeDir: 'nextflow_results/manifests') { input_gtf_name, text ->
            ["main_pipeline_outputs_${input_gtf_name}.json", text]
        }
}
