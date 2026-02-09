include { PREPROCESSING } from "./subworkflows/preprocessing"
include { ISOSEQ } from "./subworkflows/isoseq"
include { RUN_OARFISH } from "./subworkflows/oarfish"
include { SQANTI_AND_FILTER_BY_EXP } from "./subworkflows/sqanti"
include { RUN_ORFANAGE } from "./subworkflows/orfanage"
include { PREPARE_RIBOTIE } from "./subworkflows/riboseq"

workflow {    
    PREPROCESSING(params.hifi_reads_bam, params.kinnex_adapters, params.isoseq_primers, params.biosamples_csv)
    ISOSEQ(PREPROCESSING.flnc_bam, params.ref_genome_fasta)
    RUN_OARFISH(ISOSEQ.merged_sorted_collapsed_gtf)
    SQANTI_AND_FILTER_BY_EXP(params.input_rna_fastq, params.annotation_gtf, params.ref_genome_fasta, params.refTSS, params.polyA_motif_list, RUN_OARFISH.oarfish_quant, ISOSEQ.merged_sorted_collapsed_gtf)
    RUN_ORFANAGE(params.ref_genome_fasta, SQANTI_AND_FILTER_BY_EXP.final_transcripts_gtf)
    PREPARE_RIBOTIE(RUN_ORFANAGE.orfanage_gtf, SQANTI_AND_FILTER_BY_EXP.final_classification, params.annotation_gtf, SQANTI_AND_FILTER_BY_EXP.star_genomeDir, params.riboseq_unmapped_to_contaminants)
    
    // -------------------Testing PREPARE_RIBOTIE-----------------
    // orfanage_gtf = channel.of(
    //     ["low_stringency", "/scratch/nxu/astrocytes/nextflow_results/orfanage/low_stringency/orfanage.gtf"],
    //     ["high_stringency", "/scratch/nxu/astrocytes/nextflow_results/orfanage/high_stringency/orfanage.gtf"]
    // )
    // final_classification = channel.of(
    //     ["low_stringency", "/scratch/nxu/astrocytes/nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/low_stringency/final_classification.parquet"],
    //     ["high_stringency", "/scratch/nxu/astrocytes/nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/high_stringency/final_classification.parquet"]
    // )
    // star_genomeDir = channel.fromPath("/scratch/nxu/astrocytes/nextflow_results/align/star/STAR_index_v47")
    // PREPARE_RIBOTIE(orfanage_gtf, final_classification, params.annotation_gtf, star_genomeDir, params.riboseq_unmapped_to_contaminants)
    // ------------------------------------------------------------
    

    // star_riboseq_custom.out.bedGraph_UniqueMultiple
    //     .branch {
    //         unstim: it.name =~ /^merged_astro_[AB]_/
    //         stim: it.name =~ /^(merged_astro_C_|Astro_D_)/
    //     }
    //     .set { custom_bedgraphs }
    
    // star_riboseq_gencode.out.bedGraph_UniqueMultiple
    //     .branch {
    //         unstim: it.name =~ /^merged_astro_[AB]_/
    //         stim: it.name =~ /^(merged_astro_C_|Astro_D_)/
    //     }
    //     .set { gencode_bedgraphs }
    

    // merge_bg_and_convert_to_bw_custom_unstim(
    //     custom_bedgraphs.unstim.collect(),
    //     params.chrom_sizes,
    //     "custom_unstim"
    // )
    
    // Merge stimulated samples (C and D)
    // merge_bg_and_convert_to_bw_custom_stim(
    //     custom_bedgraphs.stim.collect(),
    //     params.chrom_sizes,
    //     "custom_stim"
    // )
    
    // format_gtf_for_ribotie(fixORFanageFormat.out, filter_by_expression.out.final_classification, params.annotation_gtf)
    // make_db_files(sqanti_filter.out.filtered_gtf)

    // sqanti_qc_bambu(bambu.out.supported_tx_gtf, params.annotation_gtf, params.ref_genome_fasta, params.refTSS, params.polyA_motif_list, star.out.star_aligned_bam.collect(), star.out.star_sj_tab.collect())
    // def bambu_corrected_gtf = sqanti_qc_bambu.out
    //     .map { dir -> dir / "sqanti_qc_results_corrected.gtf" }
    // def bambu_classification = sqanti_qc_bambu.out
    //     .map { dir -> dir / "sqanti_qc_results_classification.gtf" }
    // sqanti_filter_bambu(bambu_corrected_gtf, bambu_classification)
}