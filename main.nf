include { PREPROCESSING } from "./subworkflows/preprocessing"
include { ISOSEQ } from "./subworkflows/isoseq"
include { RUN_OARFISH } from "./subworkflows/oarfish"
include { SQANTI_AND_FILTER_BY_EXP } from "./subworkflows/sqanti"
include { RUN_ORFANAGE } from "./subworkflows/orfanage"
include { PREPARE_RIBOTIE } from "./subworkflows/riboseq"
include { GET_QUALITY_METRICS } from "./subworkflows/quality"
include { ISOFORMSWITCH } from "./subworkflows/IsoformSwitchAnalyzeR/main.nf"

workflow {
    channel.value(file(params.kinnex_adapters)).set { kinnex_adapters }
    channel.value(file(params.isoseq_primers)).set{ isoseq_primers }
    channel.value(file(params.biosamples_csv)).set{ biosamples_csv }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file(params.annotation_gtf)).set{ annotation_gtf }
    channel.value(file(params.refTSS)).set { refTSS }
    channel.value(file(params.polyA_motif_list)).set { polyA_motif_list }
    channel.value(params.star_genomeDir_name).set { star_genomeDir_name }
    channel.value(file(params.pfamdb)).set { pfamdb }
    channel.value(file(params.PhyloCSFpp_db)).set { PhyloCSFpp_db }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }

    PREPROCESSING(params.hifi_reads_bam, kinnex_adapters, isoseq_primers, biosamples_csv)
    ISOSEQ(PREPROCESSING.out.flnc_bam, ref_genome_fasta)
    RUN_OARFISH(ISOSEQ.out.merged_sorted_collapsed_gtf, ref_genome_fasta, PREPROCESSING.out.flnc_bam)
    SQANTI_AND_FILTER_BY_EXP(params.short_read_fastqs, annotation_gtf, ref_genome_fasta, refTSS, polyA_motif_list, params.filter_configs, RUN_OARFISH.out.oarfish_quant, ISOSEQ.out.merged_sorted_collapsed_gtf, star_genomeDir_name)
    RUN_ORFANAGE(ref_genome_fasta, SQANTI_AND_FILTER_BY_EXP.out.final_transcripts_gtf, annotation_gtf)
    PREPARE_RIBOTIE(RUN_ORFANAGE.out.orfanage_gtf, SQANTI_AND_FILTER_BY_EXP.out.final_classification, annotation_gtf, SQANTI_AND_FILTER_BY_EXP.out.star_genomeDir, params.riboseq_unmapped_to_contaminants, ref_genome_fasta)
    GET_QUALITY_METRICS(params.ribotie_training_outputs, PhyloCSFpp_db, RUN_ORFANAGE.out.orfanage_proteins, pfamdb)
    ISOFORMSWITCH(SQANTI_AND_FILTER_BY_EXP.out.final_expression, primer_to_sample, SQANTI_AND_FILTER_BY_EXP.out.corrected_fasta, RUN_ORFANAGE.out.orfanage_gtf, annotation_gtf, SQANTI_AND_FILTER_BY_EXP.out.final_classification)
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