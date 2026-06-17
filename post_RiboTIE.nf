include { AIM_2 } from "./subworkflows/local/aim_2/main.nf"
include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { ISOFORMSWITCH } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/local/ribotie_postanalysis/main.nf"
include { FILTER_RIBOTIE } from "./subworkflows/local/filter_ribotie/main.nf"
include { CDS_LENGTH_DISTRIBUTION } from "./subworkflows/local/cds_length_distribution/main.nf"
include { RUN_VEP } from "./subworkflows/local/vep/main.nf"
include { SUMMARY_TABLE } from "./subworkflows/local/summary_table/main.nf"

process fix_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/supplement_collaborator_gtf"

    input:
    path(collaborator_gtf)

    output:
    path("filtered_output_fixed.gtf"), emit: fixed_gtf

    script:
    """
    fix_collaborator_gtf.py $collaborator_gtf -o filtered_output_fixed.gtf
    """
}

process supplement_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/supplement_collaborator_gtf"

    input:
    path(collaborator_gtf)
    path(orfanage_gtf)

    output:
    path("supplemented_collaborator.gtf"), emit: supplemented_gtf

    script:
    """
    supplement_collaborator_gtf.py $collaborator_gtf $orfanage_gtf -o supplemented_collaborator.gtf
    """
}

process prepare_supplemented_files {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/supplement_collaborator_gtf"

    input:
    tuple path(supplemented_gtf), path(final_fasta), path(final_expression)

    output:
    path("supplemented_collaborator.fasta"),                    emit: supplemented_fasta
    path("supplemented_collaborator_expression.parquet"),       emit: supplemented_expression

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    prepare_supplemented_gtf_files.py \\
        $supplemented_gtf \\
        $final_fasta \\
        $final_expression \\
        --output_fasta supplemented_collaborator.fasta \\
        --output_expression supplemented_collaborator_expression.parquet
    """
}

process translate_supplemented_ORFs {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/supplement_collaborator_gtf"

    input:
    path ref_genome_fasta
    path supplemented_gtf

    output:
    path("supplemented_collaborator_proteins.fasta"), emit: supplemented_proteins

    script:
    """
    gffread -y supplemented_collaborator_proteins.fasta \\
        -g $ref_genome_fasta \\
        $supplemented_gtf
    """
}

workflow {
    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_transcripts.fasta")).set { final_fasta }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_expression.parquet")).set { final_expression }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_classification.parquet")).set { final_classification }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file("from_collaborator/filtered_output.gtf")).set { collaborator_gtf }
    channel.value(file("from_collaborator/ribotie_cpm1_3sample.csv")).set { ribotie_cpm1_3sample }
    channel.value(file("from_collaborator/concordant_all_three_exp.csv")).set { concordant_csv }
    channel.value(file(params.bigbrain_sqtl)).set { bigbrain_sqtl }
    channel.value(file(params.bigbrain_coloc)).set { bigbrain_coloc }
    channel.value(file(params.leafcutter_sig)).set { leafcutter_sig }
    channel.value(file(params.leafcutter_clu2gene)).set { leafcutter_clu2gene }
    channel.value(file("nextflow_results/orfanage/minlen/mid_stringency/orfanage.gtf")).set { orfanage_gtf }
    channel.value(file(params.Human_coding_transcripts_CDS)).set { Human_coding_transcripts_CDS }
    channel.value(file(params.Human_noncoding_transcripts_RNA)).set { Human_noncoding_transcripts_RNA }
    channel.value(file(params.Human_logitModel)).set { Human_logitModel }
    channel.value(file(params.pfamdb)).set { pfamdb }
    channel.value(file("data/study2_orfs.gtf")).set { study2_gtf }

    fix_collaborator_gtf(collaborator_gtf)

    supplement_collaborator_gtf(fix_collaborator_gtf.out.fixed_gtf, orfanage_gtf)

    supplement_collaborator_gtf.out.supplemented_gtf
        .combine(final_fasta)
        .combine(final_expression)
    | prepare_supplemented_files

    translate_supplemented_ORFs(ref_genome_fasta, supplement_collaborator_gtf.out.supplemented_gtf)

    ISOFORMSWITCH(channel.value("ALL2"), prepare_supplemented_files.out.supplemented_expression, primer_to_sample, prepare_supplemented_files.out.supplemented_fasta, supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, final_classification)

    AIM_2(supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, bigbrain_sqtl, bigbrain_coloc, ISOFORMSWITCH.out.isoform_features_csv, leafcutter_sig, leafcutter_clu2gene)

    FILTER_RIBOTIE(fix_collaborator_gtf.out.fixed_gtf, final_fasta, final_expression, final_classification, ref_genome_fasta)

    SUMMARY_TABLE(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, fix_collaborator_gtf.out.fixed_gtf, ribotie_cpm1_3sample, annotation_gtf, orfanage_gtf, final_classification, ISOFORMSWITCH.out.isoform_features_csv, study2_gtf, AIM_2.out.leafcutter_coloc, AIM_2.out.novel_coding_junction_coloc, concordant_csv)

    //TODO: MAPS analysis for variants disruption ncORFs
    // channel.value(file("/scratch/nxu/100KGP_splicing/data/gnomad/exomes/gnomad.exomes.v4.1.sites.chr16.vcf.bgz")).map { ["gnomad_exomes_chr16", it] }.set { vcf_ch }
    // channel.value(file(params.vep_uorf_data)).set { vep_uorf_data }
    // RUN_VEP(vcf_ch, vep_uorf_data, ref_genome_fasta, annotation_gtf)
}
