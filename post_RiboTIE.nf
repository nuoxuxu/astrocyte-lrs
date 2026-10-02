include { AIM_2 } from "./subworkflows/local/aim_2/main.nf"
include { ISOFORMSWITCH } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/local/ribotie_postanalysis/main.nf"
include { FILTER_RIBOTIE } from "./subworkflows/local/filter_ribotie/main.nf"
include { CDS_LENGTH_DISTRIBUTION } from "./subworkflows/local/cds_length_distribution/main.nf"
include { RUN_VEP } from "./subworkflows/local/vep/main.nf"
include { SUMMARY_TABLE } from "./subworkflows/local/summary_table/main.nf"
include { FIX_COLLABORATOR_ORF_TYPE } from "./subworkflows/local/fix_collaborator_ORF_type/main.nf"
include { SUPPLEMENT_COLLABORATOR_ORF } from "./subworkflows/local/supplement_collaborator_ORF/main.nf"
include { FRAGPIPE } from "./subworkflows/local/fragpipe/main.nf"
include { LRP2_PROTEIN } from "./subworkflows/local/lrp2_protein/main.nf"

process fix_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/supplemented_collaborator"

    input:
    path(collaborator_gtf)

    output:
    path("filtered_output_fixed.gtf"), emit: fixed_gtf

    script:
    """
    fix_collaborator_gtf.py $collaborator_gtf -o filtered_output_fixed.gtf
    """
}

process filter_ribotie_for_isoformswitch {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/${name}"

    input:
    tuple val(name), path(input_gtf), path(input_fasta), path(input_expression), path(input_classification)

    output:
    tuple val(name), path("${name}_filtered.fasta"),                            emit: fasta
    tuple val(name), path("${name}_filtered_expression.parquet"),               emit: expression
    tuple val(name), path("${name}_filtered_classification.parquet"),           emit: classification

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --input_fasta $input_fasta \\
        --input_expression $input_expression \\
        --input_classification $input_classification \\
        --output_fasta ${name}_filtered.fasta \\
        --output_expression ${name}_filtered_expression.parquet \\
        --output_classification ${name}_filtered_classification.parquet
    """
}

process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}"

    input:
    tuple val(version), path(final_expression), path(predicted_cds_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    tuple val(version), path("IsoformSwitchAnalyzeR.rds"), emit: rds
    path("isoformFeatures.csv"), emit: isoform_features_csv

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --final_fasta $final_fasta \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification \\
    """
}

process peptideTrackUCSC {
    storeDir "nextflow_results/proteomics"
    conda "/scratch/nxu/astrocytes/env"
    
    input:
    val mode
    path annotation_gtf
    path final_sample_classification
    path predicted_cds_gtf
    path protein_search_database
    path peptides

    output:
    path "peptides_${mode}.gtf"

    script:
    """
    make_peptide_gtf_file.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $final_sample_classification \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --protein_search_database $protein_search_database \\
        --peptides $peptides \\
        --output "peptides_${mode}.gtf"
    """
}

process peptideMapping {
    storeDir "nextflow_results/proteomics"
    conda "/scratch/nxu/astrocytes/env"

    input:
    path annotation_gtf
    path classification
    path protein_search_database
    path percolator_res
    
    output:
    path "peptide_mapping.parquet", emit: peptide_mapping
    path "novel_peptides.csv", emit: novel_peptides
    
    script:
    """
    peptide_mapping.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $classification \\
        --protein_search_database $protein_search_database \\
        --percolator_res $percolator_res
    """
}

process protein_differential_expression {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/proteomics/differential_expression"

    input:
    path(pg_matrix)

    output:
    path("protein_DE_stim_vs_unstim.tsv"),  emit: de_results
    path("log2_intensity_filtered.tsv"),    emit: log2_filtered
    path("log2_intensity_imputed.tsv"),     emit: log2_imputed
    path("sample_metadata.tsv"),            emit: sample_metadata
    path("protein_DE_stim_vs_unstim.pdf"),  emit: plots

    script:
    """
    proteomics_differential_expression.R \\
        --input $pg_matrix \\
        --outdir . \\
        --figdir .
    """
}

process SUMMARY {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary"

    input:
    path(ribotie_csv)
    path(isoform_features_csv)
    path(protein_de)

    output:
    path("ribotie_summary.csv"), emit: summary

    script:
    """
    make_ribotie_summary.py $ribotie_csv $isoform_features_csv $protein_de -o ribotie_summary.csv
    """
}

workflow {
    // Data flow: orfanage ORFs → RiboTIE scoring (by collaborator) → filtered_output.gtf (high-confidence RiboTIE hits)
    // supplement_collaborator_gtf re-adds orfanage ORFs filtered out by RiboTIE, creating a comprehensive ORF set

    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_transcripts.fasta")).set { final_fasta }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_expression.parquet")).set { final_expression }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_classification.parquet")).set { final_classification }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file("from_collaborator/filtered_output.gtf")).set { collaborator_gtf }
    channel.value(file("from_collaborator/ribotie_cpm1_3sample.csv")).set { ribotie_cpm1_3sample }
    channel.value(file("from_collaborator/concordant_all_three_exp.csv")).set { concordant_csv }
    channel.value(file(params.bigbrain_sqtl)).set { bigbrain_sqtl }
    channel.value(file(params.bigbrain_coloc)).set { bigbrain_coloc }
    channel.value(file(params.leafcutter_sig)).set { leafcutter_sig }
    channel.value(file(params.leafcutter_clu2gene)).set { leafcutter_clu2gene }
    channel.value(file("nextflow_results/orfanage/minlen/orfanage.gtf")).set { orfanage_gtf }
    channel.value(file(params.Human_coding_transcripts_CDS)).set { Human_coding_transcripts_CDS }
    channel.value(file(params.Human_noncoding_transcripts_RNA)).set { Human_noncoding_transcripts_RNA }
    channel.value(file(params.Human_logitModel)).set { Human_logitModel }
    channel.value(file(params.pfamdb)).set { pfamdb }
    channel.value(file("data/study2_orfs.gtf")).set { study2_gtf }

    fix_collaborator_gtf(collaborator_gtf)

    SUPPLEMENT_COLLABORATOR_ORF(
        fix_collaborator_gtf.out.fixed_gtf,
        orfanage_gtf,
        final_fasta,
        final_expression,
        ref_genome_fasta
    )

    FIX_COLLABORATOR_ORF_TYPE(
        fix_collaborator_gtf.out.fixed_gtf,
        final_classification,
        annotation_gtf,
        ribotie_cpm1_3sample
    )

    // Build per-version channel from manifests (exclude gencode)
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
        .map { name, ribotie_gtf, expr, fasta, classif ->
            tuple(name, ribotie_gtf, fasta, expr, classif)
        }
        .set { ribotie_inputs_ch }

    // Filter RiboTIE GTF to match expression data (versioned)
    filter_ribotie_for_isoformswitch(ribotie_inputs_ch)

    // Create a channel with the isoform information
    ribotie_inputs_ch.map { it.take(2) }
        .join(filter_ribotie_for_isoformswitch.out.expression)
        .map { version, ribotie_gtf, filtered_expression ->
            tuple(version, filtered_expression, ribotie_gtf)
        }
        .join(filter_ribotie_for_isoformswitch.out.classification)
        .combine(primer_to_sample)
        .join(filter_ribotie_for_isoformswitch.out.fasta)
        .combine(annotation_gtf)
        .set { my_isoform_ch }

    channel.value("supplemented_collaborator")
        .combine(SUPPLEMENT_COLLABORATOR_ORF.out.supplemented_expression)
        .combine(SUPPLEMENT_COLLABORATOR_ORF.out.supplemented_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(SUPPLEMENT_COLLABORATOR_ORF.out.supplemented_fasta)
        .combine(annotation_gtf)
        .set { collaborator_isoform_ch }

    my_isoform_ch.mix(collaborator_isoform_ch) | ISOFORMSWITCH

    FRAGPIPE(
        file(params.fragpipe_workflow),
        file(params.fragpipe_manifest),
        file(params.fragpipe_database),
        file(params.fragpipe_sif)
    )

    // Stim (+) vs Unstim (-) protein differential expression (limma) on the DIA-NN protein-group matrix
    protein_differential_expression(FRAGPIPE.out.pg_matrix)

    peptideTrackUCSC(
        "collaborator", 
        annotation_gtf, 
        final_classification, 
        fix_collaborator_gtf.out.fixed_gtf, 
        file("nextflow_results/ribotie/filtered/filtered_RiboTIE_proteins.fasta"), 
        FRAGPIPE.out.peptides
    )
    peptideMapping(
        annotation_gtf, 
        final_classification, 
        file("nextflow_results/ribotie/filtered/filtered_RiboTIE_proteins.fasta"), 
        FRAGPIPE.out.peptides
    )

    // Per-ORF table with differential splicing (supplemented_collaborator isoformFeatures is keyed by ORF_id)
    // and Stim vs Unstim protein differential expression columns
    SUMMARY(
        FIX_COLLABORATOR_ORF_TYPE.out.ribotie_res_merged_fixed_with_lncRNA,
        ISOFORMSWITCH.out.versioned_isoform_features_csv
            .filter { version, _csv -> version == "supplemented_collaborator" }
            .map { _version, csv -> csv },
        protein_differential_expression.out.de_results
    )

    LRP2_PROTEIN(
        "collaborator",
        fix_collaborator_gtf.out.fixed_gtf,
        final_fasta,
        ribotie_cpm1_3sample,
        annotation_gtf
    )

    // AIM_2(supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, bigbrain_sqtl, bigbrain_coloc, ISOFORMSWITCH.out.isoform_features_csv, leafcutter_sig, leafcutter_clu2gene)

    // FILTER_RIBOTIE(fix_collaborator_gtf.out.fixed_gtf, final_fasta, final_expression, final_classification, ref_genome_fasta)

    // SUMMARY_TABLE(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, fix_collaborator_gtf.out.fixed_gtf, ribotie_cpm1_3sample, annotation_gtf, orfanage_gtf, final_classification, ISOFORMSWITCH.out.isoform_features_csv, study2_gtf, AIM_2.out.leafcutter_coloc, AIM_2.out.novel_coding_junction_coloc, concordant_csv)

    //TODO: MAPS analysis for variants disruption ncORFs
    // channel.value(file("/scratch/nxu/100KGP_splicing/data/gnomad/exomes/gnomad.exomes.v4.1.sites.chr16.vcf.bgz")).map { ["gnomad_exomes_chr16", it] }.set { vcf_ch }
    // channel.value(file(params.vep_uorf_data)).set { vep_uorf_data }
    // RUN_VEP(vcf_ch, vep_uorf_data, ref_genome_fasta, annotation_gtf)
}
