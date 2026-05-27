include { GET_QUALITY_METRICS } from "./subworkflows/local/quality"
include { ISOFORMSWITCH as ISOFORMSWITCH_ORF } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { ISOFORMSWITCH as ISOFORMSWITCH_all } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/local/ribotie_postanalysis/main.nf"
include { FILTER_RIBOTIE } from "./subworkflows/local/filter_ribotie/main.nf"
include { CDS_LENGTH_DISTRIBUTION } from "./subworkflows/local/cds_length_distribution/main.nf"

process supplement_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/supplement_collaborator_gtf"

    input:
    path(collaborator_gtf)
    path(final_gtf)

    output:
    path("supplemented_collaborator.gtf"), emit: supplemented_gtf

    script:
    """
    supplement_collaborator_gtf.py $collaborator_gtf $final_gtf -o supplemented_collaborator.gtf
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
    channel.value(file(params.pfamdb)).set { pfamdb }
    // channel.value(file(params.PhyloCSFpp_db)).set { PhyloCSFpp_db }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_transcripts.fasta")).set { final_fasta }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_expression.parquet")).set { final_expression }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_classification.parquet")).set { final_classification }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file("from_collaborator/filtered_output.gtf")).set { collaborator_gtf }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_transcripts.gtf")).set { final_gtf }

    supplement_collaborator_gtf(collaborator_gtf, final_gtf)

    supplement_collaborator_gtf.out.supplemented_gtf
        .combine(final_fasta)
        .combine(final_expression)
    | prepare_supplemented_files

    translate_supplemented_ORFs(ref_genome_fasta, supplement_collaborator_gtf.out.supplemented_gtf)

    ISOFORMSWITCH_all(channel.value("ALL"), prepare_supplemented_files.out.supplemented_expression, primer_to_sample, prepare_supplemented_files.out.supplemented_fasta, supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, final_classification, translate_supplemented_ORFs.out.supplemented_proteins, pfamdb, file(params.Human_coding_transcripts_CDS), file(params.Human_noncoding_transcripts_RNA), file(params.Human_logitModel))
    // channel.fromPath(params.ribotie_training_outputs)
    //     .map { json_file ->
    //         def orfanage_mode = (json_file.baseName =~ /ribotie_training_outputs_(.+)/)[0][1]
    //         def entries = new groovy.json.JsonSlurper().parse(json_file)
    //         def merged_entry = entries.find { it.containsKey('ribotie_merged_gtf') }
    //         tuple(orfanage_mode, file(merged_entry.ribotie_merged_gtf))
    //     }
    //     .set { ribotie_input_ch }

    // ribotie_input_ch.filter { it[0] == "minlen" }.set{ ribotie_input_ch }

    // CDS_LENGTH_DISTRIBUTION(ribotie_input_ch, ref_genome_fasta)

    FILTER_RIBOTIE(collaborator_gtf, final_fasta, final_expression, final_classification, ref_genome_fasta)

    ISOFORMSWITCH_ORF(channel.value("ORF"), FILTER_RIBOTIE.out.filtered_RiboTIE_expression, primer_to_sample, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, collaborator_gtf, annotation_gtf, FILTER_RIBOTIE.out.filtered_RiboTIE_classification, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, file(params.Human_coding_transcripts_CDS), file(params.Human_noncoding_transcripts_RNA), file(params.Human_logitModel))
    
    // ISOFORMSWITCH_all(channel.value("ALL"), FILTER_RIBOTIE.out.filtered_RiboTIE_expression, primer_to_sample, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, collaborator_gtf, annotation_gtf, FILTER_RIBOTIE.out.filtered_RiboTIE_classification, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, file(params.Human_coding_transcripts_CDS), file(params.Human_noncoding_transcripts_RNA), file(params.Human_logitModel))
    // RIBOTIE_POSTANALYSIS(params.ribotie_training_outputs, collaborator_gtf, FILTER_RIBOTIE.out.filtered_RiboTIE_expression)
}
