process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/mid_stringency"

    input:
    tuple path(final_expression), path(orfanage_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    path("IsoformSwitchAnalyzeR.rds")

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --final_fasta $final_fasta \\
        --orfanage_gtf $orfanage_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification
    """
}

process run_cpat {
    module "python:gcc:arrow/19.0.1:rust:r/4.4.0"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/mid_stringency"

    input:
    path(Human_coding_transcripts_CDS)
    path(Human_noncoding_transcripts_RNA)
    path(Human_logitModel)
    path(final_fasta)

    output:
    path("CPAT.ORF_prob.tsv"), emit: ORF_prob_tsv
    path("CPAT.ORF_prob.best.tsv"), emit: ORF_prob_best_tsv
    path("CPAT.ORF_seqs.fa"), emit: ORF_seqs_fa
    path("CPAT.no_ORF.txt"), emit: no_ORF_txt

    script:
    """
    make_hexamer_tab -c $Human_coding_transcripts_CDS -n $Human_noncoding_transcripts_RNA > Human_Hexamer.tsv

    cpat \\
        -x Human_Hexamer.tsv \\
        -d $Human_logitModel \\
        -g $final_fasta \\
        --min-orf=50 \\
        --top-orf=50 \\
        -o CPAT \\
        1> CPAT.output \\
        2> CPAT.error
    """
}

process pfam_scan {
    conda "/scratch/nxu/astrocytes/env"
    label "long_slurm_job"
    storeDir "nextflow_results/quality/mid_stringency"

    input:
    path(translation_fasta)
    path(pfamdb)

    output:
    path("pfam_scan_results.csv")

    script:
    """
    pfam_scan.py \\
        -out pfam_scan_results.csv \\
        -outfmt csv \\
        $translation_fasta \\
        $pfamdb
    """
}

process convert_pfam_scan_results {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/mid_stringency"

    input:
    path(pfam_scan_results)

    output:
    path("pfam_scan_results_modified.txt")

    script:
    """
    convert_pfam_csv_to_txt.py $pfam_scan_results > pfam_scan_results_modified.txt
    """
}

workflow ISOFORMSWITCH {
    take:
    final_expression
    primer_to_sample
    final_fasta
    orfanage_gtf
    annotation_gtf
    final_classification
    translation_fasta
    pfamdb
    Human_coding_transcripts_CDS
    Human_noncoding_transcripts_RNA
    Human_logitModel

    main:

    final_expression
        .combine(orfanage_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(final_fasta)
        .combine(annotation_gtf)
    | IsoseqsSwitchList

    run_cpat(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, final_fasta)
    pfam_scan(translation_fasta, pfamdb)
}
