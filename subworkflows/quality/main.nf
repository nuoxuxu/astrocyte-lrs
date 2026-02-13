process pfam_scan {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/${param_set_name}"

    input:
    tuple val(param_set_name), path(translation_fasta)
    path(pfamdb)
    
    output:
    tuple val(param_set_name), path("pfam_scan_results.csv")
    
    script:
    """
    pfam_scan.py \\
        -out pfam_scan_results.csv \\
        -outfmt csv \\
        $translation_fasta \\
        $pfamdb
    """
}

process phylocsfpp {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/${param_set_name}"

    input:
    tuple val(param_set_name), val(condition), path(ribotie_output_gtf), path(phyloCSF_db)

    output:
    tuple val(param_set_name), val(condition), path("${ribotie_output_gtf.baseName}.PhyloCSF++.gtf")

    script:
    """
    phylocsf++ annotate-with-tracks $phyloCSF_db/PhyloCSF+1.bw $ribotie_output_gtf
    """
}

process run_cpat {
    module "python:gcc:arrow/19.0.1:rust:r/4.4.0"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}"
    
    input:
    path(Human_coding_transcripts_CDS)
    path(Human_noncoding_transcripts_RNA)
    path(Human_logitModel)
    tuple val(param_set_name), path(nt_fasta)

    output:
    path("CPAT.ORF_seqs.fa")

    script:
    """
    make_hexamer_tab -c $Human_coding_transcripts_CDS -n $Human_noncoding_transcripts_RNA > Human_Hexamer.tsv

    cpat \\
        -x Human_Hexamer.tsv \\
        -d $Human_logitModel \\
        -g $nt_fasta \\
        --min-orf=50 \\
        --top-orf=50 \\
        -o CPAT \\
        1> CPAT.output \\
        2> CPAT.error
    """

}

workflow GET_QUALITY_METRICS {
    take:
    ribotie_training_outputs
    PhyloCSFpp_db
    translation_fasta
    pfamdb
    nt_fasta
    Human_coding_transcripts_CDS
    Human_noncoding_transcripts_RNA
    Human_logitModel

    main:
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple(entry.param_set_name, "unstim", file(entry.ribotie_unstim_novel_gtf))
        }
        .set { unstim_input_ch }

    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple(entry.param_set_name, "stim", file(entry.ribotie_stim_novel_gtf))
        }
        .set { stim_input_ch }        

    unstim_input_ch
        .mix(stim_input_ch)
        .combine(PhyloCSFpp_db)
        | phylocsfpp

    run_cpat(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, nt_fasta)

    // pfam_scan(translation_fasta, channel.fromPath(pfamdb))
}