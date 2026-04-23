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

process split_FASTA {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"

    input:
    path(translation_fasta)

    output:
    path("chunk_*.fa")

    script:
    """
    python3 - << 'EOF'
    seqs = []
    header, seq = None, []
    with open("${translation_fasta}") as f:
        for line in f:
            if line.startswith(">"):
                if header:
                    seqs.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
    if header:
        seqs.append((header, "".join(seq)))

    n = min(int("${task.cpus}") // 8, len(seqs))
    chunk_size = max(1, -(-len(seqs) // n))
    for i in range(n):
        chunk = seqs[i * chunk_size:(i + 1) * chunk_size]
        if not chunk:
            break
        with open(f"chunk_{i+1:04d}.fa", "w") as out:
            for h, s in chunk:
                out.write(h + s)
    EOF
    """
}

process pfam_scan {
    conda "/scratch/nxu/astrocytes/env"
    label "mid_slurm_job"

    input:
    path(translation_fasta)
    path(pfamdb)

    output:
    path("pfam_scan_results.csv")

    script:
    """
    pfam_scan.py \\
        -cpu 8 \\
        -out pfam_scan_results.csv \\
        -outfmt csv \\
        $translation_fasta \\
        $pfamdb
    """
}

process merge_pfam_results {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/mid_stringency"

    input:
    path(pfam_results)

    output:
    path("pfam_scan_results.csv")

    script:
    """
    head -1 \$(ls *.csv | head -1) > pfam_scan_results.csv
    for f in *.csv; do
        tail -n +2 "\$f" >> pfam_scan_results.csv
    done
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

    split_FASTA(translation_fasta)
    pfam_scan(split_FASTA.out.flatten(), pfamdb)
    merge_pfam_results(pfam_scan.out.collect())
}
