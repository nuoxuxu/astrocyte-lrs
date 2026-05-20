process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR"

    input:
    tuple path(final_expression), path(orfanage_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf), path(cpat_results), path(pfam_results)

    output:
    path("IsoformSwitchAnalyzeR.rds"), emit: rds
    path("isoformFeatures.csv"), emit: isoform_features_csv

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --final_fasta $final_fasta \\
        --orfanage_gtf $orfanage_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification \\
        --cpat_results $cpat_results \\
        --pfam_results $pfam_results \\
    """
}

process run_cpat {
    module "python:gcc:arrow/19.0.1:rust:r/4.4.0"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/cpat"

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
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/pfam/split_fasta"

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

    n = max(1, min(int("${task.cpus}") // 8, len(seqs)))
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
    module "hmmer/3.4"
    label "mid_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/pfam"

    input:
    path(translation_fasta)
    path(pfamdb)

    output:
    path("${translation_fasta.simpleName}.csv")

    script:
    """
    pfam_scan.py \\
        -cpu 8 \\
        -out ${translation_fasta.simpleName}.csv \\
        -outfmt csv \\
        $translation_fasta \\
        $pfamdb
    """
}

process convert_cpat_format {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/cpat"

    input:
    path(orf_prob_best_tsv)

    output:
    path("cpat_results.txt")

    script:
    """
    convert_cpat_format.py $orf_prob_best_tsv -o cpat_results.txt
    """
}

process filter_cpat_results {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/cpat"

    input:
    path(cpat_results)
    path(ribotie_fasta)

    output:
    path("cpat_results_filtered.txt")

    script:
    """
    filter_cpat_results.py $cpat_results $ribotie_fasta -o cpat_results_filtered.txt
    """
}

process merge_pfam_results {
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/pfam"

    input:
    path(pfam_results)

    output:
    path("pfam_scan_results.csv")

    script:
    """
    files=($pfam_results)
    head -1 "\${files[0]}" > pfam_scan_results.csv
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> pfam_scan_results.csv
    done
    """
}

process convert_pfam_scan_results {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/pfam"

    input:
    path(pfam_scan_results)

    output:
    path("pfam_scan_results_modified.txt")

    script:
    """
    convert_pfam_csv_to_txt.py $pfam_scan_results pfam_scan_results_modified.txt
    """
}

process PlotIsoformConsequences {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/figures"

    input:
    path(rds)

    output:
    path("splicing_consequences.pdf")
    path("switch_consequences.pdf")

    script:
    """
    plot_isoform_consequences.R --rds $rds
    """
}

process volcano_plot {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/figures"

    input:
    path(rds)

    output:
    path("volcano_DGE.pdf")
    path("volcano_DTU.pdf")
    path("volcano_iso_DGE.pdf")

    script:
    """
    IsoformSwitchAnalyzeR_volcano.R --switchlist $rds
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

    run_cpat(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, final_fasta)
    convert_cpat_format(run_cpat.out.ORF_prob_best_tsv)

    split_FASTA(translation_fasta)
    pfam_scan(split_FASTA.out.flatten(), pfamdb)
    merge_pfam_results(pfam_scan.out.collect())
    convert_pfam_scan_results(merge_pfam_results.out)

    final_expression
        .combine(orfanage_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(final_fasta)
        .combine(annotation_gtf)
        .combine(convert_cpat_format.out)
        .combine(convert_pfam_scan_results.out)
    | IsoseqsSwitchList

    PlotIsoformConsequences(IsoseqsSwitchList.out.rds)
    volcano_plot(IsoseqsSwitchList.out.rds)
}
