process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}"

    input:
    tuple val(version), path(final_expression), path(predicted_cds_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    path("IsoformSwitchAnalyzeR.rds"), emit: rds
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

process run_cpat {
    module "python:gcc:arrow/19.0.1:rust:r/4.4.0"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/cpat"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/pfam/split_fasta"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/pfam"

    input:
    tuple val(version), path(translation_fasta), path(pfamdb)

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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/cpat"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/cpat"

    input:
    tuple val(version), path(cpat_results), path(ribotie_fasta)

    output:
    path("cpat_results_filtered.txt")

    script:
    """
    filter_cpat_results.py $cpat_results $ribotie_fasta -o cpat_results_filtered.txt
    """
}

process merge_pfam_results {
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/pfam"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/pfam"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/figures"

    input:
    val(version)
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
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}/figures"

    input:
    val(version)
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

process prepare_shinyApp {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "astrocyte_vis_app/data/${version}"

    input:
    val(version)
    path(rds)

    output:
    path("switchPlotFromTables.RData")

    script:
    """
    extract_rds_data.R --rds $rds --output switchPlotFromTables.RData
    """
}

workflow ISOFORMSWITCH {
    take:
    version
    final_expression
    primer_to_sample
    final_fasta
    predicted_cds_gtf
    annotation_gtf
    final_classification

    main:

    final_expression
        .combine(predicted_cds_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(final_fasta)
        .combine(annotation_gtf)
        .combine(version)
        .map { expr, cds_gtf, classif, primer, fasta, annot, ver ->
            tuple(ver, expr, cds_gtf, classif, primer, fasta, annot)
        }
    | IsoseqsSwitchList

    PlotIsoformConsequences(version, IsoseqsSwitchList.out.rds)
    volcano_plot(version, IsoseqsSwitchList.out.rds)
    prepare_shinyApp(version, IsoseqsSwitchList.out.rds)

    emit:
    isoform_features_csv = IsoseqsSwitchList.out.isoform_features_csv
}
