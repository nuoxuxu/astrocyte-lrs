
process run_cpat {
    module "python:gcc:arrow/19.0.1:rust:r/4.4.0"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/cpat"

    input:
    path(Human_coding_transcripts_CDS)
    path(Human_noncoding_transcripts_RNA)
    path(Human_logitModel)
    path(filtered_RiboTIE_fasta)

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
        -g $filtered_RiboTIE_fasta \\
        --min-orf=50 \\
        --top-orf=50 \\
        -o CPAT \\
        1> CPAT.output \\
        2> CPAT.error
    """
}

process split_FASTA {
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/pfam/split_fasta"

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
    storeDir "nextflow_results/summary_table/pfam"

    input:
    tuple path(translation_fasta), path(pfamdb)

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
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/cpat"

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
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/cpat"

    input:
    tuple path(cpat_results), path(ribotie_fasta)

    output:
    path("cpat_results_filtered.txt")

    script:
    """
    filter_cpat_results.py $cpat_results $ribotie_fasta -o cpat_results_filtered.txt
    """
}

process merge_pfam_results {
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/pfam"

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
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/pfam"

    input:
    path(pfam_scan_results)

    output:
    path("pfam_scan_results_modified.txt")

    script:
    """
    convert_pfam_csv_to_txt.py $pfam_scan_results pfam_scan_results_modified.txt
    """
}

process make_summary_table {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table"

    input:
    path(ribotie_csv)
    path(collaborator_gtf)
    path(phylocsfpp_gtf)
    path(pfam_csv)
    path(annotation_gtf)
    path(orfanage_gtf)
    path(final_classification)
    path(isoform_features_csv)
    path(study2_gtf)
    path(leafcutter_coloc)
    path(novel_coding_junction_coloc)
    path(concordant_csv)

    output:
    path("summary_table.tsv")

    script:
    """
    make_summary_table.py $ribotie_csv $collaborator_gtf $phylocsfpp_gtf $pfam_csv $annotation_gtf $orfanage_gtf $final_classification $isoform_features_csv $study2_gtf $leafcutter_coloc $novel_coding_junction_coloc $concordant_csv -o summary_table.tsv
    """
}

process filter_summary_table {
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table"

    input:
    path(summary_table)

    output:
    path("summary_table_isoform_switch.tsv")

    script:
    """
    awk -F'\\t' 'NR==1 || (\$15==1 && \$16==1 && \$17==1)' $summary_table > summary_table_isoform_switch.tsv
    """
}

process phylocsfpp {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/summary_table/phylocsfpp"

    input:
    path(ribotie_output_gtf)
    path(phyloCSF_db)

    output:
    path("${ribotie_output_gtf.baseName}.PhyloCSF++.gtf")

    script:
    """
    phylocsf++ annotate-with-tracks $phyloCSF_db/PhyloCSF+1.bw $ribotie_output_gtf
    """
}

workflow SUMMARY_TABLE {
    take:
    Human_coding_transcripts_CDS
    Human_noncoding_transcripts_RNA
    Human_logitModel
    filtered_RiboTIE_fasta
    filtered_RiboTIE_proteins
    pfamdb
    collaborator_gtf
    ribotie_cpm1_3sample
    annotation_gtf
    orfanage_gtf
    final_classification
    isoform_features_csv
    study2_gtf
    leafcutter_coloc
    novel_coding_junction_coloc
    concordant_csv

    main:
    channel.value(file(params.PhyloCSFpp_db)).set { phyloCSF_db }

    run_cpat(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, filtered_RiboTIE_fasta)
    convert_cpat_format(run_cpat.out.ORF_prob_best_tsv)

    split_FASTA(filtered_RiboTIE_proteins)
    split_FASTA.out.flatten()
        .combine(pfamdb)
    | pfam_scan
    merge_pfam_results(pfam_scan.out.collect())
    convert_pfam_scan_results(merge_pfam_results.out)
    phylocsfpp(collaborator_gtf, phyloCSF_db)
    make_summary_table(ribotie_cpm1_3sample, collaborator_gtf, phylocsfpp.out, merge_pfam_results.out, annotation_gtf, orfanage_gtf, final_classification, isoform_features_csv, study2_gtf, leafcutter_coloc, novel_coding_junction_coloc, concordant_csv)
    filter_summary_table(make_summary_table.out)

    emit:
    cpat_results = convert_cpat_format.out
    pfam_scan_results = convert_pfam_scan_results.out
    summary_table = make_summary_table.out
    summary_table_isoform_switch = filter_summary_table.out
}
