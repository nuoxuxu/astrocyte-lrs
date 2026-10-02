process sqanti_qc {
    label "mid_slurm_job"
    container "sqanti3_latest.sif"
    storeDir "nextflow_results/sqanti3/isoseq"
    input:
    path gff_file
    path annotation_gtf
    path ref_genome_fasta
    path refTSS
    path polyA_motif_list
    path star_aligned_bam
    path star_sj_tab

    output:
    path("sqanti3_qc")

    script:
    """
    export PATH=/conda/miniconda3/envs/sqanti3/bin:\$PATH
    find . -name "*.bam" > SR_bam.fofn

    sqanti3_qc.py \\
    --isoforms $gff_file \\
    --refGTF $annotation_gtf \\
    --refFasta $ref_genome_fasta \\
    --CAGE_peak $refTSS \\
    --polyA_motif_list $polyA_motif_list \\
    --report html \\
    --skipORF \\
    --SR_bam SR_bam.fofn \\
    -c "*.SJ.out.tab" \\
    -o sqanti_qc_results \\
    -d sqanti3_qc \\
    -n 10
    """
}

process sqanti_filter {
    conda "/scratch/nxu/SQANTI3/env"
    label "short_slurm_job"
    storeDir "nextflow_results/sqanti3/isoseq"

    input:
    path corrected_gtf
    path classification

    output:
    path("sqanti3_filter")

    script:
    """
    sqanti3_filter.py rules \\
    --filter_gtf $corrected_gtf \\
    --sqanti_class $classification \\
    -d sqanti3_filter \\
    -o default
    """
}

workflow SQANTI {
    take:
    annotation_gtf
    ref_genome_fasta
    refTSS
    polyA_motif_list
    merged_sorted_collapsed_gtf
    star_aligned_bam
    star_sj_tab

    main:
    sqanti_qc(merged_sorted_collapsed_gtf, annotation_gtf, ref_genome_fasta, refTSS, polyA_motif_list, star_aligned_bam.collect(), star_sj_tab.collect())
    isoseq_corrected_gtf = sqanti_qc.out
        .map { dir -> dir / "sqanti_qc_results_corrected.gtf" }
    isoseq_classification = sqanti_qc.out
        .map { dir -> dir / "sqanti_qc_results_classification.txt" }
    sqanti_filter(isoseq_corrected_gtf, isoseq_classification)
    filtered_gtf = sqanti_filter.out
        .map { dir -> dir / "default.filtered.gtf"}
    filtered_classification = sqanti_filter.out
        .map { dir -> dir / "default_RulesFilter_result_classification.txt" }
    sqanti_corrected_fasta = sqanti_qc.out
        .map { dir -> dir / "sqanti_qc_results_corrected.fasta" }

    emit:
    filtered_gtf = filtered_gtf
    filtered_classification = filtered_classification
    sqanti_corrected_fasta = sqanti_corrected_fasta
}
