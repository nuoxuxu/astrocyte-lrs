process sqanti_qc {
    label "mid_slurm_job"
    container "sqanti3_latest.sif"
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