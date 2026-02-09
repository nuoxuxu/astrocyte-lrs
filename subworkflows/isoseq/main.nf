process fofn {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/fofn"
    input:
    path(flnc_bam)

    output:
    path("*.fofn")

    script:
    """
    fofn.py
    """
}

process cluster {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/cluster"
    input:
    path fofn
    path flnc_bam

    output:
    path("${fofn.simpleName}.clustered.bam")

    script:
    """
    isoseq cluster2 \\
    $fofn \\
    ${fofn.simpleName}.clustered.bam
    """
}

process pbmm2 {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/pbmm2"
    input:
    path ref_genome_fasta
    path clustered_bam

    output:
    path("${clustered_bam.simpleName}.clustered.aligned.bam"), emit: aligned_bam
    path("${clustered_bam.simpleName}.clustered.aligned.bam.bai"), emit: aligned_bam_bai

    script:
    """
    pbmm2 align --preset ISOSEQ --sort \\
    $ref_genome_fasta \\
    $clustered_bam \\
    "${clustered_bam.baseName}.aligned.bam"
    """
}

process merge_aligned_bams {
    conda "/scratch/nxu/SQANTI3/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/merge_aligned_bams"
    input:
    path aligned_bams
    output:
    path("merged.bam")
    script:
    """
    samtools merge -o merged.bam *.aligned.bam
    """
}

process collapse_and_sort {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/discover/isoseq"
    input:
    path merged_bam

    output:
    path("${merged_bam.baseName}.sorted.collapsed.gff")

    script:
    """
    isoseq collapse \\
    --do-not-collapse-extra-5exons \\
    $merged_bam \\
    "${merged_bam.baseName}.collapsed.gff"

    pigeon sort -o "${merged_bam.baseName}.sorted.collapsed.gff" "${merged_bam.baseName}.collapsed.gff"
    """
}



workflow ISOSEQ {
    take:
    flnc_bam
    ref_genome_fasta

    main:
    flnc_bam.map { id, file -> file }.collect().set { flnc_bams }
    fofn(flnc_bams)
    cluster(fofn.out.flatten(), flnc_bams.collect())
    pbmm2(ref_genome_fasta, cluster.out)
    merge_aligned_bams(pbmm2.out.aligned_bam.collect())
    collapse_and_sort(merge_aligned_bams.out)

    emit:
    merged_sorted_collapsed_gtf = collapse_and_sort.out
}