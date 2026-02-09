process merge_flnc_bams {
    conda "/scratch/nxu/SQANTI3/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/prepare/merge_flnc_bams"
    
    input:
    path(flnc_bam)

    output:
    path("*.flnc.bam")

    script:
    """
    merge_flnc_bam.sh
    """
}

process convert_flnc_bam_to_fastqz {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/convert_ubam_to_fastqz"

    input:
    path merged_flnc_bam
    output:
    path("${merged_flnc_bam.simpleName}.fastq.gz")

    script:
    """
    samtools fastq -@ ${task.cpus} -0 ${merged_flnc_bam.simpleName}.fastq.gz \\
    -c ${params.compression_level} $merged_flnc_bam
    """
}

process minimap2_genome {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/minimap2_genome"    
    input:
    path fastqz
    output:
    path("${fastqz.simpleName}.aligned.bam")
    script:
    """
    minimap2 -ax splice:hq -t ${task.cpus} -uf ${params.ref_genome_fasta} $fastqz | samtools sort -@ ${task.cpus} \\
        -o "${fastqz.simpleName}.aligned.bam"
    """
}

process bambu {
    conda "/scratch/nxu/astrocytes/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/discover/bambu"
    input:
    path minimap_bam
    output:
    path("supportedTranscriptModels.gtf"), emit: supported_tx_gtf
    path("novelTranscripts.gtf"), emit: novel_tx_gtf
    script:
    """
    run_bambu.R ${params.annotation_gtf} "*minimap.bam" ${params.ref_genome_fasta} ${task.cpus}
    """
}

workflow BAMBU {
    take:
    flnc_bam

    main:
    flnc_bam.map { id, file -> file }.collect().set { flnc_bams }
    merge_flnc_bams(flnc_bams.collect())
    convert_flnc_bam_to_fastqz(merge_flnc_bams.out.flatten())
    minimap2_genome(convert_flnc_bam_to_fastqz.out)    
    bambu(minimap2_genome.out.collect())

    emit:
    supported_tx_gtf = bambu.out.supported_tx_gtf
    novel_tx_gtf = bambu.out.novel_tx_gtf
}