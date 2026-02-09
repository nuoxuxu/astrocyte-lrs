process extract_transcriptome {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "${params.outdir}/prepare/extract_transcriptome"

    input:
    path filtered_gtf

    output:
    path "joint_transcriptome.fasta"
    script:
    """
    gffread -w joint_transcriptome.fasta \\
    -g  ${params.ref_genome_fasta} \\
    ${filtered_gtf}
    """
}

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

process minimap2_transcriptome {
    conda "/scratch/nxu/astrocytes/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/align/minimap2_transcriptome"

    input:
    path extracted_transcriptome
    path fastqz_reads

    output:
    path("${fastqz_reads.simpleName}.aligned.bam")

    script:
    """
    minimap2 --eqx -N 100 -ax map-hifi \\
        -t ${task.cpus} $extracted_transcriptome \\
        $fastqz_reads | \\
        samtools sort -n -@ ${task.cpus} \\
        -o "${fastqz_reads.simpleName}.aligned.bam"
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

process run_oarfish {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quantify/isoseq/oarfish"

    input:
    path transcriptome_alignment
    
    output:
    path("${transcriptome_alignment.simpleName}.quant")

    script:
    """
    oarfish --threads ${task.cpus} \\
        --filter-group no-filters \\
        --model-coverage \\
        --alignments $transcriptome_alignment \\
        --output "${transcriptome_alignment.simpleName}"
    """
}

workflow RUN_OARFISH {
    take:
    filtered_gtf
    flnc_bam
    
    main:
    extract_transcriptome(filtered_gtf)
    flnc_bam.map { id, file -> file }.collect().set { flnc_bams }
    merge_flnc_bams(flnc_bams.collect())
    convert_flnc_bam_to_fastqz(merge_flnc_bams.out.flatten())
    minimap2_transcriptome(extract_transcriptome.out, convert_flnc_bam_to_fastqz.out)
    run_oarfish(minimap2_transcriptome.out)    
    
    emit:
    oarfish_quant = run_oarfish.out
}