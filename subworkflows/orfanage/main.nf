process runORFanage {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    path ref_genome_fasta
    tuple val(param_set_name), path(final_sample_gtf)

    output:
    tuple val(param_set_name), path("orfanage_with_gene_id.gtf"), emit: orfanage_gtf
    path "orfanage.stats"

    script:
    """
    gffread \\
        -g $ref_genome_fasta \\
        --adj-stop \\
        -T -F -J \\
        -o corrected.gtf \\
        $final_sample_gtf

    orfanage \\
        --reference $ref_genome_fasta \\
        --query corrected.gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        --minlen 50 \\
        --stats orfanage.stats \\
        ${params.annotation_gtf}

    gffread \\
        -g $ref_genome_fasta \\
        --adj-stop \\
        -T -F -J -C \\
        -o orfanage_with_gene_id.gtf \\
        orfanage_without_gene_id.gtf
    """
}

process fixORFanageFormat {
    label "short_slurm_job"
    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    path ref_genome_fasta
    tuple val(param_set_name), path(orfanage_gtf)

    output:
    tuple val(param_set_name), path("orfanage.gtf")
    
    script:
    """
    agat_sp_add_start_and_stop.pl --gff $orfanage_gtf --fasta $ref_genome_fasta --out "added_codons_orfanage_with_gene_id.gff3"

    agat_convert_sp_gff2gtf.pl --gff "added_codons_orfanage_with_gene_id.gff3" -o "orfanage.gtf" --gtf_version 3
    """
}

process translateORFs {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    path ref_genome_fasta
    tuple val(param_set_name), path(orfanage_gtf)

    output:
    tuple val(param_set_name), path("orfanage_proteins.fasta")

    script:
    """
    gffread -y orfanage_proteins.fasta \\
        -g $ref_genome_fasta \\
        $orfanage_gtf
    """
}

workflow RUN_ORFANAGE {
    take:
    ref_genome_fasta
    final_transcripts_gtf

    main:
    runORFanage(ref_genome_fasta, final_transcripts_gtf)
    fixORFanageFormat(ref_genome_fasta, runORFanage.out.orfanage_gtf)
    translateORFs(params.ref_genome_fasta, fixORFanageFormat.out)

    emit:
    orfanage_gtf = fixORFanageFormat.out
    orfanage_proteins = translateORFs.out
}