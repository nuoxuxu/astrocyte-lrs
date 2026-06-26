process runORFanage {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    tuple val(orfanage_mode), path(ref_genome_fasta), path(final_sample_gtf), path(annotation_gtf)

    output:
    tuple val(orfanage_mode), path("orfanage_with_gene_id.gtf"), emit: orfanage_gtf
    path "orfanage.stats"

    script:
    def minlen_opt = orfanage_mode == 'no_minlen' ? '' : '--minlen 50'
    """
    orfanage \\
        --reference $ref_genome_fasta \\
        --query $final_sample_gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        $minlen_opt \\
        --stats orfanage.stats \\
        $annotation_gtf

    gffread \\
        -g $ref_genome_fasta \\
        --adj-stop \\
        -T -F -J -C \\
        -o orfanage_with_gene_id.gtf \\
        orfanage_without_gene_id.gtf
    """
}

process addNoncodingTx {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    tuple val(orfanage_mode), path(orfanage_gtf), path(final_transcripts_gtf)

    output:
    tuple val(orfanage_mode), path("complete_orfanage.gtf")

    script:
    """
    add_noncoding_tx_to_orfanage_gtf.py \\
        -i $final_transcripts_gtf \\
        -r $orfanage_gtf \\
        -o complete_orfanage.gtf
    """
}

process fixORFanageFormat {
    label "short_slurm_job"
    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    path ref_genome_fasta
    tuple val(orfanage_mode), path(complete_orfanage_gtf)

    output:
    tuple val(orfanage_mode), path("agat_output.gtf")

    script:
    """
    agat_sp_add_start_and_stop.pl --gff $complete_orfanage_gtf --fasta $ref_genome_fasta --out "added_codons_orfanage_with_gene_id.gff3"

    agat_convert_sp_gff2gtf.pl --gff "added_codons_orfanage_with_gene_id.gff3" -o "agat_output.gtf" --gtf_version 3
    """
}

process restoreAgatRemovedTx {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    tuple val(orfanage_mode), path(complete_orfanage_gtf), path(agat_output_gtf)

    output:
    tuple val(orfanage_mode), path("orfanage.gtf")

    script:
    """
    restore_agat_removed_tx.py \\
        -i $complete_orfanage_gtf \\
        -r $agat_output_gtf \\
        -o "orfanage.gtf" \\
        -v
    """
}

process translateORFs {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    path ref_genome_fasta
    tuple val(orfanage_mode), path(orfanage_gtf)

    output:
    tuple val(orfanage_mode), path("orfanage_proteins.fasta")

    script:
    """
    gffread -y orfanage_proteins.fasta \\
        -g $ref_genome_fasta \\
        $orfanage_gtf
    """
}

process format_gtf_for_ribotie {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage/${orfanage_mode}"

    input:
    tuple val(orfanage_mode), path(orfanage_gtf), path(final_classification), path(annotation_gtf)

    output:
    tuple val(orfanage_mode), path("orfanage_numbered_exons.gtf")

    script:
    """
    format_gtf_for_ribotie.py \\
    --input_gtf ${orfanage_gtf} \\
    --final_classification ${final_classification} \\
    --annotation_gtf ${annotation_gtf} \\
    --output_gtf orfanage_numbered_exons.gtf
    """
}

workflow RUN_ORFANAGE {
    take:
    ref_genome_fasta
    final_transcripts_gtf
    final_classification
    annotation_gtf

    main:
    channel.from(["no_minlen", "minlen"])
        .combine(ref_genome_fasta)
        .combine(final_transcripts_gtf)
        .combine(annotation_gtf)
        | runORFanage

    runORFanage.out.orfanage_gtf
        .combine(final_transcripts_gtf)
        | addNoncodingTx

    fixORFanageFormat(ref_genome_fasta, addNoncodingTx.out)

    addNoncodingTx.out
        .join(fixORFanageFormat.out)
        | restoreAgatRemovedTx

    translateORFs(ref_genome_fasta, restoreAgatRemovedTx.out)

    restoreAgatRemovedTx.out
        .combine(final_classification)
        .combine(annotation_gtf)
        | format_gtf_for_ribotie

    emit:
    orfanage_gtf = format_gtf_for_ribotie.out
    orfanage_proteins = translateORFs.out
}
