process runORFanage {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    path ref_genome_fasta
    tuple val(param_set_name), path(final_sample_gtf)
    path annotation_gtf

    output:
    tuple val(param_set_name), path("orfanage_with_gene_id.gtf"), emit: orfanage_gtf
    path "orfanage.stats"

    script:
    """
    orfanage \\
        --reference $ref_genome_fasta \\
        --query $final_sample_gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        --minlen 50 \\
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
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    tuple val(param_set_name), path(orfanage_gtf), path(final_transcripts_gtf)

    output:
    tuple val(param_set_name), path("complete_orfanage.gtf")

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
    storeDir "nextflow_results/orfanage/${param_set_name}"

    input:
    path ref_genome_fasta
    tuple val(param_set_name), path(complete_orfanage_gtf)

    output:
    tuple val(param_set_name), path("fixed_complete_orfanage.gtf")
    
    script:
    """
    agat_sp_add_start_and_stop.pl --gff $complete_orfanage_gtf --fasta $ref_genome_fasta --out "added_codons_orfanage_with_gene_id.gff3"

    agat_convert_sp_gff2gtf.pl --gff "added_codons_orfanage_with_gene_id.gff3" -o "agat_output.gtf" --gtf_version 3

    restore_agat_removed_tx.py \\
        -i $complete_orfanage_gtf \\
        -r "agat_output.gtf" \\
        -o "fixed_complete_orfanage.gtf" \\
        -v
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
    annotation_gtf

    main:
    runORFanage(ref_genome_fasta, final_transcripts_gtf, annotation_gtf)
    runORFanage.out.orfanage_gtf
        .join(final_transcripts_gtf)
        | addNoncodingTx
    fixORFanageFormat(ref_genome_fasta, addNoncodingTx.out)
    translateORFs(ref_genome_fasta, fixORFanageFormat.out)

    emit:
    orfanage_gtf = fixORFanageFormat.out
    orfanage_proteins = translateORFs.out
}