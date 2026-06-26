process extract_protein_sequences {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/cds_length_distribution/${orfanage_mode}"

    input:
    tuple val(orfanage_mode), path(ribotie_gtf), path(ref_genome_fasta)

    output:
    tuple val(orfanage_mode), path("ribotie_proteins_${orfanage_mode}.fasta")

    script:
    """
    gffread -y ribotie_proteins_${orfanage_mode}.fasta \\
        -g $ref_genome_fasta \\
        $ribotie_gtf
    """
}

process plot_cds_length_distribution {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/cds_length_distribution"

    input:
    tuple val(orfanage_modes), path(protein_fastas)

    output:
    path("cds_length_distribution.png")

    script:
    """
    plot_cds_length_distribution.py \\
        --fastas ${protein_fastas.join(' ')} \\
        --modes ${orfanage_modes.join(' ')} \\
        --output_png cds_length_distribution.png
    """
}

workflow CDS_LENGTH_DISTRIBUTION {
    take:
    ribotie_gtf_ch
    ref_genome_fasta

    main:
    ribotie_gtf_ch
        .combine(ref_genome_fasta)
        | extract_protein_sequences

    extract_protein_sequences.out
        .toList()
        .map { items -> [items.collect { it[0] }, items.collect { it[1] }] }
        | plot_cds_length_distribution

    emit:
    cds_length_plot = plot_cds_length_distribution.out
}
