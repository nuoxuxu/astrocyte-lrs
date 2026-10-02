process supplement_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/supplemented_collaborator"

    input:
    path(collaborator_gtf)
    path(orfanage_gtf)

    output:
    path("supplemented_collaborator.gtf"), emit: supplemented_gtf

    script:
    """
    supplement_collaborator_gtf.py $collaborator_gtf $orfanage_gtf -o supplemented_collaborator.gtf
    """
}

process prepare_supplemented_files {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/supplemented_collaborator"

    input:
    tuple path(supplemented_gtf), path(final_fasta), path(final_expression)

    output:
    path("supplemented_collaborator.fasta"),                    emit: supplemented_fasta
    path("supplemented_collaborator_expression.parquet"),       emit: supplemented_expression

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    prepare_supplemented_gtf_files.py \\
        $supplemented_gtf \\
        $final_fasta \\
        $final_expression \\
        --output_fasta supplemented_collaborator.fasta \\
        --output_expression supplemented_collaborator_expression.parquet
    """
}

process translate_supplemented_ORFs {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/supplemented_collaborator"

    input:
    path ref_genome_fasta
    path supplemented_gtf

    output:
    path("supplemented_collaborator_proteins.fasta"), emit: supplemented_proteins

    script:
    """
    gffread -y supplemented_collaborator_proteins.fasta \\
        -g $ref_genome_fasta \\
        $supplemented_gtf
    """
}

process filter_ribotie_for_isoformswitch {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/${name}"

    input:
    tuple val(name), path(input_gtf), path(input_fasta), path(input_expression), path(input_classification)

    output:
    tuple val(name), path("${name}_filtered.fasta"),                            emit: fasta
    tuple val(name), path("${name}_filtered_expression.parquet"),               emit: expression
    tuple val(name), path("${name}_filtered_classification.parquet"),           emit: classification

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_RiboTIE_results.py \\
        --input_gtf $input_gtf \\
        --input_fasta $input_fasta \\
        --input_expression $input_expression \\
        --input_classification $input_classification \\
        --output_fasta ${name}_filtered.fasta \\
        --output_expression ${name}_filtered_expression.parquet \\
        --output_classification ${name}_filtered_classification.parquet
    """
}

workflow SUPPLEMENT_COLLABORATOR_ORF {
    take:
    fixed_collaborator_gtf
    orfanage_gtf
    final_fasta
    final_expression
    ref_genome_fasta

    main:
    // Supplement collaborator GTF with ORFs from ORFanage to create a comprehensive set of ORFs for analysis
    supplement_collaborator_gtf(fixed_collaborator_gtf, orfanage_gtf)

    // Prepare supplemented files for downstream analysis
    supplement_collaborator_gtf.out.supplemented_gtf
        .combine(final_fasta)
        .combine(final_expression)
    | prepare_supplemented_files

    // Translate supplemented ORFs to a format suitable for downstream analysis
    translate_supplemented_ORFs(ref_genome_fasta, supplement_collaborator_gtf.out.supplemented_gtf)

    emit:
    supplemented_gtf = supplement_collaborator_gtf.out.supplemented_gtf
    supplemented_fasta = prepare_supplemented_files.out.supplemented_fasta
    supplemented_expression = prepare_supplemented_files.out.supplemented_expression
    supplemented_proteins = translate_supplemented_ORFs.out.supplemented_proteins    
}