process SQANTI_PROTEIN {
    label 'short_slurm_job'
    container "sqanti3_latest.sif"
    storeDir "nextflow_results/sqanti3_protein"

    input:
    path cds_gtf
    path reference_gtf
    path sqanti_protein_script

    output:
    path("*.predicted_proteome.best_ORF_SQANTI_classification.tsv"), emit: protein_classification
    path("S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    exec > >(tee S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt) 2>&1
    source /conda/miniconda3/etc/profile.d/conda.sh
    conda activate sqanti3
    export SQANTI_PATH=\$(dirname \$(which sqanti3_qc.py))

    # Add src/utilities and utilities to PYTHONPATH for cupcake and other imports
    export PYTHONPATH=\${SQANTI_PATH}/src/utilities:\${SQANTI_PATH}/utilities:\${SQANTI_PATH}:\${PYTHONPATH:-}

    # Copy the script locally and patch it to use the platform-specific gtfToGenePred binary
    # v5.5.4 has gtfToGenePred-linux-x86_64 instead of gtfToGenePred
    cp $sqanti_protein_script ./sqanti3_protein_input_full_gtf_patched.py
    sed -i 's|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred")|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred-linux-x86_64")|g' ./sqanti3_protein_input_full_gtf_patched.py

    python ./sqanti3_protein_input_full_gtf_patched.py \\
        $cds_gtf \\
        $reference_gtf \\
        -d . \\
        -p test

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.5.4
    END_VERSIONS
    """
}

process PROTEIN_CLASSIFICATION {
    tag "$meta"
    label 'short_slurm_job'
    container "lrp2-lite_latest.sif"
    storeDir "nextflow_results/protein_classification"

    input:
    tuple val(meta), path(protein_classification), path(cds_gtf), path(corrected_fasta), path(all_orfs_mapped)
    path reference_gtf
    path protein_class_script

    output:
    tuple val(meta), path("*.predicted_proteome.best_ORF_summary.txt"), emit: protein_all_isoforms
    tuple val(meta), path("*.predicted_proteome.best_ORF.fa"), emit: protein_all_orfs_fasta
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF_hashids.txt"), emit: hashids_orf
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF.gtf"), emit: protein_gtf
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF.bed"), emit: protein_bed
    tuple val(meta), path("*_S3_PREDICTED_PROTEOME_M4_PROTEIN_CLASSIFICATION_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def min_junctions_after_stop = task.ext.min_junctions_after_stop ?: params.min_junctions_after_stop_codon
    def protein_class_keep = task.ext.protein_class_keep ?: params.protein_class_keep

    """
    exec > >(tee ${prefix}_S3_PREDICTED_PROTEOME_M4_PROTEIN_CLASSIFICATION_log.txt) 2>&1

    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript \$(pwd)/$protein_class_script \\
        --basename $prefix \\
        --gencode_gtf $reference_gtf \\
        --sample_cds_gtf $cds_gtf \\
        --sample_dna_fasta $corrected_fasta \\
        --mapped_orfs $all_orfs_mapped \\
        --protein_sqanti $protein_classification \\
        --output_dir . \\
        --min_junctions_after_stop $min_junctions_after_stop \\
        --protein_class_keep "$protein_class_keep" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.best_ORF_summary.txt
    touch ${prefix}.predicted_proteome.best_ORF.fa
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF.gtf
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF.bed
    touch ${prefix}_S3_PREDICTED_PROTEOME_M4_PROTEIN_CLASSIFICATION_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}

workflow LRP2_PROTEIN {
    take:
    name              // val: sample/run label used as PROTEIN_CLASSIFICATION prefix
    fixed_collaborator_gtf  // CDS GTF used for protein classification (fixed collaborator GTF)
    corrected_fasta   // transcript FASTA (final_transcripts.fasta)
    all_orfs_mapped   // RiboTIE ORF table (ribotie_cpm1_3sample.csv)
    annotation_gtf    // GENCODE reference GTF

    main:
    SQANTI_PROTEIN(fixed_collaborator_gtf, annotation_gtf, file("${projectDir}/bin/sqanti3_protein.py"))

    PROTEIN_CLASSIFICATION(
        channel.value(name)
            .combine(SQANTI_PROTEIN.out.protein_classification)
            .combine(fixed_collaborator_gtf)
            .combine(corrected_fasta)
            .combine(all_orfs_mapped),
        annotation_gtf,
        file("${projectDir}/bin/protein_classification.R")
    )

    emit:
    sqanti_protein_classification = SQANTI_PROTEIN.out.protein_classification
    protein_all_isoforms          = PROTEIN_CLASSIFICATION.out.protein_all_isoforms
    protein_all_orfs_fasta        = PROTEIN_CLASSIFICATION.out.protein_all_orfs_fasta
    hashids_orf                   = PROTEIN_CLASSIFICATION.out.hashids_orf
    protein_gtf                   = PROTEIN_CLASSIFICATION.out.protein_gtf
    protein_bed                   = PROTEIN_CLASSIFICATION.out.protein_bed
}
