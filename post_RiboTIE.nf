include { AIM_2 } from "./subworkflows/local/aim_2/main.nf"
include { ISOFORMSWITCH } from "./subworkflows/local/IsoformSwitchAnalyzeR/main.nf"
include { RIBOTIE_POSTANALYSIS } from "./subworkflows/local/ribotie_postanalysis/main.nf"
include { FILTER_RIBOTIE } from "./subworkflows/local/filter_ribotie/main.nf"
include { CDS_LENGTH_DISTRIBUTION } from "./subworkflows/local/cds_length_distribution/main.nf"
include { RUN_VEP } from "./subworkflows/local/vep/main.nf"
include { SUMMARY_TABLE } from "./subworkflows/local/summary_table/main.nf"
include { LABEL_ORF_TYPE_GENCODE } from "./subworkflows/local/quality" 

process fix_collaborator_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/translatome/supplemented_collaborator"

    input:
    path(collaborator_gtf)

    output:
    path("filtered_output_fixed.gtf"), emit: fixed_gtf

    script:
    """
    fix_collaborator_gtf.py $collaborator_gtf -o filtered_output_fixed.gtf
    """
}

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

process IsoseqsSwitchList {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/IsoformSwitchAnalyzeR/${version}"

    input:
    tuple val(version), path(final_expression), path(predicted_cds_gtf), path(final_classification), path(primer_to_sample), path(final_fasta), path(annotation_gtf)

    output:
    tuple val(version), path("IsoformSwitchAnalyzeR.rds"), emit: rds
    path("isoformFeatures.csv"), emit: isoform_features_csv

    script:
    """
    IsoformSwitchAnalyzeR.R \\
        --final_expression $final_expression \\
        --primer_to_sample $primer_to_sample \\
        --final_fasta $final_fasta \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --annotation_gtf $annotation_gtf \\
        --final_classification $final_classification \\
    """
}

process label_orf_type_gencode {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/${meta}"

    input:
    tuple val(meta), path(ribotie_csv), path(annotation_gtf), path(orfanage_gtf), path(final_classification)

    output:
    path("orf_type_gencode.tsv"), emit: orf_type_gencode

    script:
    """
    export POLARS_MAX_THREADS=1
    label_orf_type_gencode.py \\
        $ribotie_csv \\
        $annotation_gtf \\
        $orfanage_gtf \\
        $final_classification \\
        -o orf_type_gencode.tsv
    """
}

process add_lncRNA_to_collaborator_csv {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/collaborator"

    input:
    path(ribotie_csv)
    path(final_classification)
    path(annotation_gtf)

    output:
    path("ribotie_cpm1_3sample_with_lncRNA.csv"), emit: ribotie_csv_with_lncRNA

    script:
    """
    add_lncRNA.py $ribotie_csv $final_classification $annotation_gtf -o ribotie_cpm1_3sample_with_lncRNA.csv
    """
}

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

process peptideTrackUCSC {
    storeDir "nextflow_results/proteomic"
    conda "/scratch/nxu/astrocytes/env"
    
    input:
    val mode
    path annotation_gtf
    path final_sample_classification
    path predicted_cds_gtf
    path protein_search_database
    path peptides

    output:
    path "peptides_${mode}.gtf"

    script:
    """
    make_peptide_gtf_file.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $final_sample_classification \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --protein_search_database $protein_search_database \\
        --peptides $peptides \\
        --output "peptides_${mode}.gtf"
    """
}

process peptideMapping {
    storeDir "nextflow_results/proteomic"
    conda "/scratch/nxu/astrocytes/env"

    input:
    path annotation_gtf
    path classification
    path protein_search_database
    path percolator_res
    
    output:
    path "peptide_mapping.parquet", emit: peptide_mapping
    path "novel_peptides.csv", emit: novel_peptides
    
    script:
    """
    peptide_mapping.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $classification \\
        --protein_search_database $protein_search_database \\
        --percolator_res $percolator_res
    """
}


workflow {
    // Data flow: orfanage ORFs → RiboTIE scoring (by collaborator) → filtered_output.gtf (high-confidence RiboTIE hits)
    // supplement_collaborator_gtf re-adds orfanage ORFs filtered out by RiboTIE, creating a comprehensive ORF set

    channel.value(file(params.annotation_gtf)).set { annotation_gtf }
    channel.value(file(params.primer_to_sample)).set { primer_to_sample }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_transcripts.fasta")).set { final_fasta }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_expression.parquet")).set { final_expression }
    channel.value(file("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_classification.parquet")).set { final_classification }
    channel.value(file(params.ref_genome_fasta)).set { ref_genome_fasta }
    channel.value(file("from_collaborator/filtered_output.gtf")).set { collaborator_gtf }
    channel.value(file("from_collaborator/ribotie_cpm1_3sample.csv")).set { ribotie_cpm1_3sample }
    channel.value(file("from_collaborator/concordant_all_three_exp.csv")).set { concordant_csv }
    channel.value(file(params.bigbrain_sqtl)).set { bigbrain_sqtl }
    channel.value(file(params.bigbrain_coloc)).set { bigbrain_coloc }
    channel.value(file(params.leafcutter_sig)).set { leafcutter_sig }
    channel.value(file(params.leafcutter_clu2gene)).set { leafcutter_clu2gene }
    channel.value(file("nextflow_results/orfanage/minlen/orfanage.gtf")).set { orfanage_gtf }
    channel.value(file(params.Human_coding_transcripts_CDS)).set { Human_coding_transcripts_CDS }
    channel.value(file(params.Human_noncoding_transcripts_RNA)).set { Human_noncoding_transcripts_RNA }
    channel.value(file(params.Human_logitModel)).set { Human_logitModel }
    channel.value(file(params.pfamdb)).set { pfamdb }
    channel.value(file("data/study2_orfs.gtf")).set { study2_gtf }

    fix_collaborator_gtf(collaborator_gtf)

    supplement_collaborator_gtf(fix_collaborator_gtf.out.fixed_gtf, orfanage_gtf)

    supplement_collaborator_gtf.out.supplemented_gtf
        .combine(final_fasta)
        .combine(final_expression)
    | prepare_supplemented_files

    translate_supplemented_ORFs(ref_genome_fasta, supplement_collaborator_gtf.out.supplemented_gtf)

    // Build per-version channel from manifests (exclude gencode)
    channel.fromPath(params.main_pipeline_outputs)
        .map { f ->
            def name = f.baseName.replaceAll(/^main_pipeline_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.final_expression), file(entry.final_fasta), file(entry.final_classification))
        }
        .set { main_outputs_ch }

    channel.fromPath(params.ribotie_training_outputs)
        .filter { f -> f.baseName.replaceAll(/^ribotie_training_outputs_/, '') != 'gencode' }
        .map { f ->
            def name = f.baseName.replaceAll(/^ribotie_training_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.ribotie_merged_gtf))
        }
        .join(main_outputs_ch)
        .map { name, ribotie_gtf, expr, fasta, classif ->
            tuple(name, ribotie_gtf, fasta, expr, classif)
        }
        .set { ribotie_inputs_ch }

    // Filter RiboTIE GTF to match expression data (versioned)
    filter_ribotie_for_isoformswitch(ribotie_inputs_ch)

    // Create a channel with the isoform information
    ribotie_inputs_ch.map { it.take(2) }
        .join(filter_ribotie_for_isoformswitch.out.expression)
        .map { version, ribotie_gtf, filtered_expression ->
            tuple(version, filtered_expression, ribotie_gtf)
        }
        .join(filter_ribotie_for_isoformswitch.out.classification)
        .combine(primer_to_sample)
        .join(filter_ribotie_for_isoformswitch.out.fasta)
        .combine(annotation_gtf)
        .set { my_isoform_ch }

    channel.value("supplemented_collaborator")
        .combine(prepare_supplemented_files.out.supplemented_expression)
        .combine(supplement_collaborator_gtf.out.supplemented_gtf)
        .combine(final_classification)
        .combine(primer_to_sample)
        .combine(prepare_supplemented_files.out.supplemented_fasta)
        .combine(annotation_gtf)
        .set { collaborator_isoform_ch }

    my_isoform_ch.mix(collaborator_isoform_ch) | ISOFORMSWITCH

    channel.value("collaborator")
        .combine(ribotie_cpm1_3sample)
        .combine(annotation_gtf)
        .combine(orfanage_gtf)
        .combine(final_classification) |
        label_orf_type_gencode

    add_lncRNA_to_collaborator_csv(ribotie_cpm1_3sample, final_classification, annotation_gtf)

    peptideTrackUCSC(
        "collaborator", 
        annotation_gtf, 
        final_classification, 
        fix_collaborator_gtf.out.fixed_gtf, 
        file("nextflow_results/ribotie/filtered/filtered_RiboTIE_proteins.fasta"), 
        file("results/proteomics/peptide.tsv")
    )
    peptideMapping(
        annotation_gtf, 
        final_classification, 
        file("nextflow_results/ribotie/filtered/filtered_RiboTIE_proteins.fasta"), 
        file("results/proteomics/peptide.tsv")
    )

    SQANTI_PROTEIN(supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, file("${projectDir}/bin/sqanti3_protein.py"))

    PROTEIN_CLASSIFICATION (
        channel.value("collaborator")
            .combine(SQANTI_PROTEIN.out.protein_classification)
            .combine(fix_collaborator_gtf.out.fixed_gtf)
            .combine(final_fasta)
            .combine(ribotie_cpm1_3sample),
        annotation_gtf,
        file("${projectDir}/bin/protein_classification.R")
    )

    // AIM_2(supplement_collaborator_gtf.out.supplemented_gtf, annotation_gtf, bigbrain_sqtl, bigbrain_coloc, ISOFORMSWITCH.out.isoform_features_csv, leafcutter_sig, leafcutter_clu2gene)

    // FILTER_RIBOTIE(fix_collaborator_gtf.out.fixed_gtf, final_fasta, final_expression, final_classification, ref_genome_fasta)

    // SUMMARY_TABLE(Human_coding_transcripts_CDS, Human_noncoding_transcripts_RNA, Human_logitModel, FILTER_RIBOTIE.out.filtered_RiboTIE_fasta, FILTER_RIBOTIE.out.filtered_RiboTIE_proteins, pfamdb, fix_collaborator_gtf.out.fixed_gtf, ribotie_cpm1_3sample, annotation_gtf, orfanage_gtf, final_classification, ISOFORMSWITCH.out.isoform_features_csv, study2_gtf, AIM_2.out.leafcutter_coloc, AIM_2.out.novel_coding_junction_coloc, concordant_csv)

    //TODO: MAPS analysis for variants disruption ncORFs
    // channel.value(file("/scratch/nxu/100KGP_splicing/data/gnomad/exomes/gnomad.exomes.v4.1.sites.chr16.vcf.bgz")).map { ["gnomad_exomes_chr16", it] }.set { vcf_ch }
    // channel.value(file(params.vep_uorf_data)).set { vep_uorf_data }
    // RUN_VEP(vcf_ch, vep_uorf_data, ref_genome_fasta, annotation_gtf)
}
