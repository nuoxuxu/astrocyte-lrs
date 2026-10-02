process fix_RiboTIE_ORF_type {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "/scratch/nxu/astrocytes/from_collaborator"

    input:
    path(ribotie_config)
    path(gtf)
    path(ribotie_csv)
    path(gencode_h5_files)
    path(orfanage_h5_files)

    output:
    path("ribotie_res_merged_fixed.csv"),           emit: filtered
    path("ribotie_res_merged_fixed.redundant.csv"), emit: redundant
    path("ribotie_res_merged_fixed.novel.csv"),     emit: novel

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    fix_RiboTIE_ORF_type.py \\
        --config $ribotie_config \\
        --gtf $gtf \\
        --ribotie_csv $ribotie_csv \\
        --h5 $gencode_h5_files $orfanage_h5_files \\
        --out_prefix ribotie_res_merged_fixed
    """
}

process add_lncRNA_to_collaborator_csv {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "from_collaborator"

    input:
    path(ribotie_csv)
    path(final_classification)
    path(annotation_gtf)

    output:
    path("ribotie_res_merged_fixed_with_lncRNA.csv"), emit: ribotie_csv_with_lncRNA

    script:
    """
    add_lncRNA.py $ribotie_csv $final_classification $annotation_gtf -o ribotie_res_merged_fixed_with_lncRNA.csv
    """
}

workflow FIX_COLLABORATOR_ORF_TYPE {
    take:
    fixed_collaborator_gtf
    final_classification
    annotation_gtf
    ribotie_cpm1_3sample

    main:
    // Re-evaluate CDS overlap of collaborator RiboTIE ORFs against merged ORFanage + GENCODE CDS annotations
    gencode_h5_ch = channel.fromPath(params.ribotie_training_inputs)
        .filter { f -> f.baseName.replaceAll(/^ribotie_training_inputs_/, '') == 'gencode' }
        .map { f ->
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            file(entry.gtf_h5)
        }

    minlen_h5_ch = channel.fromPath(params.ribotie_training_inputs)
        .filter { f -> f.baseName.replaceAll(/^ribotie_training_inputs_/, '') == 'minlen' }
        .map { f ->
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            file(entry.gtf_h5)
        }

    fix_RiboTIE_ORF_type(
        file("/scratch/nxu/astrocytes/from_collaborator/template_astro.yml"),
        fixed_collaborator_gtf,
        ribotie_cpm1_3sample,
        minlen_h5_ch,
        gencode_h5_ch
    )
    add_lncRNA_to_collaborator_csv(fix_RiboTIE_ORF_type.out.filtered, final_classification, annotation_gtf)

    emit:
    ribotie_res_merged_fixed_with_lncRNA = add_lncRNA_to_collaborator_csv.out.ribotie_csv_with_lncRNA
}