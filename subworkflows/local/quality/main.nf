process phylocsfpp {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(condition), path(ribotie_output_gtf), path(phyloCSF_db)

    output:
    tuple val(condition), path("${ribotie_output_gtf.baseName}.PhyloCSF++.gtf")

    script:
    """
    phylocsf++ annotate-with-tracks $phyloCSF_db/PhyloCSF+1.bw $ribotie_output_gtf
    """
}

process label_orf_type_gencode {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(ribotie_csv), path(annotation_gtf), path(orfanage_gtf), path(final_classification)

    output:
    tuple val(name), path("${name}_orf_type_gencode.tsv"), emit: orf_type_gencode

    script:
    """
    export POLARS_MAX_THREADS=1
    label_orf_type_gencode.py \\
        $ribotie_csv \\
        $annotation_gtf \\
        $orfanage_gtf \\
        $final_classification \\
        -o ${name}_orf_type_gencode.tsv
    """
}

process convert_gtf_to_orbl {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/ORBL"

    input:
    tuple val(name), path(ribotie_gtf)

    output:
    tuple val(name), path("${name}_orbl_input.tsv"), emit: orbl_input

    script:
    """
    gtf_cds_to_orbl.py $ribotie_gtf \\
        --group-by ORF_id \\
        --extra-fields ORF_id transcript_id ribotie_score \\
        > ${name}_orbl_input.tsv
    """
}

process split_orbl_input {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/ORBL/input_chunks"
    input:
    tuple val(name), path(orbl_input)

    output:
    tuple val(name), path("${name}_chunk_*.tsv"), emit: chunks

    script:
    """
    header=\$(head -1 $orbl_input)
    total=\$(tail -n +2 $orbl_input | wc -l)
    chunk_size=\$(( (total + 9) / 10 ))
    tail -n +2 $orbl_input | split -l \$chunk_size - split_
    i=0
    for f in split_*; do
        printf '%s\\n' "\$header" | cat - \$f > ${name}_chunk_\$(printf '%02d' \$i).tsv
        i=\$((i + 1))
    done
    """
}

process run_orbl {
    // No label — runs on login node where internet access is available
    storeDir "nextflow_results/quality/ORBL/output_chunks"

    input:
    tuple val(name), path(orbl_input), val(alignment_set)

    output:
    tuple val(name), path("${orbl_input.baseName}_orbl_output.tsv"), emit: orbl_output

    script:
    """
    python /scratch/nxu/ORBL_tools/orbl.py \\
        $alignment_set \\
        $orbl_input \\
        ${orbl_input.baseName}_orbl_output.tsv \\
        --header \\
        --extraFields ORF_id,transcript_id,ribotie_score
    """
}

process merge_orbl_chunks {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/ORBL"

    input:
    tuple val(name), path(orbl_chunks)

    output:
    tuple val(name), path("${name}_orbl_output.tsv"), emit: orbl_output

    script:
    """
    first=1
    for f in \$(ls *.tsv | sort); do
        if [ "\$first" -eq 1 ]; then
            cat "\$f" > ${name}_orbl_output.tsv
            first=0
        else
            tail -n +2 "\$f" >> ${name}_orbl_output.tsv
        fi
    done
    """
}

process add_biotype_with_frameshift {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(quality_metrics), path(ribotie_csv),
          path(orfanage_gtf), path(final_classification), path(annotation_gtf)

    output:
    tuple val(name), path("${name}_quality_metrics_biotype.tsv"), emit: quality_metrics

    script:
    """
    export POLARS_MAX_THREADS=1
    add_biotype_with_frameshift.py \\
        $quality_metrics \\
        $ribotie_csv \\
        $annotation_gtf \\
        $orfanage_gtf \\
        $final_classification \\
        -o ${name}_quality_metrics_biotype.tsv
    """
}

process merge_orbl_quality {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(orf_type_gencode), path(orbl_output),
          path(ribotie_csv), path(final_classification), path(gencode_ribotie_csv, stageAs: 'gencode_ribotie_res_merged.csv')

    output:
    tuple val(name), path("${name}_quality_metrics.tsv"), emit: quality_metrics

    script:
    """
    merge_orbl_quality.py $orf_type_gencode $orbl_output \
        --ribotie-csv $ribotie_csv \
        --gencode-ribotie $gencode_ribotie_csv \
        --classification $final_classification \
        -o ${name}_quality_metrics.tsv
    """
}

workflow GET_QUALITY_METRICS {
    take:
    ribotie_training_outputs
    PhyloCSFpp_db

    main:
    channel.fromPath(ribotie_training_outputs)
        .splitJson()
        .map { entry ->
            tuple("merged", file(entry.ribotie_merged_novel_gtf))
        }
        .combine(PhyloCSFpp_db)
        | phylocsfpp
}

workflow LABEL_ORF_TYPE_GENCODE {
    take:
    main_pipeline_outputs    // glob path: nextflow_results/manifests/main_pipeline_outputs_*.json
    ribotie_training_outputs // glob path: nextflow_results/manifests/ribotie_training_outputs_*.json
    annotation_gtf
    orbl_alignment_set       // val: e.g. 'hg38_100'

    main:
    // Build channel: (name, orfanage_gtf, final_classification) from main_pipeline_outputs
    channel.fromPath(main_pipeline_outputs)
        .map { f ->
            def name = f.baseName.replaceAll(/^main_pipeline_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.orfanage_gtf), file(entry.final_classification))
        }
        .set { main_outputs_ch }

    // Build combined channel from ribotie_training_outputs (inner join filters out gencode)
    // Tuple: (name, ribotie_csv, ribotie_gtf, orfanage_gtf, final_classification)
    channel.fromPath(ribotie_training_outputs)
        .map { f ->
            def name = f.baseName.replaceAll(/^ribotie_training_outputs_/, '')
            def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
            tuple(name, file(entry.ribotie_merged_csv), file(entry.ribotie_merged_gtf))
        }
        .join(main_outputs_ch)
        .set { combined_ch }

    // GENCODE ORF type labeling
    combined_ch
        .combine(annotation_gtf)
        .map { name, csv, gtf, orfanage, clf, anno -> tuple(name, csv, anno, orfanage, clf) }
        | label_orf_type_gencode

    // ORBL scoring chain
    combined_ch
        .map { name, csv, gtf, orfanage, clf -> tuple(name, gtf) }
        | convert_gtf_to_orbl

    convert_gtf_to_orbl.out.orbl_input
        | split_orbl_input

    split_orbl_input.out.chunks
        .transpose()
        .combine(orbl_alignment_set)
        | run_orbl

    run_orbl.out.orbl_output
        .groupTuple()
        | merge_orbl_chunks

    // Extract GENCODE ribotie CSV for GENCODE label matching
    channel.fromPath(ribotie_training_outputs)
        .map { f ->
            def name = f.baseName.replaceAll(/^ribotie_training_outputs_/, '')
            if (name == 'gencode') {
                def entry = new groovy.json.JsonSlurper().parseText(f.text)[0]
                tuple(file(entry.ribotie_merged_csv))
            }
        }
        .first()
        .set { gencode_ribotie_csv }

    // Merge: join GENCODE label output, ORBL output, ribotie data, and classification by name
    label_orf_type_gencode.out.orf_type_gencode
        .join(merge_orbl_chunks.out.orbl_output)
        .join(combined_ch.map { name, csv, gtf, orfanage, clf ->
            tuple(name, csv, clf)
        })
        .combine(gencode_ribotie_csv)
        .map { name, orf_type, orbl_out, ribotie_c, clf, gencode_rib ->
            tuple(name, orf_type, orbl_out, ribotie_c, clf, gencode_rib)
        }
        | merge_orbl_quality

    // Add biotype_with_frameshift column using GENCODE-projected frame
    merge_orbl_quality.out.quality_metrics
        .join(combined_ch.map { name, csv, gtf, orfanage, clf ->
            tuple(name, csv, orfanage, clf)
        })
        .combine(annotation_gtf)
        | add_biotype_with_frameshift

    emit:
    quality_metrics  = add_biotype_with_frameshift.out.quality_metrics
    orf_type_gencode = label_orf_type_gencode.out.orf_type_gencode
}
