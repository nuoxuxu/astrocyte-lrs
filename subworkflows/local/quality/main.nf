process phylocsfpp {
    conda "/scratch/nxu/astrocyte-lrs/env"
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

process convert_gtf_to_orbl {
    conda "/scratch/nxu/astrocyte-lrs/env"
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
    conda "/scratch/nxu/astrocyte-lrs/env/"
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
    conda "/scratch/nxu/astrocyte-lrs/env"
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

process annotate_is_novel {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(quality_metrics), path(ribotie_gtf), path(final_classification), path(annotation_gtf)

    output:
    tuple val(name), path("${name}_quality_metrics_novel.tsv"), emit: quality_metrics

    script:
    """
    export POLARS_MAX_THREADS=1
    annotate_is_novel.py \\
        $quality_metrics \\
        $ribotie_gtf \\
        $final_classification \\
        $annotation_gtf \\
        -o ${name}_quality_metrics_novel.tsv
    """
}

process merge_orbl_quality {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(ribotie_csv), path(orbl_output)

    output:
    tuple val(name), path("${name}_quality_metrics.tsv"), emit: quality_metrics

    script:
    """
    merge_orbl_quality.py $ribotie_csv $orbl_output -o ${name}_quality_metrics.tsv
    """
}

process label_orf_type_gencode {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(ribotie_csv), path(orfanage_numbered_exons_gtf),
          path(final_classification), path(orfanage_plain_gtf), path(annotation_gtf)

    output:
    tuple val(name), path("${name}_orf_type_gencode.tsv")

    script:
    """
    export POLARS_MAX_THREADS=1
    label_orf_type_gencode.py \\
        $ribotie_csv \\
        $annotation_gtf \\
        $orfanage_numbered_exons_gtf \\
        $final_classification \\
        $orfanage_plain_gtf \\
        -o ${name}_orf_type_gencode.tsv
    """
}

process add_riboseq_coverage {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(quality_metrics)

    output:
    tuple val(name), path("${name}_quality_metrics_riboseq.tsv"), emit: quality_metrics

    script:
    """
    add_riboseq_coverage.py \\
        $quality_metrics \\
        ${params.riboseq_bw} \\
        --bigwig-tool ${params.bigwig_average_over_bed} \\
        -o ${name}_quality_metrics_riboseq.tsv
    """
}

process plot_orf_type_bars {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/figures"

    input:
    tuple val(name), path(orf_type_gencode), path(ribotie_gtf),
          path(orfanage_plain_gtf), path(gencode_gtf)

    output:
    tuple val(name), path("${name}_orf_type_bars_orfanage.png"),
                     path("${name}_orf_type_bars_ribotie.png")

    script:
    """
    plot_orf_type_bars.py $orf_type_gencode \\
        --name $name \\
        --prefix ${name}_orf_type_bars \\
        --ribotie-gtf $ribotie_gtf \\
        --orfanage-gtf $orfanage_plain_gtf \\
        --gencode-gtf $gencode_gtf
    """
}

process merge_orf_type_gencode {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality"

    input:
    tuple val(name), path(quality_metrics_riboseq), path(orf_type_gencode)

    output:
    tuple val(name), path("${name}_final_quality_metrics.tsv"), emit: final_quality_metrics

    script:
    """
    export POLARS_MAX_THREADS=1
    merge_orf_type_gencode.py \\
        $quality_metrics_riboseq \\
        $orf_type_gencode \\
        -o ${name}_final_quality_metrics.tsv
    """
}

process gtf_to_bigbed {
    conda "/scratch/nxu/astrocyte-lrs/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quality/UCSC"

    input:
    tuple val(track_label), path(gtf)

    output:
    path("${track_label}.bb")

    script:
    """
    gtf_to_bed12.py $gtf tmp.bed
    awk 'NR==FNR{valid[\$1]=1; next} \$1 in valid' ${params.chrom_sizes} tmp.bed \\
        | sort -k1,1 -k2,2n > sorted.bed
    /home/nxu/miniforge3/bin/bedToBigBed -type=bed12 sorted.bed ${params.chrom_sizes} ${track_label}.bb
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

    // Merge: join RiboTIE CSV and ORBL output by name
    combined_ch
        .map { name, csv, gtf, orfanage, clf -> tuple(name, csv) }
        .join(merge_orbl_chunks.out.orbl_output)
        | merge_orbl_quality

    // Annotate is_novel: 1 if any CDS segment overlaps a non-GENCODE-exon region
    merge_orbl_quality.out.quality_metrics
        .join(combined_ch.map { name, csv, gtf, orfanage, clf -> tuple(name, gtf, clf) })
        .combine(annotation_gtf)
        | annotate_is_novel

    // Add Ribo-seq coverage over novel regions
    annotate_is_novel.out.quality_metrics
        | add_riboseq_coverage

    // GTF → bigBed: RiboTIE and ORFanage (mix into one channel; gencode excluded by combined_ch inner join)
    combined_ch
        .map { name, csv, gtf, orfanage, clf -> tuple("${name}_ribotie", gtf) }
        .mix(combined_ch.map { name, csv, gtf, orfanage, clf -> tuple("${name}_orfanage", orfanage) })
        | gtf_to_bigbed

    // Label ORF types against GENCODE canonical CDS
    channel.of(
        tuple('no_minlen', file(params.orfanage_gtf_no_minlen)),
        tuple('minlen',    file(params.orfanage_gtf_minlen))
    ).set { orfanage_plain_ch }

    combined_ch
        .map { name, csv, gtf, orfanage, clf -> tuple(name, csv, orfanage, clf) }
        .join(orfanage_plain_ch)
        .combine(annotation_gtf)
        | label_orf_type_gencode

    // Bar chart: in-frame ORF type distribution per version, with GENCODE CDS bin overlay
    label_orf_type_gencode.out
        .join(combined_ch.map { name, csv, gtf, orfanage, clf -> tuple(name, gtf) })
        .join(orfanage_plain_ch)
        .combine(annotation_gtf)
        | plot_orf_type_bars

    // Merge ORF_type_ORFanage + ORF_type_RiboTIE columns into the riboseq quality metrics
    add_riboseq_coverage.out.quality_metrics
        .join(label_orf_type_gencode.out)
        | merge_orf_type_gencode

    emit:
    quality_metrics = add_riboseq_coverage.out.quality_metrics
    final_quality_metrics = merge_orf_type_gencode.out.final_quality_metrics
}
