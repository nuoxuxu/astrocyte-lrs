process novelCodingJunctions {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/aim_2"

    input:
    path predicted_cds_gtf
    path gencode_gtf

    output:
    path "novel_coding_junctions.tsv", emit: novel_junctions

    script:
    """
    novel_coding_junctions.py \\
        $predicted_cds_gtf \\
        $gencode_gtf \\
        -o novel_coding_junctions.tsv
    """
}

process matchNovelJunctionsSqtl {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/aim_2"

    input:
    path predicted_cds_gtf
    path novel_junctions
    path sqtl_file
    path iso_file

    output:
    path "novel_coding_junction_sqtl_matches.tsv", emit: sqtl_matches

    script:
    """
    match_novel_junctions_sqtl.py \\
        $predicted_cds_gtf \\
        $novel_junctions \\
        $sqtl_file \\
        $iso_file \\
        -o novel_coding_junction_sqtl_matches.tsv
    """
}

process matchLeafcutterSqtl {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/aim_2"

    input:
    path sig_file
    path clu2gene
    path sqtl_file

    output:
    path "leafcutter_sqtl_matches.tsv", emit: sqtl_matches

    script:
    """
    match_leafcutter_sqtl.py \\
        $sig_file \\
        $clu2gene \\
        $sqtl_file \\
        -o leafcutter_sqtl_matches.tsv
    """
}

// Receives a [label, sqtl_matches] tuple; label is used as the output filename prefix.
process matchSqtlColoc {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/aim_2"

    input:
    tuple val(label), path(sqtl_matches)
    path coloc_file

    output:
    tuple val(label), path("${label}_sqtl_coloc_matches.tsv"), emit: coloc_matches

    script:
    """
    match_novel_junctions_sqtl_coloc.py \\
        $sqtl_matches \\
        $coloc_file \\
        -o ${label}_sqtl_coloc_matches.tsv
    """
}

// Receives a [label, coloc_matches] tuple; label is used as the output filename prefix.
process plotColocPie {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocyte-lrs/env"
    storeDir "nextflow_results/aim_2/figures"

    input:
    tuple val(label), path(coloc_matches)
    path bigbrain_coloc

    output:
    path "${label}_coloc_pie.pdf", emit: coloc_pie

    script:
    """
    plot_novel_junction_coloc_pie.py \\
        $coloc_matches \\
        $bigbrain_coloc \\
        -o ${label}_coloc_pie.pdf
    """
}

workflow AIM_2 {
    take:
    predicted_cds_gtf
    gencode_gtf
    sqtl_file
    coloc_file
    iso_file
    leafcutter_sig
    leafcutter_clu2gene

    main:
    novelCodingJunctions(predicted_cds_gtf, gencode_gtf)
    matchNovelJunctionsSqtl(predicted_cds_gtf, novelCodingJunctions.out.novel_junctions, sqtl_file, iso_file)
    matchLeafcutterSqtl(leafcutter_sig, leafcutter_clu2gene, sqtl_file)

    // Tag each sqtl_matches output with its label, then process both through the shared coloc chain
    sqtl_for_coloc = matchNovelJunctionsSqtl.out.sqtl_matches.map { f -> tuple("novel_coding_junction", f) }
        .mix(matchLeafcutterSqtl.out.sqtl_matches.map { f -> tuple("leafcutter", f) })

    matchSqtlColoc(sqtl_for_coloc, coloc_file)
    plotColocPie(matchSqtlColoc.out.coloc_matches, coloc_file)

    // Split the labelled coloc-match tuples back into per-source channels
    leafcutter_coloc = matchSqtlColoc.out.coloc_matches
        .filter { it[0] == "leafcutter" }.map { it[1] }
    novel_coding_junction_coloc = matchSqtlColoc.out.coloc_matches
        .filter { it[0] == "novel_coding_junction" }.map { it[1] }

    emit:
    novel_junctions         = novelCodingJunctions.out.novel_junctions
    novel_sqtl_matches      = matchNovelJunctionsSqtl.out.sqtl_matches
    leafcutter_sqtl_matches = matchLeafcutterSqtl.out.sqtl_matches
    leafcutter_coloc        = leafcutter_coloc
    novel_coding_junction_coloc = novel_coding_junction_coloc
}
