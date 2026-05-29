process novelCodingJunctions {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocytes/env"
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
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/aim_2"

    input:
    path predicted_cds_gtf
    path novel_junctions
    path sqtl_file

    output:
    path "novel_coding_junction_sqtl_matches.tsv", emit: sqtl_matches

    script:
    """
    match_novel_junctions_sqtl.py \\
        $predicted_cds_gtf \\
        $novel_junctions \\
        $sqtl_file \\
        -o novel_coding_junction_sqtl_matches.tsv
    """
}

process matchNovelJunctionsSqtlColoc {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/aim_2"

    input:
    path sqtl_matches
    path coloc_file
    path iso_file

    output:
    path "novel_coding_junction_sqtl_coloc_matches.tsv", emit: coloc_matches

    script:
    """
    match_novel_junctions_sqtl_coloc.py \\
        $sqtl_matches \\
        $coloc_file \\
        $iso_file \\
        -o novel_coding_junction_sqtl_coloc_matches.tsv
    """
}

workflow AIM_2 {
    take:
    predicted_cds_gtf
    gencode_gtf
    sqtl_file
    coloc_file
    iso_file

    main:
    novelCodingJunctions(predicted_cds_gtf, gencode_gtf)
    matchNovelJunctionsSqtl(predicted_cds_gtf, novelCodingJunctions.out.novel_junctions, sqtl_file)
    matchNovelJunctionsSqtlColoc(matchNovelJunctionsSqtl.out.sqtl_matches, coloc_file, iso_file)

    emit:
    novel_junctions = novelCodingJunctions.out.novel_junctions
    sqtl_matches    = matchNovelJunctionsSqtl.out.sqtl_matches
    coloc_matches   = matchNovelJunctionsSqtlColoc.out.coloc_matches
}
