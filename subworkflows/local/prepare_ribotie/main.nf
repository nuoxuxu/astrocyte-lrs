process format_gtf_for_ribotie {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orf_prediction/orfanage/${param_set_name}"
    
    input:
    tuple val(param_set_name), path(orfanage_gtf), path(final_classification), path(annotation_gtf)

    output:
    tuple val(param_set_name), path("orfanage_numbered_exons.gtf")

    script:
    """
    format_gtf_for_ribotie.py \\
    --input_gtf ${orfanage_gtf} \\
    --final_classification ${final_classification} \\
    --annotation_gtf ${annotation_gtf} \\
    --output_gtf orfanage_numbered_exons.gtf
    """
}

process rename_chromosomes {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/other"

    input:
    path(chrom_sizes)
    path(chromAlias)

    output:
    path("hg38.p13.GENCODE_chrom.size")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    make_gencode_chrom_sizes.py \\
        --chromAlias $chromAlias \\
        --chrom_sizes $chrom_sizes > hg38.p13.GENCODE_chrom.size
    """
}
process star_riboseq {
    module "StdEnv/2023:star/2.7.11b:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/align/riboseq/${param_set_name}"
    
    input:
    tuple path(star_genomeDir), path(riboseq_unmapped_to_contaminants), val(param_set_name), path(sjdbGTFfile)
    
    output:
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Aligned.sortedByCoord.out.bam"), path("${riboseq_unmapped_to_contaminants.simpleName}.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bam
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.SJ.out.tab"), emit: sj_tab
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Log.final.out"), emit: log_final_out
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Log.out"), emit: log_out
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Signal.Unique.str1.out.bg"), emit: bedGraph_Unique
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Signal.UniqueMultiple.str1.out.bg"), emit: bedGraph_UniqueMultiple
    
    script:
    """
    STAR --runThreadN ${task.cpus} \\
    --genomeDir $star_genomeDir \\
    --readFilesIn ${riboseq_unmapped_to_contaminants} \\
    --outFileNamePrefix "${riboseq_unmapped_to_contaminants.simpleName}." \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMattributes All \\
    --limitBAMsortRAM 31568141173 \\
    --sjdbGTFfile $sjdbGTFfile \\
    --quantMode TranscriptomeSAM \\
    --outFilterMultimapNmax 10 \\
    --outMultimapperOrder Random \\
    --outFilterMismatchNmax 2 \\
    --seedSearchStartLmaxOverLread 0.5 \\
    --alignEndsType EndToEnd \\
    --outWigType bedGraph \\
    --outWigStrand Unstranded

    samtools index -@ ${task.cpus} "${riboseq_unmapped_to_contaminants.simpleName}.Aligned.sortedByCoord.out.bam"
    """
}

//TO DO: Fix the --samples argument so that it takes in sample JSON file
process generate_ribotie_yml {
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/prepare_ribotie/${param_set_name}/${mode}"

    input:
    tuple val(param_set_name), path(gtf_path), path(transcriptome_bam), val(mode), path(ref_genome_fasta)
    output:
    tuple val(param_set_name), val(mode), path("RiboTIE.yml")

    script:
    """
    generate_ribotie_yml.py \\
    --gtf $gtf_path \\
    --fa $ref_genome_fasta \\
    --bam-glob "*.bam" \\
    --h5 ribotie_res.h5 \\
    --out-prefix ribotie_res \\
    --samples "Unstim:merged_astro_A,merged_astro_B" "Stim:merged_astro_C,Astro_D" \\
    --mode $mode \\
    -o RiboTIE.yml
    """
}

process generate_ribotie_db {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare_ribotie/${param_set_name}/${mode}"

    input:
    tuple val(param_set_name), path(gtf_path), path(transcriptome_bam), val(mode), path(ribotie_yml), path(ref_genome_fasta)
    
    output:
    tuple val(param_set_name), path("${gtf_path.baseName}.h5"), path("ribotie_res.h5"), val(mode)

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotie $ribotie_yml --data
    """
}

process merge_bg_and_convert_to_bw {
    module "StdEnv/2023:bedtools/2.31.0:kent_tools/486"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/riboseq/${param_set_name}"
    
    input:
    tuple val(param_set_name), path(bedgraph_files), path(chrom_sizes)

    output:
    path("merged_riboseq.bw")
    
    script:
    """
    # Merge BedGraphs
    bedtools unionbedg -i ${bedgraph_files.join(' ')} -filler 0 | awk '{sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' OFS='\\t' > merged.bedgraph
    
    # Convert merged BedGraph to bigWig
    sort -k1,1 -k2,2n merged.bedgraph > merged_sorted.bedgraph
    bedGraphToBigWig merged_sorted.bedgraph ${chrom_sizes} merged_riboseq.bw
    """
}

workflow PREPARE_RIBOTIE {
    take:
    orfanage_gtf
    final_classification
    annotation_gtf
    star_genomeDir
    riboseq_unmapped_to_contaminants
    ref_genome_fasta
    chrom_sizes
    chromAlias

    main:

    // Format ORFanage GTF for use with RiboTIE 
    orfanage_gtf
        .join(final_classification)
        .combine(annotation_gtf)
        | format_gtf_for_ribotie
    
    orfanage_numbered_exons_gtf = format_gtf_for_ribotie.out
        .filter { param_set_name, _orfanage_numbered_exons_gtf ->
            param_set_name == "mid_stringency"
        }

    // Align Ribo-seq reads to transcriptome
    star_genomeDir
        .combine(channel.fromPath(riboseq_unmapped_to_contaminants))
        .combine(orfanage_numbered_exons_gtf)
        | star_riboseq

    // Rename UCSC chromosome names to GENCODE chromosome names
    rename_chromosomes(chrom_sizes, chromAlias)

    // Merge all bedGraph files generated by STAR and convert to bigWig
    star_riboseq.out.bedGraph_Unique
        .groupTuple()
        .combine(rename_chromosomes.out)
        | merge_bg_and_convert_to_bw

    // Generate the yml file that contains RiboTIE input arguments
    // For now only run on the "merged" mode
    orfanage_numbered_exons_gtf
        .join(star_riboseq.out.transcriptome_bam.groupTuple())
        .combine(channel.of("merged"))
        .combine(ref_genome_fasta)
        | generate_ribotie_yml

    // Generate RiboTIE h5 databse
    orfanage_numbered_exons_gtf
        .join(star_riboseq.out.transcriptome_bam.groupTuple())
        .cross(generate_ribotie_yml.out)
        .combine(ref_genome_fasta)
        .map { v -> tuple(v[0][0], v[0][1], v[0][2], v[1][1], v[1][2], v[2]) }
        | generate_ribotie_db

    emit:
    ribotie_db = generate_ribotie_db.out
    genome_bam = star_riboseq.out.genome_bam
    transcriptome_bam = star_riboseq.out.transcriptome_bam
}