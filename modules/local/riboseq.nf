process star_riboseq {
    module "StdEnv/2023:star/2.7.11b"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/riboseq/"
    input:
    path star_genomeDir
    path riboseq_unmapped_to_contaminants
    path sjdbGTFfile
    val outputDir
    output:
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Aligned.sortedByCoord.out.bam"), emit: genome_bam
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.SJ.out.tab"), emit: sj_tab
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Log.final.out"), emit: log_final_out
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Log.out"), emit: log_out
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Signal.Unique.str1.out.bg"), emit: bedGraph_Unique
    path("${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}.Signal.UniqueMultiple.str1.out.bg"), emit: bedGraph_UniqueMultiple
    script:
    """
    STAR --runThreadN ${task.cpus} \\
    --genomeDir $star_genomeDir \\
    --readFilesIn ${riboseq_unmapped_to_contaminants} \\
    --outFileNamePrefix "${outputDir}/${riboseq_unmapped_to_contaminants.simpleName}." \\
    --outSAMtype BAM SortedByCoordinate \\
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
    """
}

process merge_bg_and_convert_to_bw {
    module "StdEnv/2023:bedtools/2.31.0:kent_tools/486"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/riboseq/"
    
    input:
    path bedgraph_files
    path chrom_sizes
    val prefix

    output:
    path "${prefix}_merged_riboseq.bw"
    
    script:
    """
    # Merge BedGraphs
    bedtools unionbedg -i ${bedgraph_files.join(' ')} -filler 0 | awk '{sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' OFS='\\t' > merged.bedgraph
    
    # Convert merged BedGraph to bigWig
    sort -k1,1 -k2,2n merged.bedgraph > merged_sorted.bedgraph
    bedGraphToBigWig merged_sorted.bedgraph ${chrom_sizes} ${prefix}_merged_riboseq.bw
    """
}