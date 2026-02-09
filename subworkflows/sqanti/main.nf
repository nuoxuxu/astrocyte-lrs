process extract_transcriptome {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "${params.outdir}/prepare/extract_transcriptome"

    input:
    path filtered_gtf

    output:
    path "joint_transcriptome.fasta"
    script:
    """
    gffread -w joint_transcriptome.fasta \\
    -g  ${params.ref_genome_fasta} \\
    ${filtered_gtf}
    """
}

process star_genomeGenerate {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/"

    input:
    path ref_genome_fasta
    path annotation_gtf
    val outputDir

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir $outputDir \\
        --genomeFastaFiles $ref_genome_fasta \\
        --sjdbGTFfile $annotation_gtf \\
        --sjdbOverhang ReadLength-1
    """

    output:
    path("${outputDir}")
}

process star_sr_genome {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/sr_genome"
    input:
    path star_genomeDir
    tuple val(sample_id), path(fastq_files)
    path annotation_gtf

    output:
    path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: star_aligned_bam
    path("${sample_id}.SJ.out.tab"), emit: star_sj_tab

    script:
    """
    STAR --runThreadN ${task.cpus} \\
    --genomeDir $star_genomeDir \\
    --readFilesIn $fastq_files \\
    --readFilesCommand gunzip -c \\
    --outFileNamePrefix "${sample_id}." \\
    --outSAMtype BAM SortedByCoordinate \\
    --sjdbGTFfile $annotation_gtf
    """
}

process sqanti_qc {
    label "mid_slurm_job"
    container "sqanti3_latest.sif"
    input:
    path gff_file
    path annotation_gtf
    path ref_genome_fasta
    path refTSS
    path polyA_motif_list
    path star_aligned_bam
    path star_sj_tab

    output:
    path("sqanti3_qc")

    script:
    """
    export PATH=/conda/miniconda3/envs/sqanti3/bin:\$PATH
    find . -name "*.bam" > SR_bam.fofn

    sqanti3_qc.py \\
    --isoforms $gff_file \\
    --refGTF $annotation_gtf \\
    --refFasta $ref_genome_fasta \\
    --CAGE_peak $refTSS \\
    --polyA_motif_list $polyA_motif_list \\
    --report html \\
    --skipORF \\
    --SR_bam SR_bam.fofn \\
    -c "*.SJ.out.tab" \\
    -o sqanti_qc_results \\
    -d sqanti3_qc \\
    -n 10
    """
}

process sqanti_filter {
    conda "/scratch/nxu/SQANTI3/env"
    label "short_slurm_job"
    
    input:
    path corrected_gtf
    path classification
    
    output:
    path("sqanti3_filter")
    
    script:
    """
    sqanti3_filter.py rules \\
    --filter_gtf $corrected_gtf \\
    --sqanti_class $classification \\
    -d sqanti3_filter \\
    -o default
    """
}

process filter_by_expression {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/sqanti3/isoseq/sqanti3_filter/${param_set_name}"

    input:
    path oarfish_quant_files
    path filtered_classification
    path filtered_gtf
    val param_set_name
    val min_reads
    val min_n_sample

    output:
    tuple val(param_set_name), path("final_classification.parquet"), emit: final_classification
    tuple val(param_set_name), path("final_transcripts.gtf"), emit: final_transcripts_gtf
    tuple val(param_set_name), path("final_expression.parquet"), emit: final_expression
    
    script:
    """
    get_counts_from_oarfish.py \\
        --min_reads ${min_reads} \\
        --min_n_sample ${min_n_sample}
    """
}

workflow SQANTI_AND_FILTER_BY_EXP {
    take:
    input_rna_fastq
    annotation_gtf
    ref_genome_fasta
    refTSS
    polyA_motif_list
    oarfish_quant
    merged_sorted_collapsed_gtf
    star_genomeGenerate_outputDir

    main:
    channel.fromFilePairs(input_rna_fastq).set { short_read_fastqs }
    star_genomeGenerate(ref_genome_fasta, annotation_gtf, star_genomeGenerate_outputDir)
    star_sr_genome(star_genomeGenerate.out, short_read_fastqs, annotation_gtf)
    sqanti_qc(merged_sorted_collapsed_gtf, annotation_gtf, ref_genome_fasta, refTSS, polyA_motif_list, star_sr_genome.out.star_aligned_bam.collect(), star_sr_genome.out.star_sj_tab.collect())    
    def isoseq_corrected_gtf = sqanti_qc.out
        .map { dir -> dir / "sqanti_qc_results_corrected.gtf" }
    def isoseq_classification = sqanti_qc.out
        .map { dir -> dir / "sqanti_qc_results_classification.txt" }
    sqanti_filter(isoseq_corrected_gtf, isoseq_classification)
    def filtered_gtf = sqanti_filter.out
        .map {dir -> dir / "default.filtered.gtf"}
    def filtered_classification = sqanti_filter.out
        .map {dir -> dir / "default_RulesFilter_result_classification.txt"}
    isoform_exp_filter_params = channel.fromList(params.filter_configs)        
    filter_by_expression_input_ch = oarfish_quant.collect().map { [it] }
        .combine(filtered_classification)
        .combine(filtered_gtf)
        .collect()
        .combine(isoform_exp_filter_params)
    filter_by_expression(filter_by_expression_input_ch)

    emit:
    final_classification = filter_by_expression.out.final_classification
    final_transcripts_gtf = filter_by_expression.out.final_transcripts_gtf
    final_expression = filter_by_expression.out.final_expression
    star_genomeDir = star_genomeGenerate.out
}