include { SUBREAD_FEATURECOUNTS } from "../../../modules/nf-core/subread/featurecounts"

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
    storeDir "nextflow_results/align/short_read/gencode"
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

workflow SHORT_READ {
    take:
    short_read_fastqs
    annotation_gtf
    ref_genome_fasta
    star_genomeGenerate_outputDir

    main:
    channel.fromFilePairs(short_read_fastqs).set { short_read_fastqs }

    star_genomeGenerate(ref_genome_fasta, annotation_gtf, star_genomeGenerate_outputDir)
    star_sr_genome(star_genomeGenerate.out, short_read_fastqs, annotation_gtf)

    star_sr_genome.out.star_aligned_bam
        .map { bam ->
            def meta = [id: bam.baseName.replaceAll(/\.Aligned.*/, ''), single_end: false]
            [meta, bam]
        }
        .combine(annotation_gtf)
    | SUBREAD_FEATURECOUNTS

    emit:
    star_genomeDir = star_genomeGenerate.out
    star_aligned_bam = star_sr_genome.out.star_aligned_bam
    star_sj_tab = star_sr_genome.out.star_sj_tab
}
