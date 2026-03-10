process make_combined_gtf_for_ribotie {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare_ribotie/${param_set_name}"

    input:
    tuple val(param_set_name), path(final_transcripts_gtf), path(final_classification), path(annotation_gtf)

    output:
    tuple val(param_set_name), path("combined_for_ribotie.gtf")

    script:
    """
    make_combined_gtf_for_ribotie.py \\
        --final_transcripts_gtf ${final_transcripts_gtf} \\
        --final_classification ${final_classification} \\
        --annotation_gtf ${annotation_gtf} \\
        --output_gtf combined_for_ribotie.gtf
    """
}


process sort_gtf {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare_ribotie/${param_set_name}"

    input:
    tuple val(param_set_name), path(input_gtf)

    output:
    tuple val(param_set_name), path("${input_gtf.baseName}.sorted.gtf")

    script:
    """
    grep "^#" ${input_gtf} > ${input_gtf.baseName}.sorted.gtf || true
    grep -v "^#" ${input_gtf} | sort -k1,1 -k4,4n -k5,5n >> ${input_gtf.baseName}.sorted.gtf
    """
}

process star_riboseq {
    module "StdEnv/2023:star/2.7.11b"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/riboseq/${param_set_name}"
    input:
    tuple path(star_genomeDir), path(riboseq_unmapped_to_contaminants), val(param_set_name), path(sjdbGTFfile)
    output:
    tuple val(param_set_name), path("${riboseq_unmapped_to_contaminants.simpleName}.Aligned.sortedByCoord.out.bam"), emit: genome_bam
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

process remove_unplaced_chromosomes_from_bg {
    module "StdEnv/2023:star/2.7.11b"
    label "short_slurm_job"
    storeDir "nextflow_results/align/star/riboseq/${param_set_name}"

    input:
    tuple val(param_set_name), path(bedGraph_Unique)

    output:
    tuple val(param_set_name), path("filtered_${bedGraph_Unique}")

    script:
    """
    awk '\$1 ~ /^chr[0-9XYM]+\$/ {print \$0}' $bedGraph_Unique > "filtered_${bedGraph_Unique}"
    """
}

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
    tuple val(param_set_name), path(bedgraph_files), path(chrom_sizes), val(prefix)

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

workflow PREPARE_RIBOTIE {
    take:
    final_transcripts_gtf
    final_classification
    annotation_gtf
    star_genomeDir
    riboseq_unmapped_to_contaminants
    ref_genome_fasta

    main:

    // Combine Isoseq novel trqanscripts (no CDS)
    // with GENCODE transcripts that correspond to known Isoseq transcripts
    // as training input for RiboTIE
    
    final_transcripts_gtf
        .join(final_classification)
        .combine(annotation_gtf)
        | make_combined_gtf_for_ribotie

    make_combined_gtf_for_ribotie.out
        | sort_gtf

    sjdbGTFfile_tuples = sort_gtf.out
        .mix(
            channel.value("gencode")
                .concat(annotation_gtf)
                .collect()
        )

    star_genomeDir
        .combine(channel.fromPath(riboseq_unmapped_to_contaminants))
        .combine(sjdbGTFfile_tuples)        
        | star_riboseq

    filtered_bg = star_riboseq.out
        .bedGraph_Unique
        | remove_unplaced_chromosomes_from_bg

    unstim_filtered_bg = filtered_bg
        .filter{ prefix, bedGraph_Unique -> (prefix == "mid_stringency") && (bedGraph_Unique.toString().contains("merged_astro_A") || bedGraph_Unique.toString().contains("merged_astro_B")) }
        .groupTuple()
        .combine(channel.of(file(params.chrom_sizes)))
        .combine(channel.of("unstim"))

    filtered_bg
        .filter{ prefix, bedGraph_Unique -> (prefix == "mid_stringency") && (bedGraph_Unique.toString().contains("merged_astro_A") || bedGraph_Unique.toString().contains("merged_astro_B")) }
        .groupTuple()
        .combine(channel.of(file(params.chrom_sizes)))
        .combine(channel.of("stim"))
        .mix(unstim_filtered_bg)
        | merge_bg_and_convert_to_bw

    sjdbGTFfile_tuples
        .join(star_riboseq.out.transcriptome_bam.groupTuple())
        .combine(channel.of("merged", "separate"))
        .combine(ref_genome_fasta)
        | generate_ribotie_yml

    sjdbGTFfile_tuples
        .join(star_riboseq.out.transcriptome_bam.groupTuple())
        .cross(generate_ribotie_yml.out)
        .combine(ref_genome_fasta)
        .map { v -> tuple(v[0][0], v[0][1], v[0][2], v[1][1], v[1][2], v[2]) }
        .filter { param_set_name, _gtf_path, _transcriptome_bam, mode, _ribotie_yml, _ref_genome_fasta ->
            param_set_name == "mid_stringency" && mode == "separate"
        }
        | generate_ribotie_db

    emit:
    ribotie_db = generate_ribotie_db.out
}