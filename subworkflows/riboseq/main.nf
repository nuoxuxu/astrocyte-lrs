process format_gtf_for_ribotie {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage/${param_set_name}"
    
    input:
    tuple val(param_set_name), path(orfanage_gtf)
    tuple val(param_set_name_2), path(final_classification)
    path annotation_gtf

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

// Subset GENCODE annotation to only include transcripts that are present in the custom GTF
process subset_GENCODE_tx {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/"

    input:
    tuple val(param_set_name), path(final_classification), path(tmap), path(annotation_gtf)

    output:
    tuple val("gencode_low"), path("gencode_supported_transcripts.gtf")
    
    script:
    """
    subset_gencode_by_tmap.py \\
        --tmap $tmap \\
        --gtf $annotation_gtf \\
        --output gencode_supported_transcripts.gtf
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

process generate_ribotie_yml {
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}"

    input:
    tuple val(param_set_name), path(gtf_path)
    path ref_genome_fasta
    tuple val(param_set_name_2), path(transcriptome_bam)

    output:
    tuple val(param_set_name), path("RiboTIE.yml")

    script:
    """
    generate_ribotie_yml.py \
    --gtf $gtf_path \
    --fa $ref_genome_fasta \
    --bam-glob "*.bam" \
    --h5 ribotie_res.h5 \
    --out-prefix ribotie_res \
    --samples "Unstim:merged_astro_A,merged_astro_B" "Stim:merged_astro_C,Astro_D" \
    -o RiboTIE.yml
    """
}

process generate_ribotie_db {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotie/${param_set_name}"

    input:
    tuple val(param_set_name), path(gtf_path), path(transcriptome_bam), path(ribotie_yml), path(ref_genome_fasta)
    
    output:
    tuple val(param_set_name), path("${gtf_path.baseName}.h5"), path("ribotie_res.h5")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotie $ribotie_yml --data
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
    tmap

    main:
    final_classification
        .filter { label, file -> 
            label == "low_stringency" 
        }
        .join(tmap)
        .combine(annotation_gtf)
        | subset_GENCODE_tx

    channel.fromPath(riboseq_unmapped_to_contaminants).set { riboseq_unmapped_to_contaminants }
    
    custom_gtf = format_gtf_for_ribotie(orfanage_gtf, final_classification, annotation_gtf)
    
    sjdbGTFfile_tuples = custom_gtf
        .mix(
            channel.value("gencode")
                .concat(annotation_gtf)
                .collect()
        )
        .mix(
            subset_GENCODE_tx.out
        )

    star_genomeDir
        .combine(riboseq_unmapped_to_contaminants)
        .combine(sjdbGTFfile_tuples)        
        | star_riboseq

    generate_ribotie_yml(sjdbGTFfile_tuples, ref_genome_fasta, star_riboseq.out.transcriptome_bam.groupTuple())
    
    sjdbGTFfile_tuples
        .join(star_riboseq.out.transcriptome_bam.groupTuple())
        .join(generate_ribotie_yml.out)
        .combine(ref_genome_fasta)
        | generate_ribotie_db

    generate_ribotie_db.out
        .toList()
        .map { results ->
            def outputs = results.collect { param_name, gtf_h5, ribotie_h5 ->
                [
                    param_set_name: param_name,
                    gtf_h5: "nextflow_results/ribotie/${param_name}/${gtf_h5.name}",
                    ribotie_h5: "nextflow_results/ribotie/${param_name}/${ribotie_h5.name}",
                    base_dir: "nextflow_results/ribotie/${param_name}",
                    ribotie_yml: "nextflow_results/ribotie/${param_name}/RiboTIE.yml"
                ]
            }
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputs))
        }
        .collectFile(name: 'ribotie_training_inputs.json', storeDir: 'nextflow_results/manifests')
}