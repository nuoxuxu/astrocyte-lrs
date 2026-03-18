process ribotish_quality {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotish/quality/${param_set_name}"

    input:
    tuple val(param_set_name), path(genome_bam), path(genome_bam_index), path(orfanage_gtf)

    output:
    tuple val(param_set_name), path("${genome_bam.baseName}_qual.pdf"), emit: pdf
    tuple val(param_set_name), path("${genome_bam.baseName}_qual.txt"), emit: txt
    tuple val(param_set_name), path("${genome_bam}.para.py"), emit: para_py

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotish quality -b $genome_bam -g $orfanage_gtf
    """    
}

process ribotish_predict {
    module "python:gcc:arrow/19.0.1:rust"
    label "mid_slurm_job"
    storeDir "nextflow_results/ribotish/predict/${param_set_name}"
    input:
    tuple val(param_set_name), path(genome_bam), path(genome_bam_index), path(para_py), path(final_transcripts_gtf), path(ref_genome_fasta)
    
    output:
    tuple val(param_set_name), path("pred.txt"), emit: pred_txt

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    ribotish predict -b ${genome_bam.join(',')} \\
        -g $final_transcripts_gtf -f $ref_genome_fasta \\
        --longest -o pred.txt --ribopara ${para_py.join(',')} \\
        --fpth 0.8 --minpth 0.8 --fspth 0.8 --fsqth 1
    """
}

process ORFannotate {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotish/ORFannotate/${param_set_name}"

    input:
    tuple val(param_set_name), path(pred_txt_all), path(final_transcripts_gtf), path(ref_genome_fasta)

    output:
    tuple val(param_set_name), path("transcripts.fa"), emit: transcripts_fa
    tuple val(param_set_name), path("cds.fa"), emit: cds_fa
    tuple val(param_set_name), path("protein.fa"), emit: proteins_fa
    tuple val(param_set_name), path("ORFannotate_annotated.gtf"), emit: ORFannotate_annotated_gtf
    tuple val(param_set_name), path("ORFannotate_all_orfs.gtf"), emit: ORFannotate_all_orfs
    tuple val(param_set_name), path("ORFannotate_summary.tsv"), emit: ORFannotate_summary_tsv
    tuple val(param_set_name), path("ORFannotate.log"), emit: ORFannotate_log

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    python /home/nxu/tools/ORFannotate/ORFannotate.py --gtf $final_transcripts_gtf --fa $ref_genome_fasta --ribotish $pred_txt_all --outdir ./ --pvalue-cutoff 1
    """
}

process trim_stop_codons {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "nextflow_results/ribotish/ORFannotate/${param_set_name}"

    input:
    tuple val(param_set_name), path(ORFannotate_all_orfs)
    
    output:
    tuple val(param_set_name), path("ribotish.gtf")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    trim_ribotish_stop_codon.py --input_gtf $ORFannotate_all_orfs --output_gtf ribotish.gtf
    """
}

process save_to_parquet {
    module "python:gcc:arrow/19.0.1:rust"
    label "short_slurm_job"
    storeDir "export/local"
    
    input:
    tuple val(param_set_name), path(input_gtf)

    output:
    path("${input_gtf.baseName}.parquet")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    convert_gtf_to_parquet.py --input_gtf $input_gtf --output_gtf "${input_gtf.baseName}.parquet"
    """
}

workflow RUN_RIBOTISH {
    take:
    genome_bam
    final_transcripts_gtf
    orfanage_gtf
    ref_genome_fasta

    main:
    genome_bam
        .filter { param_set_name, _genome_bam, _genome_bam_index ->
            param_set_name == "mid_stringency"
        }
        .combine(orfanage_gtf, by: 0)
        | ribotish_quality

    genome_bam
        .groupTuple()
        .combine(ribotish_quality.out.para_py.groupTuple(), by: 0)\
        .combine(final_transcripts_gtf, by: 0)
        .combine(ref_genome_fasta)
        | ribotish_predict

    ribotish_predict.out
        .pred_txt
        .join(final_transcripts_gtf)
        .combine(ref_genome_fasta)
        | ORFannotate

    trim_stop_codons(ORFannotate.out.ORFannotate_all_orfs)

    trim_stop_codons.out
        | save_to_parquet

    emit:
    trimmed_gtf = trim_stop_codons.out
}