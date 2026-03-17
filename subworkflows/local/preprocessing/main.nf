process skera {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/skera"
    input:
    tuple val(sample_id), path(hifi_bam)
    path kinnex_adapters

    output:
    tuple val(sample_id), path("${sample_id}.skera.bam")

    script:
    """
    skera split \\
    ${hifi_bam} \\
    ${kinnex_adapters} \\
    ${sample_id}.skera.bam
    """
}

process lima {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/lima"
    input:
    tuple val(sample_id), path(skera_bam)
    path primer_fasta
    path biosamples_csv

    output:
    tuple val(sample_id), path("${sample_id}.fl.*.bam"), emit: demultiplexed_bam

    script:
    """
    lima \\
    ${skera_bam} \\
    ${primer_fasta} \\
    ${sample_id}.fl.bam \\
    --isoseq --peek-guess --overwrite-biosample-names
    """
}

process refine {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/refine"
    input:
    tuple val(sample_id), path(lima_bam)
    path primer_fasta

    output:
    tuple val(sample_id), path("${sample_id}.*.flnc.bam"), emit: flnc_bam
    path("${sample_id}.*.flnc.report.csv"), emit: flnc_report_csv

    script:
    """
    for input_file in ${lima_bam.join(' ')}; do
        temp_name="\${input_file#*fl.}"
        middle_part="\${temp_name%.bam}"
        output_file="\${middle_part}"
        isoseq3 refine \\
            \${input_file} \\
            ${primer_fasta} \\
            ${sample_id}.\${output_file}.flnc.bam
        done
    """
}

workflow PREPROCESSING {
    take:
    hifi_reads_bam
    kinnex_adapters
    isoseq_primers
    biosamples_csv

    main:
    channel.fromPath(params.hifi_reads_bam)
        .map { file -> 
            def sample_id = file.simpleName
            return tuple(sample_id, file) 
        }
        .set { hifi_reads_bam }
        
    skera(hifi_reads_bam, kinnex_adapters)
    lima(skera.out, isoseq_primers, biosamples_csv)
    refine(lima.out.demultiplexed_bam, isoseq_primers)

    emit:
    flnc_bam = refine.out.flnc_bam
}