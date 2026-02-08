include { sqanti_qc as sqanti_qc_bambu } from "./modules/local/sqanti3"
include { sqanti_filter as sqanti_filter_bambu } from "./modules/local/sqanti3"
include { sqanti_qc as sqanti_qc_isoseq } from "./modules/local/sqanti3"
include { sqanti_filter as sqanti_filter_isoseq } from "./modules/local/sqanti3"
include { star_riboseq as star_riboseq_custom } from "./modules/local/riboseq"
include { star_riboseq as star_riboseq_gencode } from "./modules/local/riboseq"
include { merge_bg_and_convert_to_bw as merge_bg_and_convert_to_bw_custom_unstim } from "./modules/local/riboseq"
include { merge_bg_and_convert_to_bw as merge_bg_and_convert_to_bw_custom_stim } from "./modules/local/riboseq"

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

process merge_flnc_bams {
    conda "/scratch/nxu/SQANTI3/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/prepare/merge_flnc_bams"
    
    input:
    path(flnc_bam)

    output:
    path("*.flnc.bam")

    script:
    """
    merge_flnc_bam.sh
    """
}

process fofn {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/fofn"
    input:
    path(flnc_bam)

    output:
    path("*.fofn")

    script:
    """
    fofn.py
    """
}

process cluster {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/cluster"
    input:
    path fofn
    path flnc_bam

    output:
    path("${fofn.simpleName}.clustered.bam")

    script:
    """
    isoseq cluster2 \\
    $fofn \\
    ${fofn.simpleName}.clustered.bam
    """
}

process convert_flnc_bam_to_fastqz {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/convert_ubam_to_fastqz"

    input:
    path merged_flnc_bam
    output:
    path("${merged_flnc_bam.simpleName}.fastq.gz")

    script:
    """
    samtools fastq -@ ${task.cpus} -0 ${merged_flnc_bam.simpleName}.fastq.gz \\
    -c ${params.compression_level} $merged_flnc_bam
    """
}

process minimap2_genome {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/minimap2_genome"    
    input:
    path fastqz
    output:
    path("${fastqz.simpleName}.aligned.bam")
    script:
    """
    minimap2 -ax splice:hq -t ${task.cpus} -uf ${params.ref_genome_fasta} $fastqz | samtools sort -@ ${task.cpus} \\
        -o "${fastqz.simpleName}.aligned.bam"
    """
}

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
process minimap2_transcriptome {
    conda "/scratch/nxu/astrocytes/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/align/minimap2_transcriptome"

    input:
    path extracted_transcriptome
    path fastqz_reads

    output:
    path("${fastqz_reads.simpleName}.aligned.bam")

    script:
    """
    minimap2 --eqx -N 100 -ax map-hifi \\
        -t ${task.cpus} $extracted_transcriptome \\
        $fastqz_reads | \\
        samtools sort -n -@ ${task.cpus} \\
        -o "${fastqz_reads.simpleName}.aligned.bam"
    """
}

process run_oarfish {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/quantify/isoseq/oarfish"

    input:
    path transcriptome_alignment
    
    output:
    path("${transcriptome_alignment.simpleName}.quant")

    script:
    """
    oarfish --threads ${task.cpus} \\
        --filter-group no-filters \\
        --model-coverage \\
        --alignments $transcriptome_alignment \\
        --output "${transcriptome_alignment.simpleName}"
    """

}
process bambu {
    conda "/scratch/nxu/astrocytes/env"
    label "mid_slurm_job"
    storeDir "nextflow_results/discover/bambu"
    input:
    path minimap_bam
    output:
    path("supportedTranscriptModels.gtf"), emit: supported_tx_gtf
    path("novelTranscripts.gtf"), emit: novel_tx_gtf
    script:
    """
    run_bambu.R ${params.annotation_gtf} "*minimap.bam" ${params.ref_genome_fasta} ${task.cpus}
    """
}

process pbmm2 {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/align/pbmm2"
    input:
    path ref_genome_fasta
    path clustered_bam

    output:
    path("${clustered_bam.simpleName}.clustered.aligned.bam"), emit: aligned_bam
    path("${clustered_bam.simpleName}.clustered.aligned.bam.bai"), emit: aligned_bam_bai

    script:
    """
    pbmm2 align --preset ISOSEQ --sort \\
    $ref_genome_fasta \\
    $clustered_bam \\
    "${clustered_bam.baseName}.aligned.bam"
    """
}

process merge_aligned_bams {
    conda "/scratch/nxu/SQANTI3/env"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/merge_aligned_bams"
    input:
    path aligned_bams
    output:
    path("merged.bam")
    script:
    """
    samtools merge -o merged.bam *.aligned.bam
    """
}

process collapse_and_sort {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/discover/isoseq"
    input:
    path merged_bam

    output:
    path("${merged_bam.baseName}.sorted.collapsed.gff")

    script:
    """
    isoseq collapse \\
    --do-not-collapse-extra-5exons \\
    $merged_bam \\
    "${merged_bam.baseName}.collapsed.gff"

    pigeon sort -o "${merged_bam.baseName}.sorted.collapsed.gff" "${merged_bam.baseName}.collapsed.gff"
    """
}

process filter_by_expression {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/sqanti3/isoseq/sqanti3_filter/"

    input:
    path oarfish_quant_files
    path filtered_classification
    path filtered_gtf
    val min_reads
    val min_n_sample

    output:
    path("final_expression_${min_reads}_${min_n_sample}/final_classification.parquet"), emit: final_classification
    path("final_expression_${min_reads}_${min_n_sample}/final_transcripts.gtf"), emit: final_transcripts_gtf
    path("final_expression_${min_reads}_${min_n_sample}/final_expression.parquet"), emit: final_expression

    script:
    """
    get_counts_from_oarfish.py \\
        --min_reads ${min_reads} \\
        --min_n_sample ${min_n_sample}
    """
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

process make_db_files {
    conda "/scratch/nxu/SQANTI3/env"
    label "short_slurm_job"
    storeDir "nextflow_results/"

    input:
    path filtered_gtf

    output:
    path "sqanti3.db", emit: filtered_gtf_db

    script:
    """
    create_db_files.py --input $filtered_gtf \\
    --output sqanti3.db
    """
}
process isoquant {
    conda "/scratch/nxu/SQANTI3/env"
    label "long_slurm_job"
    storeDir "nextflow_results/isoquant_output/OUT"
    input:
    path annotation_gtf
    path ref_genome_fasta
    path ref_genome_index
    path aligned_bam
    path aligned_bam_bai

    output:
    path("*read_assignments.tsv.gz"), emit: read_assignments
    path("*transcript_counts.tsv"), emit: transcript_counts

    script:
    """
    export HOME="\$SCRATCH"
    isoquant.py \\
    --genedb $annotation_gtf \\
    --data_type pacbio \\
    --reference $ref_genome_fasta \\
    --index $ref_genome_index \\
    --complete_genedb \\
    --threads ${task.cpus} \\
    --bam ${aligned_bam} \\
    --polya_requirement never \\
    --no_model_construction
    """
}

process IsoAnnotLite {
    conda "/scratch/nxu/SQANTI3/env"
    input:
    path corrected_gtf
    path classification
    path junctions_txt
    output:
    path "IsoAnnotLite*"
    script:
    """
    python /scratch/nxu/SQANTI3/src/utilities/IsoAnnotLite_SQ3.py \\
    ${corrected_gtf} \\
    ${classification} \\
    ${junctions_txt}
    """
}

process runORFanage {
    label "short_slurm_job"
    conda "/scratch/nxu/astrocytes/env"
    storeDir "nextflow_results/orfanage"

    input:
    path ref_genome_fasta
    path final_sample_gtf

    output:
    path "orfanage_with_gene_id.gtf", emit: orfanage_gtf
    path "orfanage.stats"

    script:
    """
    gffread \\
        -g $ref_genome_fasta \\
        --adj-stop \\
        -T -F -J \\
        -o corrected.gtf \\
        $final_sample_gtf

    orfanage \\
        --reference $ref_genome_fasta \\
        --query corrected.gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        --minlen 50 \\
        --stats orfanage.stats \\
        ${params.annotation_gtf}

    gffread \\
        -g $ref_genome_fasta \\
        --adj-stop \\
        -T -F -J -C \\
        -o orfanage_with_gene_id.gtf \\
        orfanage_without_gene_id.gtf
    """
}

process fixORFanageFormat {
    label "short_slurm_job"
    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"
    storeDir "nextflow_results/orfanage"

    input:
    path ref_genome_fasta
    path orfanage_gtf

    output:
    path "orfanage.gtf"
    
    script:
    """
    agat_sp_add_start_and_stop.pl --gff $orfanage_gtf --fasta $ref_genome_fasta --out "added_codons_orfanage_with_gene_id.gff3"

    agat_convert_sp_gff2gtf.pl --gff "added_codons_orfanage_with_gene_id.gff3" -o "orfanage.gtf" --gtf_version 3
    """
}

process salmon_index {
    module "StdEnv/2023:gcc/12.3:openmpi/4.1.5:salmon/1.10.2"
    label "short_slurm_job"
    storeDir "nextflow_results/prepare/salmon_index"
    input:
    path transcriptome_fasta

    output:
    path "salmon_index"

    script:
    """
    salmon index -t $transcriptome_fasta -i salmon_index -k 31
    """
}

process salmon_quant {
    module "StdEnv/2023:gcc/12.3:openmpi/4.1.5:salmon/1.10.2"
    label "short_slurm_job"
    storeDir "nextflow_results/quantify/salmon_quant"
    input:
    path salmon_index
    tuple val(sample_id), path(fastq_files)

    output:
    path("${sample_id}_salmon_quant"), emit: salmon_quant_dir

    script:
    """
    salmon quant -i $salmon_index \\
    -l A \\
    -1 ${fastq_files[0]} \\
    -2 ${fastq_files[1]} \\
    -p ${task.cpus} \\
    --validateMappings \\
    -o ${sample_id}_salmon_quant
    """
}

process translateORFs {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage"

    input:
    path ref_genome_fasta
    path orfanage_gtf

    output:
    path "orfanage_proteins.fasta"

    script:
    """
    gffread -y orfanage_proteins.fasta \\
        -g $ref_genome_fasta \\
        $orfanage_gtf
    """
}

process format_gtf_for_ribotie {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/orfanage"
    input:
    path orfanage_gtf
    path final_classification
    path annotation_gtf

    output:
    path("orfanage_numbered_exons.gtf")

    script:
    """
    format_gtf_for_ribotie.py \\
    --input_gtf ${orfanage_gtf} \\
    --final_classification ${final_classification} \\
    --annotation_gtf ${annotation_gtf} \\
    --output_gtf orfanage_numbered_exons.gtf
    """
}

workflow {
    channel.fromPath("data/long_read/pacbio/*/*/*/hifi_reads/*.hifi_reads.bcM0001.bam")
        .map { file -> 
            def sample_id = file.simpleName
            return tuple(sample_id, file) 
        }
        .set { hifi_bams }
    skera(hifi_bams, params.kinnex_adapters)
    lima(skera.out, params.isoseq_primers, params.biosamples_csv)
    refine(lima.out.demultiplexed_bam, params.isoseq_primers)
    refine.out.flnc_bam.map { id, file -> file }.collect().set { flnc_bams }
    fofn(flnc_bams)
    cluster(fofn.out.flatten(), flnc_bams.collect())
    merge_flnc_bams(flnc_bams.collect())
    convert_flnc_bam_to_fastqz(merge_flnc_bams.out.flatten())
    minimap2_genome(convert_flnc_bam_to_fastqz.out)
    // bambu(minimap2_genome.out.collect())
    pbmm2(params.ref_genome_fasta, cluster.out)
    merge_aligned_bams(pbmm2.out.aligned_bam.collect())
    collapse_and_sort(merge_aligned_bams.out)
    // isoquant(params.annotation_gtf, params.ref_genome_fasta, params.ref_genome_index, pbmm2.out.aligned_bam.collect(), pbmm2.out.aligned_bam_bai.collect())
    channel.fromFilePairs(params.input_rna_fastq).set { short_read_fastqs }
    star_sr_genome(params.star_genomeDir, short_read_fastqs, params.annotation_gtf)
    sqanti_qc_isoseq(collapse_and_sort.out, params.annotation_gtf, params.ref_genome_fasta, params.refTSS, params.polyA_motif_list, star_sr_genome.out.star_aligned_bam.collect(), star_sr_genome.out.star_sj_tab.collect())
    def isoseq_corrected_gtf = sqanti_qc_isoseq.out
        .map { dir -> dir / "sqanti_qc_results_corrected.gtf" }
    def isoseq_classification = sqanti_qc_isoseq.out
        .map { dir -> dir / "sqanti_qc_results_classification.txt" }
    sqanti_filter_isoseq(isoseq_corrected_gtf, isoseq_classification)
    def filtered_gtf = sqanti_filter_isoseq.out
        .map {dir -> dir / "default.filtered.gtf"}
    def filtered_classification = sqanti_filter_isoseq.out
        .map {dir -> dir / "default_RulesFilter_result_classification.txt"}        
    extract_transcriptome(filtered_gtf)
    minimap2_transcriptome(extract_transcriptome.out, convert_flnc_bam_to_fastqz.out)
    run_oarfish(minimap2_transcriptome.out)
    isoform_exp_filter_params = channel.of(
            ['low_stringency', 1, 2],
            ['high_stringency', 5, 3]
        )
    filter_by_expression_input_ch = run_oarfish.out.collect().map { [it] }
        .combine(filtered_classification)
        .combine(filtered_gtf)
        .collect()
        .combine(isoform_exp_filter_params)
    filter_by_expression(filter_by_expression_input_ch)
    runORFanage(params.ref_genome_fasta, filter_by_expression.out.final_transcripts_gtf)
    fixORFanageFormat(params.ref_genome_fasta, runORFanage.out.orfanage_gtf)
    translateORFs(params.ref_genome_fasta, fixORFanageFormat.out)
    salmon_index(extract_transcriptome.out)
    channel.fromPath(params.riboseq_unmapped_to_contaminants).set{ riboseq_unmapped_to_contaminants }
    star_riboseq_custom(params.star_genomeDir, riboseq_unmapped_to_contaminants, filter_by_expression.out.final_transcripts_gtf, "custom")
    star_riboseq_gencode(params.star_genomeDir, riboseq_unmapped_to_contaminants, params.annotation_gtf, "gencode")
    
    // Split bedGraph files based on sample prefix
    star_riboseq_custom.out.bedGraph_UniqueMultiple
        .branch {
            unstim: it.name =~ /^merged_astro_[AB]_/
            stim: it.name =~ /^(merged_astro_C_|Astro_D_)/
        }
        .set { custom_bedgraphs }
    
    star_riboseq_gencode.out.bedGraph_UniqueMultiple
        .branch {
            unstim: it.name =~ /^merged_astro_[AB]_/
            stim: it.name =~ /^(merged_astro_C_|Astro_D_)/
        }
        .set { gencode_bedgraphs }
    
    // Merge unstimulated samples (A and B)
    merge_bg_and_convert_to_bw_custom_unstim(
        custom_bedgraphs.unstim.collect(),
        params.chrom_sizes,
        "custom_unstim"
    )
    
    // Merge stimulated samples (C and D)
    merge_bg_and_convert_to_bw_custom_stim(
        custom_bedgraphs.stim.collect(),
        params.chrom_sizes,
        "custom_stim"
    )
    
    format_gtf_for_ribotie(fixORFanageFormat.out, filter_by_expression.out.final_classification, params.annotation_gtf)
    // make_db_files(sqanti_filter.out.filtered_gtf)

    // sqanti_qc_bambu(bambu.out.supported_tx_gtf, params.annotation_gtf, params.ref_genome_fasta, params.refTSS, params.polyA_motif_list, star.out.star_aligned_bam.collect(), star.out.star_sj_tab.collect())
    // def bambu_corrected_gtf = sqanti_qc_bambu.out
    //     .map { dir -> dir / "sqanti_qc_results_corrected.gtf" }
    // def bambu_classification = sqanti_qc_bambu.out
    //     .map { dir -> dir / "sqanti_qc_results_classification.gtf" }
    // sqanti_filter_bambu(bambu_corrected_gtf, bambu_classification)
}