process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    module "python:gcc:arrow/19.0.1:rust"
    beforeScript 'source /scratch/nxu/astrocytes/pytorch/bin/activate'

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path ("quants/*")
    val quant_type
    val id
    val extra

    output:
    tuple val(meta), path("*tx2gene.tsv"), emit: tx2gene
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tx2gene.py'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
