process fragpipe {
    // Apptainer is invoked directly (not via the `container` directive) because FragPipe
    // needs a custom --home and writable binds over its install's cache/ and jobs/ dirs.
    module "apptainer"
    label "mid_slurm_job"
    storeDir "nextflow_results/proteomics"

    input:
    path(workflow_file)
    path(manifest)
    path(database)
    path(fragpipe_sif)

    output:
    path("fragpipe"),                                         emit: workdir
    path("fragpipe/dia-quant-output/report.pg_matrix.tsv"),   emit: pg_matrix
    path("fragpipe/peptide.tsv"),                             emit: peptides
    path("fragpipe/psm.tsv"),                                 emit: psms

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

    # Point the workflow at the staged protein database
    sed "s|^database.db-path=.*|database.db-path=\$(readlink -f $database)|" $workflow_file > run.workflow

    mkdir -p fragpipe
    apptainer run \\
        --home ${params.fragpipe_home}:/home/nxu \\
        -B /home/nxu/tools:/home/nxu/tools \\
        -B /home/nxu/SCRATCH:/home/nxu/SCRATCH \\
        -B /scratch/nxu \\
        -B /cvmfs/soft.computecanada.ca \\
        -B ${params.fragpipe_writable}/cache:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/cache \\
        -B ${params.fragpipe_writable}/jobs:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/jobs \\
        $fragpipe_sif \\
        /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe \\
        --headless \\
        --workflow \$PWD/run.workflow \\
        --manifest \$(readlink -f $manifest) \\
        --workdir \$PWD/fragpipe \\
        --threads ${task.cpus} \\
        --config-tools-folder ${params.fragpipe_tools} \\
        --config-python /bin/python3
    """
}

workflow FRAGPIPE {
    take:
    workflow_file   // FragPipe .workflow file (e.g. assets/diaPASEF.workflow)
    manifest        // FragPipe .fp-manifest listing raw .d files and experiment/condition
    database        // Protein FASTA with decoys (overrides database.db-path in the workflow)
    fragpipe_sif    // FragPipe Apptainer image

    main:
    fragpipe(workflow_file, manifest, database, fragpipe_sif)

    emit:
    workdir   = fragpipe.out.workdir
    pg_matrix = fragpipe.out.pg_matrix
    peptides  = fragpipe.out.peptides
    psms      = fragpipe.out.psms
}
