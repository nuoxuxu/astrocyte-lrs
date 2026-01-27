#!/scratch/nxu/astrocytes/pytorch/bin/python
from transcript_transformer.processing import csv_to_gtf, construct_output_table, create_multiqc_reports
from collections import namedtuple
import numpy as np
import os
from pathlib import Path
import polars as pl
"We need to rerun the output table construction and gtf conversion steps because they needed gene_name column to run"

def add_gene_name(df, classification, annotation_gtf):
    from bin.src.utils import read_gtf
    import polars as pl
    id_to_gene_name = read_gtf(annotation_gtf, attributes=["gene_id", "gene_name"])\
        .filter(pl.col("feature")=="gene")\
        .select("gene_id", "gene_name")\
        .unique()
    classification = pl.read_parquet(classification)\
        .join(
            id_to_gene_name,
            left_on="associated_gene",
            right_on="gene_id",
            how="left"
        )
    return df.join(
        classification.select(["isoform", "gene_name"]),
        left_on="transcript_id", 
        right_on="isoform"
    )

Args = namedtuple("Args", [
    "conf", "debug", "data", "overwrite_data", "overwrite_models", "overwrite_preds",
    "pretrain", "gtf_path", "fa_path", "h5_path", "ribo_paths", "samples", "parallel",
    "no_backup", "backup_path", "offsets", "low_memory", "cores", "out_prefix",
    "prob_cutoff", "start_codons", "min_ORF_len", "include_invalid_TTS",
    "keep_duplicates", "return_ORF_coords", "no_correction", "distance",
    "num_workers", "max_memory", "accelerator", "strategy", "devices", "lr",
    "decay_rate", "warmup_steps", "max_epochs", "patience", "transfer_checkpoint",
    "folds", "val_frac", "test_frac", "cond", "strict_validation", "leaky_frac",
    "min_seq_len", "max_seq_len", "max_transcripts_per_batch", "use_seq",
    "input_type", "exp_path", "y_path", "seqn_path", "id_path", "num_tokens",
    "dim", "depth", "heads", "dim_head", "local_attn_heads", "nb_features",
    "feature_redraw_interval", "no_generalized_attention", "reversible",
    "ff_chunks", "use_scalenorm", "use_rezero", "ff_glu", "emb_dropout",
    "ff_dropout", "attn_dropout", "local_window_size", "mlm", "mask_frac",
    "rand_frac", "model_dir", "missing_models", "use_ribo", "grouped_ribo_ids",
])

args = Args(
    conf=["ribotie_isoseq_ORFanage.yml"],
    debug=False,
    data=False,
    overwrite_data=False,
    overwrite_models=False,
    overwrite_preds=False,
    pretrain=False,
    gtf_path="nextflow_results/orfanage/orfanage_numbered_exons.gtf",
    fa_path="data/GRCh38.primary_assembly.genome.fa",
    h5_path="ribotie_res/isoseq_ORFanage/custom.h5",
    ribo_paths={
        "Unstim": [
            ["astro_A", "nextflow_results/align/riboseq/merged_astro_A_Unmapped.Aligned.toTranscriptome.out.bam"],
            ["astro_B", "nextflow_results/align/riboseq/merged_astro_B_Unmapped.Aligned.toTranscriptome.out.bam"],
        ],
        "Stim": [
            ["astro_C", "nextflow_results/align/riboseq/merged_astro_C_Unmapped.Aligned.toTranscriptome.out.bam"],
            ["astro_D", "nextflow_results/align/riboseq/Astro_D_Unmapped.Aligned.toTranscriptome.out.bam"],
        ],
    },
    samples={"Unstim": ["astro_A", "astro_B"], "Stim": ["astro_C", "astro_D"]},
    parallel=False,
    no_backup=False,
    backup_path=None,
    offsets=None,
    low_memory=False,
    cores=6,
    out_prefix="ribotie_res/isoseq_ORFanage/custom",
    prob_cutoff=0.125,
    start_codons=".*TG$",
    min_ORF_len=15,
    include_invalid_TTS=False,
    keep_duplicates=False,
    return_ORF_coords=False,
    no_correction=False,
    distance=9,
    num_workers=2,
    max_memory=30000,
    accelerator="gpu",
    strategy="auto",
    devices=1,
    lr=0.001,
    decay_rate=0.96,
    warmup_steps=1500,
    max_epochs=60,
    patience=4,
    transfer_checkpoint=None,
    folds=None,
    val_frac=0.2,
    test_frac=0.5,
    cond={
        "global": {"transcript_len": lambda x: x},  # placeholder if you need a specific lambda
        "grouped": {
            "Unstim": {"num_reads": lambda x: x},
            "Stim": {"num_reads": lambda x: x},
        },
    },
    strict_validation=False,
    leaky_frac=0.05,
    min_seq_len=0,
    max_seq_len=30000,
    max_transcripts_per_batch=2000,
    use_seq=False,
    input_type="hdf5",
    exp_path="transcript",
    y_path="tis",
    seqn_path="seqname",
    id_path="transcript_id",
    num_tokens=8,
    dim=42,
    depth=8,
    heads=6,
    dim_head=16,
    local_attn_heads=4,
    nb_features=80,
    feature_redraw_interval=1000,
    no_generalized_attention=True,
    reversible=False,
    ff_chunks=1,
    use_scalenorm=False,
    use_rezero=False,
    ff_glu=False,
    emb_dropout=0.1,
    ff_dropout=0.1,
    attn_dropout=0.1,
    local_window_size=256,
    mlm=False,
    mask_frac=False,
    rand_frac=False,
    model_dir=None,
    missing_models=True,
    use_ribo=True,
    grouped_ribo_ids={"Unstim": ["astro_A", "astro_B"], "Stim": ["astro_C", "astro_D"]},
)
output_sets = [str(k) for k in args.grouped_ribo_ids.keys()]
for output in output_sets:
    out = np.load(f"{args.out_prefix}_{output}.npy", allow_pickle=True)
    out_prefix = f"{args.out_prefix}_{output}"
    df, df_filt, df_novel = construct_output_table(
        h5_path=args.h5_path,
        out_prefix=out_prefix,
        prob_cutoff=args.prob_cutoff,
        correction=not args.no_correction,
        dist=args.distance,
        start_codons=args.start_codons,
        min_ORF_len=args.min_ORF_len,
        remove_duplicates=not args.keep_duplicates,
        exclude_invalid_TTS=not args.include_invalid_TTS,
        ribo_output=out,
        grouped_ribo_ids=args.grouped_ribo_ids,
        parallel=args.parallel,
        return_ORF_coords=args.return_ORF_coords,
    )
    if df is not None:
        ids = ["ribotie_all", "ribotie"]
        names = ["RiboTIE_redundant", "RiboTIE"]
        paths = [out_prefix + ".redundant", out_prefix]
        multiqc_path = os.path.join(os.path.dirname(args.out_prefix), "multiqc")
        os.makedirs(multiqc_path, exist_ok=True)
        for df, id, name, path in zip([df, df_filt], ids, names, paths):
            classification = "nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_classification.parquet"
            annotation_gtf = "/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf"
            df = add_gene_name(df, classification, annotation_gtf)
            csv_to_gtf(args.h5_path, df, path, "RiboTIE")
            out = os.path.join(multiqc_path, os.path.basename(path))
            create_multiqc_reports(df, out, id, name)
        csv_to_gtf(args.h5_path, df_novel, out_prefix + ".novel", "RiboTIE")
    else:
        print("No positive predictions found")

    if args.pretrain:
        print(
            f"""\n
            !!! 
            Do not use the pre-trained predictions as is.
            RiboTIE is meant to fine-tune on individual samples after pre-training.
            Run RiboTIE without the --pretrain flag and with newly created yml file, e.g.,
            'ribotie {' '.join(args.conf)} {args.out_prefix}_pretrain_params.rt.yml'.
            !!!
            """
        )