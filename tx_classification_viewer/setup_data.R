#!/usr/bin/env Rscript
# Setup script: Convert original data files to optimized RDS format
# Run once: Rscript setup_data.R

library(rtracklayer)
library(dplyr)
library(readr)

data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)

# Define paths to original files
BASE_PATH <- "/scratch/nxu/astrocytes"
gencode_gtf <- file.path(BASE_PATH, "data/gencode.v47.annotation.gtf")
ribotie_gtf <- file.path(BASE_PATH, "nextflow_results/translatome/supplemented_collaborator/filtered_output_fixed.gtf")
orfanage_gtf <- file.path(BASE_PATH, "nextflow_results/orfanage/minlen/orfanage.gtf")
pbid_mapping <- file.path(BASE_PATH, "nextflow_results/sqanti3_protein/test.predicted_proteome.best_ORF_SQANTI_classification.tsv")
orf_type_file <- file.path(BASE_PATH, "nextflow_results/quality/collaborator/orf_type_gencode.tsv")
pclass_file <- file.path(BASE_PATH, "nextflow_results/protein_classification/collaborator.predicted_proteome.best_ORF_summary.txt")

cat("Loading and converting GTF files...\n")

# Load GENCODE GTF
if (!file.exists(file.path(data_dir, "gencode.rds"))) {
    cat("  Converting GENCODE GTF...\n")
    gencode <- import(gencode_gtf) %>%
        as_tibble() %>%
        select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
        mutate(source = "GENCODE")
    saveRDS(gencode, file.path(data_dir, "gencode.rds"))
    cat("  ✓ Saved to data/gencode.rds\n")
}

# Load RiboTIE GTF
if (!file.exists(file.path(data_dir, "ribotie.rds"))) {
    cat("  Converting RiboTIE GTF...\n")
    ribotie <- import(ribotie_gtf) %>%
        as_tibble() %>%
        select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
        mutate(source = "RiboTIE")
    saveRDS(ribotie, file.path(data_dir, "ribotie.rds"))
    cat("  ✓ Saved to data/ribotie.rds\n")
}

# Load ORFanage GTF
if (!file.exists(file.path(data_dir, "orfanage.rds"))) {
    cat("  Converting ORFanage GTF...\n")
    orfanage <- import(orfanage_gtf) %>%
        as_tibble() %>%
        select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
        mutate(source = "ORFanage")
    saveRDS(orfanage, file.path(data_dir, "orfanage.rds"))
    cat("  ✓ Saved to data/orfanage.rds\n")
}

cat("\nLoading and converting mapping/metadata files...\n")

# Load pbid to pr_transcripts mapping
if (!file.exists(file.path(data_dir, "pbid_to_pr_transcripts.rds"))) {
    cat("  Converting pbid_to_pr_transcripts...\n")
    pbid_to_pr_transcripts <- read_tsv(pbid_mapping) %>%
        select(isoform_id, pr_transcripts)
    saveRDS(pbid_to_pr_transcripts, file.path(data_dir, "pbid_to_pr_transcripts.rds"))
    cat("  ✓ Saved to data/pbid_to_pr_transcripts.rds\n")
}

# Load pbid to orfanage_template mapping
if (!file.exists(file.path(data_dir, "pbid_to_orfanage_template.rds"))) {
    cat("  Converting pbid_to_orfanage_template...\n")
    pbid_to_orfanage_template <- import(orfanage_gtf) %>%
        as_tibble() %>%
        filter(type == "transcript") %>%
        select(transcript_id, orfanage_template)
    saveRDS(pbid_to_orfanage_template, file.path(data_dir, "pbid_to_orfanage_template.rds"))
    cat("  ✓ Saved to data/pbid_to_orfanage_template.rds\n")
}

# Load ORF type lookup
if (!file.exists(file.path(data_dir, "orf_type_lookup.rds"))) {
    cat("  Converting ORF type lookup...\n")
    orf_type <- read_tsv(orf_type_file) %>%
        select(ORF_id, ORF_type_ORFanage)
    orf_type_lookup <- setNames(orf_type$ORF_type_ORFanage, orf_type$ORF_id)
    saveRDS(orf_type_lookup, file.path(data_dir, "orf_type_lookup.rds"))
    cat("  ✓ Saved to data/orf_type_lookup.rds\n")
}

# Load pclass lookup
if (!file.exists(file.path(data_dir, "pclass_lookup.rds"))) {
    cat("  Converting pclass lookup...\n")
    pclass <- read_tsv(pclass_file) %>%
        select(transcript_id, pclass)
    pclass_lookup <- setNames(pclass$pclass, pclass$transcript_id)
    saveRDS(pclass_lookup, file.path(data_dir, "pclass_lookup.rds"))
    cat("  ✓ Saved to data/pclass_lookup.rds\n")
}

cat("\n✓ All data files converted successfully!\n")
cat("You can now run: shiny::runApp()\n")
