#!/usr/bin/env Rscript
# Optimized setup: Pre-filter to only used transcripts to reduce memory footprint
# This version keeps only data relevant to the RiboTIE transcripts shown in the app

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

cat("Loading and optimizing GTF files...\n")

# Load RiboTIE first - this is the base set of transcripts we care about
cat("  Loading RiboTIE GTF...\n")
ribotie <- import(ribotie_gtf) %>%
    as_tibble() %>%
    select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
    mutate(source = "RiboTIE")

# Get unique RiboTIE transcript IDs for filtering
ribotie_tx_ids <- unique(ribotie$transcript_id)
cat("    Found", length(ribotie_tx_ids), "unique RiboTIE transcripts\n")

# Save RiboTIE (it's already small)
saveRDS(ribotie, file.path(data_dir, "ribotie.rds"))
cat("  ✓ Saved to data/ribotie.rds\n")

# Load ORFanage - filter to only RiboTIE templates + transcripts
cat("  Loading ORFanage GTF...\n")
orfanage_full <- import(orfanage_gtf) %>% as_tibble()

# Get mapping of RiboTIE -> ORFanage templates
pbid_to_orfanage_template <- orfanage_full %>%
    filter(type == "transcript") %>%
    select(transcript_id, orfanage_template)

orfanage_template_ids <- unique(pbid_to_orfanage_template$transcript_id)
cat("    Found", length(orfanage_template_ids), "ORFanage template transcripts\n")

# Keep only ORFanage data for: the templates + the RiboTIE-matched transcripts
orfanage <- orfanage_full %>%
    filter(transcript_id %in% c(orfanage_template_ids, ribotie_tx_ids)) %>%
    select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
    mutate(source = "ORFanage")

saveRDS(orfanage, file.path(data_dir, "orfanage.rds"))
cat("  ✓ Saved to data/orfanage.rds (filtered)\n")

# Save the template mapping
saveRDS(pbid_to_orfanage_template, file.path(data_dir, "pbid_to_orfanage_template.rds"))

# Load GENCODE - only keep genes that contain RiboTIE transcripts or their matches
cat("  Loading GENCODE GTF (largest file - this may take a moment)...\n")
gencode_full <- import(gencode_gtf) %>% as_tibble()

# Load the SQANTI mapping to get GENCODE matches
pbid_to_pr_transcripts <- read_tsv(pbid_mapping) %>%
    select(isoform_id, pr_transcripts)

# Get all GENCODE transcript IDs that are referenced
gencode_tx_ids <- c(
    pbid_to_pr_transcripts$pr_transcripts[pbid_to_pr_transcripts$pr_transcripts != "novel"],
    orfanage_template_ids
)
gencode_tx_ids <- unique(na.omit(gencode_tx_ids))
cat("    Found", length(gencode_tx_ids), "GENCODE transcripts to keep\n")

# Filter GENCODE to only these transcripts
gencode <- gencode_full %>%
    filter(transcript_id %in% gencode_tx_ids) %>%
    select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
    mutate(source = "GENCODE")

saveRDS(gencode, file.path(data_dir, "gencode.rds"))
cat("  ✓ Saved to data/gencode.rds (filtered to relevant transcripts)\n")

# Save the mapping
saveRDS(pbid_to_pr_transcripts, file.path(data_dir, "pbid_to_pr_transcripts.rds"))

cat("\nLoading and converting metadata files...\n")

# Load ORF type - filter to RiboTIE ORFs
cat("  Converting ORF type lookup...\n")
orf_type <- read_tsv(orf_type_file) %>%
    select(ORF_id, ORF_type_ORFanage)
orf_type_lookup <- setNames(orf_type$ORF_type_ORFanage, orf_type$ORF_id)
saveRDS(orf_type_lookup, file.path(data_dir, "orf_type_lookup.rds"))
cat("  ✓ Saved to data/orf_type_lookup.rds\n")

# Load pclass - filter to RiboTIE ORFs
cat("  Converting pclass lookup...\n")
pclass <- read_tsv(pclass_file) %>%
    select(transcript_id, pclass)
pclass_lookup <- setNames(pclass$pclass, pclass$transcript_id)
saveRDS(pclass_lookup, file.path(data_dir, "pclass_lookup.rds"))
cat("  ✓ Saved to data/pclass_lookup.rds\n")

# Summary
cat("\n=== Optimization Summary ===\n")
cat("RiboTIE transcripts: ", length(ribotie_tx_ids), "\n", sep = "")
cat("GENCODE filtered from all to: ", length(gencode_tx_ids), " transcripts\n", sep = "")
cat("ORFanage filtered to: ", nrow(orfanage) / max(table(orfanage$transcript_id)), " transcripts\n", sep = "")
cat("\n✓ All optimized data files created!\n")
cat("✓ Memory footprint significantly reduced\n")
cat("You can now run: shiny::runApp()\n")
