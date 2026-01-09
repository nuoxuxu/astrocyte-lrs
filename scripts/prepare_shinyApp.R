library(arrow)
library(dplyr)
library(rtracklayer)
library(tidyr)
library(readr)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC",
    "Other"                    = "Other"
)

gene_id_mapping <- rtracklayer::import("/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf") %>% 
  as_tibble() %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  dplyr::rename(associated_gene = gene_id)

classification <- read_parquet("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_classification.parquet") %>% 
  mutate(
    structural_category = if_else(structural_category %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"), structural_category, "Other"),
    structural_category = structural_category_labels[structural_category],
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
  ) %>% 
  left_join(
    gene_id_mapping,
    by = "associated_gene"
  ) %>%
  mutate(
    associated_gene = case_when(
      is.na(gene_name) ~ associated_gene,
      TRUE ~ gene_name
    ),
  ) %>% select(-gene_name)

lr_read_count <- read_parquet("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_expression.parquet")

# Get lr_log2_cpm
lr_log2_cpm <- lr_read_count %>% 
  mutate(
    across(
      where(is.numeric), 
      ~ log2((.x / sum(.x) * 1e6) + 1)
    )
  )

# Rename lr_log2_cpm column names from BC0X to Astro_X
mapping_df <- read.csv("data/primer_to_sample.csv")
indices <- match(colnames(lr_log2_cpm), mapping_df$Index.primer)
colnames(lr_log2_cpm) <- mapping_df$Sample_ID[indices]
lr_log2_cpm <- lr_log2_cpm %>%
  dplyr::rename(isoform = 1)

# Calculate mean over stim and unstim and add structural category
Stim_columns <- mapping_df %>% filter(condition == "Stim") %>% pull(Sample_ID)
Unstim_columns <- mapping_df %>% filter(condition == "Unstim") %>% pull(Sample_ID)

lr_log2_cpm <- lr_log2_cpm %>% 
  mutate(
    mean_Stim = rowMeans(pick(all_of(Stim_columns)), na.rm=TRUE),
    mean_Unstim = rowMeans(pick(all_of(Unstim_columns)), na.rm=TRUE)
  ) %>% 
  select(isoform, mean_Stim, mean_Unstim) %>%
  pivot_longer(
    cols = c("mean_Stim", "mean_Unstim"),
    names_to = "condition",
    values_to = "abundance"
  ) %>%
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = "isoform"
  )

orfanage_gtf <- rtracklayer::import("nextflow_results/orfanage/orfanage.gtf") %>% as.data.frame()
full_gtf <- rtracklayer::import("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_transcripts.gtf") %>% as.data.frame()
tx_no_CDS <- setdiff(pull(distinct(full_gtf, transcript_id), transcript_id), pull(distinct(orfanage_gtf, transcript_id), transcript_id))
gtf <- bind_rows(orfanage_gtf, full_gtf %>% filter(transcript_id %in% tx_no_CDS))

gtf <- gtf %>%
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = c("transcript_id" = "isoform")
  )

gtf <- gtf %>%
  select(
    -c(width, source, score, phase, gene_id, ID, orfanage_duplicity, orfanage_status, orfanage_template, orfanage_template_source, Parent, original_biotype)
  ) %>%
  filter(!(type %in% c("five_prime_utr", "three_prime_utr", "start_codon", "stop_codon", "gene")))

lr_log2_cpm %>% write.csv("transcript_vis_app/data/lr_log2_cpm.csv", row.names = FALSE)
gtf %>% write.csv("transcript_vis_app/data/gtf.csv", row.names = FALSE)