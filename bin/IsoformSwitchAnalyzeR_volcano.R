#!/usr/bin/env Rscript
library(IsoformSwitchAnalyzeR)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(argparse)

parser <- ArgumentParser(description='Plot volcano plots for DGE, DTE and DTU results')
parser$add_argument('--switchlist', type='character', required=TRUE, help='Path to switch list object RDS file')
args <- parser$parse_args()
obj <- readRDS(args$switchlist)

# ── Volcano plots ─────────────────────────────────────────────────────────────

feat <- obj$isoformFeatures

ALPHA      <- 0.05
LFC_CUTOFF <- 1.0
DIF_CUTOFF <- 0.1

# Gene-level volcano (DGE)
gene_df <- feat %>%
    dplyr::select(gene_id, gene_name, gene_log2_fold_change, gene_q_value) %>%
    distinct() %>%
    mutate(
        neglog10q = -log10(pmax(gene_q_value, 1e-300)),
        sig = gene_q_value < ALPHA & abs(gene_log2_fold_change) >= LFC_CUTOFF,
        direction = case_when(
            sig & gene_log2_fold_change >  LFC_CUTOFF ~ "Up in Stim",
            sig & gene_log2_fold_change < -LFC_CUTOFF ~ "Down in Stim",
            TRUE                                       ~ "NS"
        )
    )

top_gene_ids <- bind_rows(
    gene_df %>% filter(direction == "Up in Stim")   %>% slice_max(gene_log2_fold_change, n = 20),
    gene_df %>% filter(direction == "Down in Stim")  %>% slice_min(gene_log2_fold_change, n = 20)
) %>% pull(gene_id)

gene_df <- gene_df %>%
    mutate(label = ifelse(gene_id %in% top_gene_ids & !is.na(gene_name), gene_name, NA_character_))

n_up   <- sum(gene_df$direction == "Up in Stim",   na.rm = TRUE)
n_down <- sum(gene_df$direction == "Down in Stim",  na.rm = TRUE)

pal <- c("Up in Stim" = "#e41a1c", "Down in Stim" = "#377eb8", "NS" = "grey70")

p_dge <- ggplot(gene_df, aes(x = gene_log2_fold_change, y = neglog10q, colour = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_text_repel(
        aes(label = label),
        size = 2.5, max.overlaps = Inf,
        min.segment.length = 0.2, box.padding = 0.3,
        na.rm = TRUE
    ) +
    scale_colour_manual(values = pal, name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("Up: ", n_up, " genes"),   colour = "#e41a1c", size = 3.5) +
    annotate("text", x =  -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Down: ", n_down, " genes"), colour = "#377eb8", size = 3.5) +
    labs(
        title    = "Differential Gene Expression  (Stim vs Unstim)",
        subtitle = paste0("DESeq2  |  padj < ", ALPHA, "  &  |log2FC| >= ", LFC_CUTOFF),
        x        = "log2 Fold Change",
        y        = "-log10(adjusted p-value)"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

# Isoform-level volcano (DTU)
iso_df <- feat %>%
    dplyr::select(isoform_id, gene_name, dIF, isoform_switch_q_value) %>%
    distinct() %>%
    mutate(
        neglog10q = -log10(pmax(isoform_switch_q_value, 1e-300)),
        sig = isoform_switch_q_value < ALPHA & abs(dIF) >= DIF_CUTOFF,
        direction = case_when(
            sig & dIF >  DIF_CUTOFF ~ "Increased in Stim",
            sig & dIF < -DIF_CUTOFF ~ "Decreased in Stim",
            TRUE                    ~ "NS"
        )
    )

top_iso_ids <- bind_rows(
    iso_df %>% filter(direction == "Increased in Stim") %>% slice_max(dIF, n = 20),
    iso_df %>% filter(direction == "Decreased in Stim") %>% slice_min(dIF, n = 20)
) %>% pull(isoform_id)

iso_df <- iso_df %>%
    mutate(label = ifelse(isoform_id %in% top_iso_ids & !is.na(gene_name), gene_name, NA_character_))

n_inc <- sum(iso_df$direction == "Increased in Stim", na.rm = TRUE)
n_dec <- sum(iso_df$direction == "Decreased in Stim", na.rm = TRUE)

pal2 <- c("Increased in Stim" = "#e41a1c", "Decreased in Stim" = "#377eb8", "NS" = "grey70")

p_dtu <- ggplot(iso_df, aes(x = dIF, y = neglog10q, colour = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-DIF_CUTOFF, DIF_CUTOFF), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_text_repel(
        aes(label = label),
        size = 2.5, max.overlaps = Inf,
        min.segment.length = 0.2, box.padding = 0.3,
        na.rm = TRUE
    ) +
    scale_colour_manual(values = pal2, name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("Increased: ", n_inc, " isoforms"), colour = "#e41a1c", size = 3.5) +
    annotate("text", x =  -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Decreased: ", n_dec, " isoforms"), colour = "#377eb8", size = 3.5) +
    labs(
        title    = "Differential Transcript Usage  (Stim vs Unstim)",
        subtitle = paste0("DEXSeq  |  q < ", ALPHA, "  &  |dIF| >= ", DIF_CUTOFF),
        x        = "Delta Isoform Fraction (dIF)",
        y        = "-log10(q-value)"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

# Isoform-level volcano (isoform DGE)
iso_dge_df <- feat %>%
    dplyr::select(isoform_id, gene_name, iso_log2_fold_change, iso_q_value) %>%
    distinct() %>%
    mutate(
        neglog10q = -log10(pmax(iso_q_value, 1e-300)),
        sig = iso_q_value < ALPHA & abs(iso_log2_fold_change) >= LFC_CUTOFF,
        direction = case_when(
            sig & iso_log2_fold_change >  LFC_CUTOFF ~ "Up in Stim",
            sig & iso_log2_fold_change < -LFC_CUTOFF ~ "Down in Stim",
            TRUE                                       ~ "NS"
        )
    )

top_iso_dge_ids <- bind_rows(
    iso_dge_df %>% filter(direction == "Up in Stim")   %>% slice_max(iso_log2_fold_change, n = 20),
    iso_dge_df %>% filter(direction == "Down in Stim")  %>% slice_min(iso_log2_fold_change, n = 20)
) %>% pull(isoform_id)

iso_dge_df <- iso_dge_df %>%
    mutate(label = ifelse(isoform_id %in% top_iso_dge_ids & !is.na(gene_name), gene_name, NA_character_))

n_iso_up   <- sum(iso_dge_df$direction == "Up in Stim",   na.rm = TRUE)
n_iso_down <- sum(iso_dge_df$direction == "Down in Stim",  na.rm = TRUE)

p_iso_dge <- ggplot(iso_dge_df, aes(x = iso_log2_fold_change, y = neglog10q, colour = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_text_repel(
        aes(label = label),
        size = 2.5, max.overlaps = Inf,
        min.segment.length = 0.2, box.padding = 0.3,
        na.rm = TRUE
    ) +
    scale_colour_manual(values = pal, name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("Up: ", n_iso_up, " isoforms"),   colour = "#e41a1c", size = 3.5) +
    annotate("text", x =  -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Down: ", n_iso_down, " isoforms"), colour = "#377eb8", size = 3.5) +
    labs(
        title    = "Differential Isoform Expression  (Stim vs Unstim)",
        subtitle = paste0("DESeq2  |  padj < ", ALPHA, "  &  |log2FC| >= ", LFC_CUTOFF),
        x        = "log2 Fold Change",
        y        = "-log10(adjusted p-value)"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave("volcano_DGE.pdf",     p_dge,     width = 7, height = 6)
ggsave("volcano_DTU.pdf",     p_dtu,     width = 7, height = 6)
ggsave("volcano_iso_DGE.pdf", p_iso_dge, width = 7, height = 6)

cat("Saved: volcano_DGE.pdf\n")
cat("Saved: volcano_DTU.pdf\n")
cat("Saved: volcano_iso_DGE.pdf\n")

cat("\nDGE summary (padj <", ALPHA, "& |LFC| >=", LFC_CUTOFF, "):\n")
cat("  Up in Stim:   ", n_up,   "\n")
cat("  Down in Stim: ", n_down, "\n")

cat("\nDTU summary (q <", ALPHA, "& |dIF| >=", DIF_CUTOFF, "):\n")
cat("  Increased in Stim: ", n_inc, "\n")
cat("  Decreased in Stim: ", n_dec, "\n")

cat("\nIsoform DGE summary (padj <", ALPHA, "& |LFC| >=", LFC_CUTOFF, "):\n")
cat("  Up in Stim:   ", n_iso_up,   "\n")
cat("  Down in Stim: ", n_iso_down, "\n")
