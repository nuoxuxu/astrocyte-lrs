#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(arrow)
library(stringr)

ALPHA      <- 0.05
DIF_CUTOFF <- 0.1

CAT_LABELS <- c(
    "full-splice_match"       = "FSM",
    "incomplete-splice_match" = "ISM",
    "novel_in_catalog"        = "NIC",
    "novel_not_in_catalog"    = "NNC",
    "fusion"                  = "Fusion",
    "genic"                   = "Genic",
    "antisense"               = "Antisense",
    "intergenic"              = "Intergenic",
    "genic_intron"            = "Genic intron"
)

CAT_ORDER <- c("FSM", "ISM", "NIC", "NNC", "Fusion", "Genic", "Antisense", "Intergenic", "Genic intron")

CAT_COLORS <- c(
    "FSM"         = "#009E73",
    "ISM"         = "#0072B2",
    "NIC"         = "#D55E00",
    "NNC"         = "#E69F00",
    "Fusion"      = "#CC79A7",
    "Genic"       = "#F0E442",
    "Antisense"   = "#000000",
    "Intergenic"  = "#000000",
    "Genic intron"= "#000000"
)

obj <- readRDS("nextflow_results/IsoformSwitchAnalyzeR/mid_stringency/IsoformSwitchAnalyzeR.rds")

classification <- read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/mid_stringency/final_classification.parquet"
) %>%
    dplyr::select(isoform, structural_category) %>%
    mutate(category = recode(structural_category, !!!CAT_LABELS))

# ISA isoform IDs have a _<number> suffix not present in the classification
feat <- obj$isoformFeatures %>%
    dplyr::select(isoform_id, dIF, isoform_switch_q_value) %>%
    distinct() %>%
    mutate(isoform = str_replace(isoform_id, "_[0-9]+$", "")) %>%
    left_join(classification, by = "isoform") %>%
    filter(!is.na(category))

all_df <- feat %>%
    dplyr::count(category, name = "n") %>%
    mutate(pct = n / sum(n) * 100, group = "All tested")

sig_df <- feat %>%
    filter(isoform_switch_q_value < ALPHA, abs(dIF) >= DIF_CUTOFF) %>%
    dplyr::count(category, name = "n") %>%
    mutate(pct = n / sum(n) * 100, group = "Differential (DTU)")

# Ensure all categories present in all_df also appear in sig_df (with 0) for
# consistent dodging
sig_df <- all_df %>%
    dplyr::select(category) %>%
    left_join(sig_df, by = "category") %>%
    mutate(
        n     = replace(n,   is.na(n),   0L),
        pct   = replace(pct, is.na(pct), 0),
        group = "Differential (DTU)"
    )

plot_df <- bind_rows(all_df, sig_df) %>%
    mutate(
        category = factor(category, levels = CAT_ORDER),
        group    = factor(group, levels = c("All tested", "Differential (DTU)"))
    )

n_all <- sum(all_df$n)
n_sig <- sum(sig_df$n)

# Keep only categories that appear in the data
present_cats <- levels(droplevels(plot_df$category[plot_df$n > 0]))
plot_df <- plot_df %>% filter(category %in% present_cats) %>%
    mutate(category = factor(category, levels = present_cats))

p <- ggplot(plot_df, aes(x = category, y = pct, fill = category, alpha = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = CAT_COLORS, guide = "none") +
    scale_alpha_manual(
        values = c("All tested" = 0.35, "Differential (DTU)" = 1),
        name   = NULL,
        labels = c(
            paste0("All tested (n=", n_all, ")"),
            paste0("Differential DTU (n=", n_sig, ")")
        )
    ) +
    labs(
        title    = "Structural category of differential transcripts (DTU)",
        subtitle = paste0("q < ", ALPHA, "  &  |dIF| >= ", DIF_CUTOFF, "  |  Stim vs Unstim"),
        x        = "Structural category",
        y        = "Percentage of transcripts (%)"
    ) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.text.x     = element_text(angle = 35, hjust = 1)
    )

ggsave("figures/dtu_structural_categories.pdf", p, width = 7, height = 5)
cat("Saved: figures/dtu_structural_categories.pdf\n")

cat("\nCategory breakdown of differential DTU transcripts:\n")
print(sig_df %>% arrange(desc(pct)) %>% dplyr::select(category, n, pct))
