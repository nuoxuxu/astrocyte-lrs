#!/usr/bin/env Rscript
# Vignette figure: TCF4 x SCZ colocalization (PP4 = 0.977)
# Shows the novel astrocyte coding isoform (PB.19585.198) alongside the
# GENCODE MANE-Select canonical isoform (ENST00000354452.8 / TCF4-201),
# with the lead sQTL variant rs1319637 (chr18:55,540,094) marked in red.
# The novel junction (chr18:55,585,352-55,588,037) exon-skips a canonical
# exon (55,587,045-55,587,136) present in TCF4-201.
#
# Outputs:
#   figures/tcf4_scz_vignette.pdf       -- zoomed to 3' junction region
#   figures/tcf4_scz_vignette_full.pdf  -- full canonical isoform length
#
# Run from repo root:
#   Rscript scripts/aim_2/plot_tcf4_scz_vignette.R

library(ggtranscript)
library(ggplot2)
library(dplyr)
library(patchwork)

OUT_ZOOM  <- "figures/tcf4_scz_vignette.pdf"
OUT_FULL  <- "figures/tcf4_scz_vignette_full.pdf"
SNP_POS   <- 55540094   # rs1319637 (G->A, SCZ risk allele)
ZOOM_MIN  <- 55585000
ZOOM_MAX  <- 55588500

# ── Exon tables ───────────────────────────────────────────────────────────────

# GENCODE v47 ENST00000354452.8 (TCF4-201, MANE Select, minus strand)
gencode_exons <- data.frame(
  seqnames = "chr18",
  start = c(
    55222185, 55228221, 55228847, 55232509, 55234548, 55254497,
    55257315, 55259949, 55261466, 55269831, 55275619, 55279551,
    55350359, 55350874, 55403454, 55461019, 55464076,
    55585280, 55587045, 55588038
  ),
  end = c(
    55228030, 55228361, 55229076, 55232671, 55234683, 55254700,
    55257391, 55260027, 55261533, 55269963, 55275752, 55279656,
    55350408, 55351003, 55403518, 55461115, 55464137,
    55585352, 55587136, 55588192
  ),
  strand     = "-",
  transcript = "TCF4-201  (GENCODE v47, MANE Select)",
  exon_type  = "canonical"
)

# Novel astrocyte isoform PB.19585.198 (sqanti3 filtered, minus strand)
novel_exons <- data.frame(
  seqnames = "chr18",
  start = c(
    55227649, 55228221, 55228847, 55232521, 55234548, 55254497,
    55257315, 55259949, 55261466, 55269831, 55275619, 55279551,
    55350359, 55350874, 55403454, 55461019, 55464076,
    55585280, 55588038
  ),
  end = c(
    55228030, 55228361, 55229076, 55232671, 55234683, 55254700,
    55257391, 55260027, 55261533, 55269963, 55275752, 55279656,
    55350408, 55351003, 55403518, 55461115, 55464137,
    55585352, 55588161
  ),
  strand     = "-",
  transcript = "PB.19585.198  (novel astrocyte isoform)",
  exon_type  = "canonical"
)

# Tag exons flanking the novel junction
novel_exons <- novel_exons |>
  mutate(exon_type = if_else(
    start %in% c(55585280, 55588038),
    "junction", exon_type
  ))

all_exons <- bind_rows(gencode_exons, novel_exons) |>
  mutate(transcript = factor(transcript, levels = c(
    "TCF4-201  (GENCODE v47, MANE Select)",
    "PB.19585.198  (novel astrocyte isoform)"
  )))

all_introns <- to_intron(all_exons, group_var = "transcript") |>
  mutate(transcript = factor(transcript, levels = levels(all_exons$transcript)))

# ── CDS tables (thick boxes drawn on top of thin UTR exon boxes) ──────────────

# GENCODE ENST00000354452.8 CDS (from gencode.v47.annotation.gtf, feature==CDS)
gencode_cds <- data.frame(
  seqnames = "chr18",
  start = c(
    55228228, 55228847, 55232509, 55234548, 55254497,
    55257315, 55259949, 55261466, 55269831, 55275619, 55279551,
    55350359, 55350874, 55403454, 55461019, 55464076,
    55585280, 55587045
  ),
  end = c(
    55228361, 55229076, 55232671, 55234683, 55254700,
    55257391, 55260027, 55261533, 55269963, 55275752, 55279656,
    55350408, 55351003, 55403518, 55461115, 55464137,
    55585352, 55587116
  ),
  strand     = "-",
  transcript = "TCF4-201  (GENCODE v47, MANE Select)",
  exon_type  = "canonical"
)

# PB.19585.198 CDS (from orfanage/minlen/mid_stringency/orfanage.gtf, feature==CDS)
# Start codon at 55585350-55585352; 5' UTR is exon 55588038-55588161 (entire)
# 3' UTR: exon 55227649-55228030 (entire) + 55228221-55228224 (partial)
novel_cds <- data.frame(
  seqnames = "chr18",
  start = c(
    55228228, 55228847, 55232521, 55234548, 55254497,
    55257315, 55259949, 55261466, 55269831, 55275619, 55279551,
    55350359, 55350874, 55403454, 55461019, 55464076,
    55585280
  ),
  end = c(
    55228361, 55229076, 55232671, 55234683, 55254700,
    55257391, 55260027, 55261533, 55269963, 55275752, 55279656,
    55350408, 55351003, 55403518, 55461115, 55464137,
    55585352
  ),
  strand     = "-",
  transcript = "PB.19585.198  (novel astrocyte isoform)",
  exon_type  = "canonical"
)

# Tag junction-flanking CDS exon (55585280-55585352 is CDS in the novel isoform)
novel_cds <- novel_cds |>
  mutate(exon_type = if_else(start == 55585280, "junction", exon_type))

all_cds <- bind_rows(gencode_cds, novel_cds) |>
  mutate(transcript = factor(transcript, levels = levels(all_exons$transcript)))

# ── Shared base layers ────────────────────────────────────────────────────────

shared_theme <- theme_bw(base_size = 11) + theme(
  legend.position    = "bottom",
  panel.grid.major.y = element_blank(),
  panel.grid.minor   = element_blank(),
  axis.text.y        = element_text(size = 10),
  plot.title         = element_text(face = "bold", size = 12),
  plot.subtitle      = element_text(size = 9, colour = "grey40")
)

shared_scales <- list(
  scale_fill_manual(
    values = c(canonical = "grey50", junction = "#2166ac"),
    labels = c(canonical = "exon", junction  = "novel-junction exon"),
    name   = NULL
  ),
  scale_x_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    name   = "chr18 position (GRCh38)"
  ),
  labs(y = NULL)
)

base_plot <- ggplot(all_exons, aes(xstart = start, xend = end, y = transcript)) +
  geom_range(aes(fill = exon_type), height = 0.15) +       # thin: UTR extent
  geom_range(data = all_cds, aes(fill = exon_type),
             height = 0.30) +                               # thick: CDS only
  geom_intron(data = all_introns, aes(strand = strand),
              arrow.min.intron.length = 1e5) +
  geom_vline(xintercept = SNP_POS, colour = "red", linewidth = 0.9) +
  annotate("rect",
    xmin = 55585280, xmax = 55588192, ymin = 0.45, ymax = 2.55,
    fill = NA, colour = "grey40", linewidth = 0.4, linetype = "dashed"
  ) +
  shared_scales +
  shared_theme

# ── Zoomed version (3' junction region) ──────────────────────────────────────

p_zoom <- base_plot +
  annotate("text",
    x = SNP_POS + 1800, y = 2.55,
    label = "rs1319637\n(SCZ risk allele)",
    colour = "red", size = 3, hjust = 0, lineheight = 0.9
  ) +
  coord_cartesian(xlim = c(ZOOM_MIN, ZOOM_MAX)) +
  labs(
    title    = expression(italic("TCF4") ~ "x SCZ  --  novel astrocyte coding isoform (PP4 = 0.977)"),
    subtitle = paste0(
      "Zoomed to 3' junction region  |  ",
      "novel junction skips canonical exon 55,587,045-55,587,136  |  ",
      "gene reads right to left (minus strand)"
    )
  )

# ── Full-length version ───────────────────────────────────────────────────────

full_range  <- 55588192 - 55222185          # ~366 kb
label_nudge <- full_range * 0.01            # 1% of range

p_full <- base_plot +
  annotate("text",
    x = SNP_POS + label_nudge, y = 2.55,
    label = "rs1319637\n(SCZ risk allele)",
    colour = "red", size = 3, hjust = 0, lineheight = 0.9
  ) +
  labs(
    title    = expression(italic("TCF4") ~ "x SCZ  --  novel astrocyte coding isoform (PP4 = 0.977)"),
    subtitle = paste0(
      "Full canonical isoform length (~366 kb)  |  ",
      "dashed box = novel junction region  |  ",
      "gene reads right to left (minus strand)"
    )
  )

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(OUT_ZOOM, p_zoom, width = 10, height = 4, device = pdf)
message("Written: ", OUT_ZOOM)

ggsave(OUT_FULL, p_full, width = 12, height = 4, device = pdf)
message("Written: ", OUT_FULL)
