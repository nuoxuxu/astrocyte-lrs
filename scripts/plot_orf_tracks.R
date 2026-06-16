#!/usr/bin/env Rscript
# Plot two-track gene annotation figure for a given ORF_id:
#   Track 1 – single input transcript from the collaborator (RiboTIE) GTF
#   Track 2 – GENCODE v47 protein-coding transcripts overlapping the region
# Zoomed to the ORF genomic location from summary_table_gencode_orftype_mismatch.tsv.
#
# Usage (interactive):
#   source("scripts/plot_orf_tracks.R")
#   load_track_data()                          # run once; loads GTFs into cache
#   plot_orf_tracks("PB.100.1_78")
#   ggsave("fig.pdf", plot_orf_tracks("PB.100.1_78"), width = 10, height = 6)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtranscript)
  library(rtracklayer)
  library(patchwork)
  library(glue)
})

# ── Default paths (relative to project root) ─────────────────────────────────
# .MISMATCH_TABLE   <- "nextflow_results/summary_table/summary_table.tsv"
.MISMATCH_TABLE   <- "nextflow_results/summary_table/summary_table_gencode_orftype_mismatch.tsv"
.COLLABORATOR_GTF <- "from_collaborator/filtered_output.gtf"
.GENCODE_GTF      <- "/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf"

# ── Private cache (populated by load_track_data) ─────────────────────────────
.cache <- new.env(parent = emptyenv())

# ── load_track_data() ─────────────────────────────────────────────────────────
#' Load GTFs and the mismatch table into memory once.
#' Must be called before plot_orf_tracks().
load_track_data <- function(
    mismatch_table   = .MISMATCH_TABLE,
    collaborator_gtf = .COLLABORATOR_GTF,
    gencode_gtf      = .GENCODE_GTF
) {
  message("Loading mismatch table ...")
  .cache$tbl <- read.delim(mismatch_table, stringsAsFactors = FALSE, check.names = FALSE)

  message("Loading collaborator GTF ...")
  .cache$collab <- as.data.frame(import(collaborator_gtf))

  message("Loading GENCODE GTF (large file, may take ~1-2 min) ...")
  gencode_full <- as.data.frame(import(gencode_gtf))

  .cache$gencode <- gencode_full %>%
    filter(!is.na(transcript_type), transcript_type == "protein_coding") %>%
    mutate(y_label = coalesce(transcript_name, transcript_id))

  # Set of Ensembl_canonical transcript IDs. rtracklayer collapses repeated
  # `tag` attributes (keeping only one), so the canonical flag is unreliable in
  # the imported data frame — derive it from the raw GTF text instead.
  message("Identifying Ensembl_canonical transcripts ...")
  canon_lines <- system(
    sprintf("grep Ensembl_canonical %s | grep -oE 'transcript_id \"[^\"]+\"'",
            shQuote(gencode_gtf)),
    intern = TRUE)
  .cache$canonical_tx <- unique(sub('transcript_id "([^"]+)"', "\\1", canon_lines))

  # Set of all GENCODE splice junctions (across every transcript biotype) used to
  # decide whether a collaborator junction is novel. A junction is keyed by
  # (seqnames, donor, acceptor) where, matching ggtranscript::to_intron, the
  # intron spans (previous exon end, next exon start).
  message("Computing GENCODE splice junctions for novelty annotation ...")
  .cache$gencode_junctions <- gencode_full %>%
    filter(type == "exon") %>%
    arrange(transcript_id, start) %>%
    group_by(transcript_id) %>%
    mutate(j_start = end, j_end = dplyr::lead(start)) %>%
    ungroup() %>%
    filter(!is.na(j_end)) %>%
    distinct(seqnames, j_start, j_end)

  message("Done. Call plot_orf_tracks(orf_id) to plot.")
  invisible(.cache)
}

# ── plot_orf_tracks() ─────────────────────────────────────────────────────────
#' @param orf_id  ORF_id from the mismatch table, e.g. "PB.100.1_78"
#' @param padding bp padding around the ORF span; NULL → max(1 kb, 20 % of span)
plot_orf_tracks <- function(orf_id, padding = NULL) {

  if (!exists("tbl", envir = .cache))
    stop("Run load_track_data() first.")

  tbl     <- .cache$tbl
  collab  <- .cache$collab
  gencode <- .cache$gencode

  # ── 1. Look up ORF ──────────────────────────────────────────────────────────
  row <- tbl[tbl$ORF_id == orf_id, , drop = FALSE]
  if (nrow(row) == 0)
    stop(sprintf("ORF_id '%s' not found in mismatch table", orf_id))
  row <- row[1L, ]

  tx_id     <- row$transcript_id
  chrom     <- sub(":.*", "", row$location)
  coords    <- as.integer(strsplit(sub(".*:", "", row$location), "-")[[1]])
  orf_start <- coords[1L]
  orf_end   <- coords[2L]

  if (is.null(padding))
    padding <- max(1000L, as.integer((orf_end - orf_start) * 0.2))
  view_start <- orf_start - padding
  view_end   <- orf_end   + padding

  message(sprintf("%s  |  %s:%d-%d  |  view: %d-%d",
                  orf_id, chrom, orf_start, orf_end, view_start, view_end))

  # ── 2. Track 1: single transcript from collaborator GTF ─────────────────────
  # In the collaborator GTF, transcript_id == ORF_id (e.g. "PB.100.1_78"),
  # which differs from transcript_id in the ribotie CSV ("PB.100.1").
  tx_features  <- collab %>% filter(transcript_id == orf_id)
  collab_exons <- tx_features %>% filter(type == "exon") %>% mutate(y_label = transcript_id)
  collab_cds   <- tx_features %>% filter(type == "CDS")  %>% mutate(y_label = transcript_id)
  target_cds   <- collab_cds  %>% filter(!is.na(ORF_id), ORF_id == orf_id)

  # Flag splice junctions of the collaborator transcript that are novel, i.e. the
  # intron (donor, acceptor) is absent from every GENCODE transcript.
  novel_junctions <- NULL
  if (nrow(collab_exons) > 1) {
    collab_introns <- to_intron(collab_exons, group_var = "y_label")
    if (nrow(collab_introns) > 0) {
      gj   <- .cache$gencode_junctions
      gkey <- paste(gj$seqnames, gj$j_start, gj$j_end)
      ikey <- paste(collab_introns$seqnames, collab_introns$start, collab_introns$end)
      novel_junctions <- collab_introns[!(ikey %in% gkey), , drop = FALSE]
    }
  }
  n_novel <- if (is.null(novel_junctions)) 0L else nrow(novel_junctions)

  # ── 3. Track 2: canonical GENCODE transcript(s) overlapping the view region ───
  # Restrict to Ensembl_canonical isoform(s) only (see .cache$canonical_tx).
  # Include all exons of the kept transcripts so intron arrows are drawn
  # correctly even when exons extend outside the view.
  tx_in_region <- gencode %>%
    filter(seqnames == chrom, start <= view_end, end >= view_start,
           transcript_id %in% .cache$canonical_tx) %>%
    pull(transcript_id) %>%
    unique()

  gencode_all   <- gencode %>% filter(transcript_id %in% tx_in_region)
  gencode_exons <- gencode_all %>% filter(type == "exon")
  gencode_cds   <- gencode_all %>% filter(type == "CDS")

  # ── 4. Panel builder ─────────────────────────────────────────────────────────
  make_panel <- function(exons, cds, title, highlight = NULL,
                         novel_junctions = NULL, subtitle = NULL) {
    if (nrow(exons) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No features in region",
                   size = 4, colour = "grey40") +
          labs(title = title) + theme_void(base_size = 9)
      )
    }

    introns <- to_intron(exons, group_var = "y_label")

    p <- ggplot(exons, aes(xstart = start, xend = end, y = y_label)) +
      geom_range(height = 0.25, fill = "grey80", colour = "grey50", linewidth = 0.3) +
      geom_intron(data = introns, aes(strand = strand),
                  arrow.min.intron.length = 200) +
      coord_cartesian(xlim = c(view_start, view_end)) +
      scale_x_continuous(labels = scales::label_comma()) +
      labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        axis.text.y        = element_text(size = 7),
        plot.title         = element_text(size = 9, face = "bold"),
        plot.subtitle      = element_text(size = 7, colour = "firebrick4")
      )

    if (nrow(cds) > 0)
      p <- p + geom_range(data = cds, height = 0.5,
                          fill = "steelblue", colour = "steelblue4", linewidth = 0.3)

    if (!is.null(highlight) && nrow(highlight) > 0)
      p <- p + geom_range(data = highlight, height = 0.5,
                          fill = "firebrick", colour = "firebrick4", linewidth = 0.3)

    # Overlay novel splice junctions (introns not present in GENCODE) in red.
    if (!is.null(novel_junctions) && nrow(novel_junctions) > 0)
      p <- p + geom_intron(data = novel_junctions, aes(strand = strand),
                           colour = "firebrick", linewidth = 0.8,
                           arrow.min.intron.length = 200)
    p
  }

  # ── 5. Assemble ──────────────────────────────────────────────────────────────
  p_collab <- make_panel(
    exons     = collab_exons,
    cds       = collab_cds,
    title     = sprintf("Collaborator (RiboTIE)  •  %s  •  ORF_type (orfanage): %s",
                        orf_id, row$ORF_type_orfanage),
    highlight = if (nrow(target_cds) > 0) target_cds else NULL,
    novel_junctions = novel_junctions,
    subtitle  = if (n_novel > 0)
                  sprintf("Red introns: %d novel splice junction%s not in GENCODE",
                          n_novel, if (n_novel == 1) "" else "s")
                else "No novel splice junctions (all introns present in GENCODE)"
  )

  p_gencode <- make_panel(
    exons = gencode_exons,
    cds   = gencode_cds,
    title = sprintf("GENCODE v47 (canonical)  •  %s:%d–%d (%s)  •  ORF_type (GENCODE): %s",
                    chrom, orf_start, orf_end, row$strand, row$ORF_type_gencode)
  )

  p_collab / p_gencode +
    plot_layout(heights = c(1, 1)) &
    theme(plot.margin = margin(4, 8, 4, 8))
}

load_track_data()

isoform_of_interest <- "PB.20515.11_56"

plot_orf_tracks(isoform_of_interest)
ggsave(glue("figures/inspect_orf_type_calling/{isoform_of_interest}.pdf"), plot_orf_tracks(isoform_of_interest), width = 10, height = 6)
print(glue("Saved figure for figures/inspect_orf_type_calling/{isoform_of_interest}.pdf"))
