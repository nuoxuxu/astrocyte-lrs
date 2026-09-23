library(Gviz)
library(GenomicFeatures)
library(GenomicRanges)
library(glue)
library(tidyr)
library(rtracklayer)
library(dplyr)
library(arrow)

# Settings
options(stringsAsFactors = FALSE)
options(Gviz.scheme = "myScheme")
options(ucscChromosomeNames = FALSE)

scheme <- getScheme()
scheme$GeneRegionTrack$col <- NULL
addScheme(scheme, "myScheme")

# Current gene
current_gene <- "RAB2A"
gene_id_to_plot <- "ENSG00000104388.15"
tx_to_plot <- "ENST00000262646.12"
transcript_id <- "PB.10505.6"
ORF_id <- "PB.10505.6_172"

# GENCODE track
# GENCODE_GRList <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf") %>%
#   makeTxDbFromGFF(format = "gtf") %>%
#   exonsBy(by = "tx", use.names = TRUE)

# gencode_gtf <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
#   as_tibble()

# gencode_track <- GENCODE_GRList[tx_to_plot, ] %>%
#   unlist() %>%
#   reduce() %>%
#   GeneRegionTrack(name = "Known")

# ranges(gencode_track)$transcript <- rep("transcript", length(ranges(gencode_track)))

gr_txdb <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf"))
gencode_track <- GeneRegionTrack(gr_txdb, name = "GENCODE", collapseTranscripts = "longest", transcriptAnnotation = "gene")
gencode_track <- gencode_track[transcript(gencode_track) == tx_to_plot]

displayPars(gencode_track) <- list(
  stacking = "squish",
  background.panel = "transparent",
  fill = "#009E73",
  col = "#009E73",
  lwd = 0.3,
  col.line = "black",
  fontcolor.title = "black",
  background.title = "#d1861d",
  thinBoxFeature = c("utr3", "utr5")
)

# Peptide track
peptide_gr <- import("/scratch/nxu/astrocytes_lrs_UCSC/hg38/pep_output.bb")
peptide_track <- GeneRegionTrack(peptide_gr, group = "transcript", name = "Peptides")

displayPars(peptide_track) <- list(
  stacking = "squish",
  background.panel = "transparent",
  fill = "black",
  col = "black",
  lwd = 0.3,
  col.line = "black",
  showId = TRUE,
  transcriptAnnotation = "transcript",
  fontcolor.title = "black",
  background.title = "#ef6548"
)

# ORFanage track
orfanage_txdb <- "nextflow_results/orfanage/minlen/orfanage.gtf" %>%
  makeTxDbFromGFF(format = "gtf")
orfanage_track <- GeneRegionTrack(orfanage_txdb, name = "ORFanage")
orfanage_track <- orfanage_track[transcript(orfanage_track) == transcript_id]

displayPars(orfanage_track) <- list(
  stacking = "squish",
  fill = "#E69F00",
  col = "#E69F00",
  lwd = 0.3,
  col.line = "black",
  showId = TRUE,
  background.panel = "transparent",
  transcriptAnnotation = "transcript",
  fontcolor.title = "black",
  background.title = "#045a8d"
)

# Collaborator translatome track
collaborator_txdb <- "nextflow_results/translatome/supplemented_collaborator/filtered_output_fixed.gtf" %>% 
  makeTxDbFromGFF(format = "gtf")
collaborator_track <- GeneRegionTrack(collaborator_txdb, name = "Translatome")
collaborator_track <- collaborator_track[transcript(collaborator_track) == ORF_id]

displayPars(collaborator_track) <- list(
  stacking = "squish",
  fill = "#E69F00",
  col = "#E69F00",
  lwd = 0.3,
  col.line = "black",
  showId = TRUE,
  background.panel = "transparent",
  transcriptAnnotation = "transcript",
  fontcolor.title = "black",
  background.title = "#045a8d"
)

# Ribo-seq coverage track
dtrack <- DataTrack(range = "nextflow_results/align/riboseq/minlen/merged_riboseq.bw", genome = "hg38", chromosome = "chr1", name = "Ribo-seq", type = "h")

# Combine all tracks
chr <- as.character(seqnames(ranges(gencode_track)))[1]
leftmost <- min(start(ranges(gencode_track)))
rightmost <- max(end(ranges(gencode_track)))
extra <- (rightmost - leftmost) * 0.05

# pdf(glue("figures/figure_2/genome_track_{current_gene}.pdf"), width = 9, height = 3)
track_grob <- grid.grabExpr(
  plotTracks(
    list(
      gencode_track, orfanage_track, collaborator_track, dtrack, peptide_track
    ),
    chromosome = chr,
    from = leftmost - extra,
    to = rightmost + extra,
    sizes = c(
      1, 1, 1, 1, 1
    ),
    groupAnnotation = "group", 
    just.group = "right",
    fontsize = 3
  )
)
source("/scratch/nxu/astrocytes/scripts/figures/figure_2/figure_2.R")
p / wrap_elements(full = track_grob) +
  plot_layout(heights = c(1, 1, 1, 1), widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 9))
ggsave("figures/figure_2/figure_2.pdf", width = 180, height = 140, units = "mm")

gencode_track <- GeneRegionTrack(gr_txdb, name = "GENCODE", transcriptAnnotation = "transcript", chromosome = chr, start = leftmost - extra, end = rightmost + extra)
gencode_track <- gencode_track[transcript(gencode_track) == tx_to_plot]
displayPars(gencode_track) <- list(
  thinBoxFeature = c("utr5", "utr3")
)

plotTracks(gencode_track, chromosome = chr, from = leftmost - extra, to = rightmost + extra)
