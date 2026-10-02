#!/usr/bin/env Rscript
# Protein differential expression: Stimulated (+) vs Unstimulated (-) astrocytes.
#
# Input is the DIA-NN protein-group matrix from the FragPipe DIA workflow
# (timsTOF diaPASEF). DIA-NN 1.8 reports MaxLFQ intensities that are already
# cross-run normalised, so no extra normalisation is applied by default.
#
# Sample columns are raw file paths such as
#   .../ME20260612_BN_Kalish_Astrocytes_8_F+_600ng_S3-A8_1_10440.d
# where "F" is the culture/donor and "+"/"-" is the condition. C, E and F are
# measured in both conditions, so the culture is used as a blocking factor
# (limma duplicateCorrelation) unless --no-block is given.
#
# Steps:
#   1. log2-transform intensities (0 -> NA)
#   2. keep protein groups with >= --min-valid values in at least one condition
#   3. impute remaining NAs Perseus-style (left-censored normal: downshift 1.8 SD,
#      width 0.3 SD per sample) — appropriate for timsTOF DIA where missingness is
#      mostly low-abundance (MNAR)
#   4. limma (trend + robust eBayes) with contrast Stim - Unstim
#
# Called by the protein_differential_expression process in post_RiboTIE.nf.
# Usage:
#   proteomics_differential_expression.R \
#     [--input nextflow_results/proteomics/fragpipe/dia-quant-output/report.pg_matrix.tsv] \
#     [--outdir nextflow_results/proteomics/differential_expression] \
#     [--figdir figures/proteomics] [--min-valid 2] [--no-impute] \
#     [--median-normalize] [--no-block]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option("--input", default = "nextflow_results/proteomics/fragpipe/dia-quant-output/report.pg_matrix.tsv"),
  make_option("--outdir", default = "nextflow_results/proteomics/differential_expression"),
  make_option("--figdir", default = "figures/proteomics"),
  make_option("--min-valid", type = "integer", default = 2,
              help = "Min non-missing values required in at least one condition [%default]"),
  make_option("--no-impute", action = "store_true", default = FALSE,
              help = "Skip imputation; limma is run on proteins with missing values as-is"),
  make_option("--median-normalize", action = "store_true", default = FALSE,
              help = "Median-centre log2 intensities per sample (DIA-NN already normalises)"),
  make_option("--no-block", action = "store_true", default = FALSE,
              help = "Ignore culture pairing (A/C/E/F) and run an unpaired comparison"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 1,
              help = "|log2FC| threshold used for labelling significant proteins [%default]"),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load and parse sample metadata ------------------------------------------
pg <- fread(opt$input)
annot_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                "First.Protein.Description")
sample_cols <- setdiff(names(pg), annot_cols)

sample_re <- "Astrocytes_(\\d+)_([A-Za-z]+)([+-])_"
if (!all(grepl(sample_re, sample_cols))) {
  stop("Could not parse condition from columns:\n",
       paste(sample_cols[!grepl(sample_re, sample_cols)], collapse = "\n"))
}
m <- regmatches(basename(sample_cols), regexec(sample_re, basename(sample_cols)))
meta <- data.table(
  column    = sample_cols,
  run       = as.integer(sapply(m, `[`, 2)),
  culture   = sapply(m, `[`, 3),
  sign      = sapply(m, `[`, 4)
)
meta[, condition := factor(ifelse(sign == "+", "Stim", "Unstim"), levels = c("Unstim", "Stim"))]
meta[, sample := paste0(culture, sign)]
setorder(meta, condition, culture)
print(meta[, .(sample, run, culture, condition)])

# ---- Log2 matrix --------------------------------------------------------------
mat <- as.matrix(pg[, meta$column, with = FALSE])
mat[mat <= 0] <- NA
mat <- log2(mat)
colnames(mat) <- meta$sample
rownames(mat) <- pg$Protein.Group

if (opt$`median-normalize`) {
  mat <- sweep(mat, 2, apply(mat, 2, median, na.rm = TRUE)) + median(mat, na.rm = TRUE)
}

# ---- Valid-value filter -------------------------------------------------------
n_valid <- sapply(levels(meta$condition), function(g) {
  rowSums(!is.na(mat[, meta$condition == g, drop = FALSE]))
})
keep <- apply(n_valid, 1, max) >= opt$`min-valid`
message(sprintf("Protein groups: %d total, %d pass valid-value filter (>= %d in one condition)",
                nrow(mat), sum(keep), opt$`min-valid`))
mat <- mat[keep, , drop = FALSE]
n_valid <- n_valid[keep, , drop = FALSE]
mat_raw <- mat

# ---- Imputation (Perseus-style left-censored) ---------------------------------
if (!opt$`no-impute`) {
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    na <- is.na(x)
    if (!any(na)) next
    mu <- mean(x, na.rm = TRUE) - 1.8 * sd(x, na.rm = TRUE)
    s  <- 0.3 * sd(x, na.rm = TRUE)
    mat[na, j] <- rnorm(sum(na), mu, s)
  }
}

# ---- limma --------------------------------------------------------------------
condition <- meta$condition
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
contr <- makeContrasts(Stim_vs_Unstim = Stim - Unstim, levels = design)

block <- NULL
corr <- NULL
if (!opt$`no-block` && any(duplicated(meta$culture))) {
  dc <- duplicateCorrelation(mat, design, block = meta$culture)
  corr <- dc$consensus.correlation
  block <- meta$culture
  message(sprintf("Blocking on culture; consensus within-culture correlation = %.3f", corr))
}

fit <- lmFit(mat, design, block = block, correlation = corr)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "Stim_vs_Unstim", number = Inf, sort.by = "none")

# ---- Assemble results ---------------------------------------------------------
res <- data.table(Protein.Group = rownames(tt))
res <- merge(res, pg[, ..annot_cols], by = "Protein.Group", sort = FALSE)
res[, `:=`(
  log2FC      = tt$logFC,
  AveExpr     = tt$AveExpr,
  t           = tt$t,
  P.Value     = tt$P.Value,
  adj.P.Val   = tt$adj.P.Val,
  n_valid_Unstim = n_valid[Protein.Group, "Unstim"],
  n_valid_Stim   = n_valid[Protein.Group, "Stim"]
)]
res[, n_imputed := if (opt$`no-impute`) 0L else
      (sum(meta$condition == "Unstim") - n_valid_Unstim) + (sum(meta$condition == "Stim") - n_valid_Stim)]
res[, direction := fifelse(adj.P.Val < opt$fdr & log2FC >=  opt$lfc, "Up in Stim",
                   fifelse(adj.P.Val < opt$fdr & log2FC <= -opt$lfc, "Down in Stim", "NS"))]
setorder(res, P.Value)

fwrite(res, file.path(opt$outdir, "protein_DE_stim_vs_unstim.tsv"), sep = "\t")
fwrite(data.table(Protein.Group = rownames(mat_raw), mat_raw),
       file.path(opt$outdir, "log2_intensity_filtered.tsv"), sep = "\t", na = "NA")
fwrite(data.table(Protein.Group = rownames(mat), mat),
       file.path(opt$outdir, "log2_intensity_imputed.tsv"), sep = "\t")
fwrite(meta[, .(sample, column, run, culture, condition)],
       file.path(opt$outdir, "sample_metadata.tsv"), sep = "\t")

message("Significant at FDR < ", opt$fdr, ", |log2FC| >= ", opt$lfc, ":")
print(res[, .N, by = direction])

# ---- Figures ------------------------------------------------------------------
pdf(file.path(opt$figdir, "protein_DE_stim_vs_unstim.pdf"), width = 7, height = 6)

# Intensity distributions (pre-imputation)
long <- melt(data.table(Protein.Group = rownames(mat_raw), mat_raw),
             id.vars = "Protein.Group", variable.name = "sample", value.name = "log2_int",
             na.rm = TRUE)
long <- merge(long, meta[, .(sample, condition)], by = "sample")
print(ggplot(long, aes(sample, log2_int, fill = condition)) +
        geom_boxplot(outlier.size = 0.3) +
        labs(title = "log2 MaxLFQ intensity (before imputation)", x = NULL, y = "log2 intensity") +
        theme_bw())

# PCA
pca <- prcomp(t(mat), scale. = FALSE)
ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc <- data.table(pca$x[, 1:2], sample = rownames(pca$x))
pc <- merge(pc, meta[, .(sample, culture, condition)], by = "sample")
print(ggplot(pc, aes(PC1, PC2, colour = condition, shape = culture)) +
        geom_point(size = 3) +
        geom_text_repel(aes(label = sample), show.legend = FALSE) +
        labs(title = "PCA of protein log2 intensities",
             x = sprintf("PC1 (%s%%)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2])) +
        theme_bw())

# Volcano
cols <- c("Up in Stim" = "#D55E00", "Down in Stim" = "#0072B2", "NS" = "grey70")
res[, label := fifelse(!is.na(Genes) & nzchar(as.character(Genes)), as.character(Genes), Protein.Group)]
top_lab <- res[direction != "NS"][order(P.Value)][1:min(.N, 20)]
print(ggplot(res, aes(log2FC, -log10(P.Value), colour = direction)) +
        geom_point(size = 0.8, alpha = 0.7) +
        geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = "dashed", colour = "grey40") +
        geom_text_repel(data = top_lab, aes(label = label), size = 2.5,
                        max.overlaps = 30, show.legend = FALSE) +
        scale_colour_manual(values = cols) +
        labs(title = "Stimulated vs Unstimulated",
             subtitle = sprintf("limma%s; FDR < %s, |log2FC| >= %s",
                                if (is.null(block)) "" else " (blocked on culture)",
                                opt$fdr, opt$lfc),
             x = "log2 fold change (Stim / Unstim)", y = "-log10 P", colour = NULL) +
        theme_bw())

# P-value histogram
print(ggplot(res, aes(P.Value)) + geom_histogram(bins = 40, boundary = 0) +
        labs(title = "P-value distribution") + theme_bw())

invisible(dev.off())
message("Results written to ", opt$outdir, " and ", opt$figdir)
