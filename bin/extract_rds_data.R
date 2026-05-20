#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
    make_option("--rds",    type = "character", help = "Path to IsoformSwitchAnalyzeR .rds file"),
    make_option("--output", type = "character", default = "switchPlotFromTables.RData",
                help = "Output .RData file path [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rds)) stop("--rds is required")

message("Reading RDS: ", opt$rds)
obj <- readRDS(opt$rds)

message("Extracting tables...")
isoformFeatures <- obj$isoformFeatures
conditions      <- obj$conditions
exons           <- obj$exons
orfAnalysis     <- obj$orfAnalysis
domainAnalysis  <- obj$domainAnalysis

message("Saving to: ", opt$output)
save(isoformFeatures, conditions, exons, orfAnalysis, domainAnalysis,
     file = opt$output)

message("Done. Objects saved: isoformFeatures, conditions, exons, orfAnalysis, domainAnalysis")
message("  isoformFeatures: ", nrow(isoformFeatures), " rows")
message("  conditions:      ", nrow(conditions), " rows")
message("  exons:           ", length(exons), " ranges")
message("  orfAnalysis:     ", if (is.null(orfAnalysis)) "NULL" else paste(nrow(orfAnalysis), "rows"))
message("  domainAnalysis:  ", if (is.null(domainAnalysis)) "NULL" else paste(nrow(domainAnalysis), "rows"))
