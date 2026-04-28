#!/usr/bin/env Rscript
library(IsoformSwitchAnalyzeR)
library(argparse)

parser <- ArgumentParser(description='Plot splicing and switch consequence enrichments')
parser$add_argument('--rds', type='character', required=TRUE, help='Path to IsoformSwitchAnalyzeR RDS file')
args <- parser$parse_args()

IsoseqsSwitchList <- readRDS(args$rds)

pdf(file='splicing_consequences.pdf', onefile=FALSE, height=6, width=9)
extractSplicingEnrichment(
    IsoseqsSwitchList,
    localTheme=theme_bw(base_size=14),
    returnResult=TRUE
)
dev.off()

pdf(file='switch_consequences.pdf', onefile=FALSE, height=6, width=9)
extractConsequenceEnrichment(
    IsoseqsSwitchList,
    analysisOppositeConsequence=TRUE,
    localTheme=theme_bw(base_size=14),
    returnResult=TRUE
)
dev.off()
