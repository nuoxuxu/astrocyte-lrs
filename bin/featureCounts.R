#!/usr/bin/env Rscript
library(Rsubread)
library(dplyr)
library(readr)
library(argparse)

parser <- ArgumentParser(description='Count reads aligned to genomic features using Rsubread')
parser$add_argument('--annotation-gtf', type='character', required=TRUE, help='Path to annotation GTF file')
parser$add_argument('--output-csv', type='character', required=TRUE, help='Output CSV file path')
parser$add_argument('--is-paired-end', type='logical', required=TRUE, help='TRUE for paired-end, FALSE for single-end')
parser$add_argument('--feature-type', type='character', default='exon', help='GTF feature type to count (default: exon)')
parser$add_argument('--attr-type', type='character', default='gene_id', help='GTF attribute type for grouping (default: gene_id)')
parser$add_argument('--nthreads', type='integer', default=4, help='Number of threads to use (default: 4)')
args <- parser$parse_args()

# Find all genome-mapped BAM files in current directory (staged by Nextflow)
# Only include .Aligned.sortedByCoord.out.bam files (STAR genome alignment output)
bam_files <- list.files(".", pattern="Aligned\\.sortedByCoord\\.out\\.bam$", full.names=TRUE)
if (length(bam_files) == 0) {
  stop("No genome-mapped BAM files found (.Aligned.sortedByCoord.out.bam)")
}
annotation_gtf <- normalizePath(args$`annotation_gtf`, mustWork=TRUE)
output_csv <- normalizePath(args$`output_csv`, mustWork=FALSE)

cat("BAM files to process:\n")
for (f in bam_files) cat("  ", f, "\n")
cat("Annotation GTF:", annotation_gtf, "\n")
cat("Output CSV:", output_csv, "\n\n")

# featureCounts — fastest approach
fc <- featureCounts(
  files      = bam_files,
  annot.ext  = annotation_gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = args$`feature_type`,
  GTF.attrType        = args$`attr_type`,
  isPairedEnd  = args$`is_paired_end`,
  nthreads     = args$`nthreads`
)

counts <- fc$counts
counts %>% write_csv(output_csv)

cat("\nfeatureCounts completed. Output written to:", output_csv, "\n")
