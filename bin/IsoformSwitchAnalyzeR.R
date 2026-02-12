#!/usr/bin/env Rscript
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(arrow)
library(rtracklayer)
library(readr)
library(argparse)

parser <- ArgumentParser(description='Create IsoformSwitchAnalyzeR object')
parser$add_argument('--final_expression', type='character', required=TRUE, help='Path to final expression parquet file')
parser$add_argument('--primer_to_sample', type='character', required=TRUE, help='Path to primer to sample CSV file')
parser$add_argument('--corrected_fasta', type='character', required=TRUE, help='Path to corrected FASTA file')
parser$add_argument('--orfanage_gtf', type='character', required=TRUE, help='Path to orfanage GTF file')
parser$add_argument('--annotation_gtf', type='character', required=TRUE, help='Path to annotation GTF file')
parser$add_argument('--final_classification', type='character', required=TRUE, help='Path to final classification parquet file')
args <- parser$parse_args()

final_expression <- read_parquet(args$final_expression)

primer_to_sample <- read_csv(args$primer_to_sample) %>% 
    rename(index = `Index primer`, sampleID = `Sample_ID`)

indices <- match(colnames(final_expression), primer_to_sample$index)
colnames(final_expression) <- primer_to_sample$sampleID[indices]
myDesign <- primer_to_sample %>% select(-index)

final_expression <- final_expression %>%
    rename(isoform_id = 1)

orfanage_tx <- rtracklayer::import(args$orfanage_gtf) %>%
    as_tibble() %>%
    filter(
        type == 'transcript'
    ) %>% 
    pull(transcript_id)

final_expression <- final_expression %>%
    filter(
        isoform_id %in% orfanage_tx
    )

aSwitchList <- importRdata(
    isoformCountMatrix = final_expression,
    designMatrix = myDesign,
    isoformExonAnnoation = args$orfanage_gtf,
    isoformNtFasta = args$corrected_fasta,
    addAnnotatedORFs = TRUE,
    fixStringTieAnnotationProblem = FALSE
)

aSwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE    
)

IsoseqsSwitchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchListFiltered,
    reduceToSwitchingGenes = FALSE,
    showProgress = TRUE
)

gene_id_mapping <- rtracklayer::import(args$annotation_gtf) %>% 
  as_tibble() %>%
  filter(type == 'gene') %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  dplyr::rename(associated_gene = gene_id)

classification <- read_parquet(args$final_classification) %>% 
    left_join(
        gene_id_mapping,
        by = 'associated_gene'
    )

IsoseqsSwitchList$isoformFeatures <- IsoseqsSwitchList$isoformFeatures %>% 
    left_join(
        dplyr::rename(dplyr::select(classification, c(isoform, gene_name)), isoform_id = isoform),
        by = 'isoform_id',
    ) %>% 
    select(-gene_name.x) %>% 
    dplyr::rename(gene_name=gene_name.y)

saveRDS(IsoseqsSwitchList, "IsoformSwitchAnalyzeR.rds")