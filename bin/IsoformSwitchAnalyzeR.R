#!/usr/bin/env Rscript
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(arrow)
library(rtracklayer)
library(readr)
library(argparse)

parser <- ArgumentParser(description='Create IsoformSwitchAnalyzeR object')
parser$add_argument('--final_expression',     type='character', required=TRUE, help='Path to final expression parquet file')
parser$add_argument('--primer_to_sample',     type='character', required=TRUE, help='Path to primer to sample CSV file')
parser$add_argument('--final_fasta',          type='character', required=TRUE, help='Path to corrected FASTA file')
parser$add_argument('--orfanage_gtf',         type='character', required=TRUE, help='Path to orfanage GTF file')
parser$add_argument('--annotation_gtf',       type='character', required=TRUE, help='Path to annotation GTF file')
parser$add_argument('--final_classification', type='character', required=TRUE, help='Path to final classification parquet file')
parser$add_argument('--cpat_results',         type='character', required=TRUE, help='Path to filtered CPAT results file')
parser$add_argument('--pfam_results',         type='character', required=TRUE, help='Path to filtered pfam results file')
args <- parser$parse_args()

# Rename final_expression column names from index primers to sample IDs
final_expression <- read_parquet(args$final_expression)
primer_to_sample <- read_csv(args$primer_to_sample) %>% 
    rename(index = `Index primer`, sampleID = `Sample_ID`)
indices <- match(colnames(final_expression), primer_to_sample$index)
colnames(final_expression) <- primer_to_sample$sampleID[indices]
final_expression <- final_expression %>%
    rename(isoform_id = 1)

# Make sure "Unstim" is the reference level for condition in the design matrix
myDesign <- primer_to_sample %>%
    select(-index) %>%
    mutate(condition = relevel(factor(condition), ref = "Unstim"))

# Make sure that final_expression and orfanage_gtf have the same set of transcripts
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

# Get Ensembl ID to gene name mapping from annotation GTF
gene_id_mapping <- rtracklayer::import(args$annotation_gtf) %>%
  as_tibble() %>%
  filter(type == 'gene') %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  dplyr::rename(associated_gene = gene_id)

id_to_name <- setNames(gene_id_mapping$gene_name, gene_id_mapping$associated_gene)

resolve_gene_name <- function(compound_id) {
    parts <- strsplit(compound_id, "_")[[1]]
    names <- id_to_name[parts]
    names[is.na(names)] <- parts[is.na(names)]
    paste(names, collapse = "_")
}

# Creat GRanges object for isoformExonAnnoation
isoformExonAnnoation <- rtracklayer::import(args$orfanage_gtf)
isoformExonAnnoation <- subset(isoformExonAnnoation, mcols(isoformExonAnnoation)$type=="exon")
mcols(isoformExonAnnoation)$type <- NULL

classification <- read_parquet(args$final_classification) %>%
    mutate(gene_name = sapply(associated_gene, resolve_gene_name))

mcols(isoformExonAnnoation)$gene_name <- classification$gene_name[match(
        sub("_\\d+$", "", mcols(isoformExonAnnoation)$transcript_id),
        classification$isoform
    )]

mcols(isoformExonAnnoation)$gene_id <- mcols(isoformExonAnnoation)$gene_name
names(mcols(isoformExonAnnoation))[names(mcols(isoformExonAnnoation)) == "transcript_id"] <- "isoform_id"

# Create switchAnalyzeRlist object
IsoseqsSwitchList <- importRdata(
    isoformCountMatrix = final_expression,
    designMatrix = myDesign,
    isoformExonAnnoation = isoformExonAnnoation,
    isoformNtFasta = args$final_fasta,
    addAnnotatedORFs = TRUE,
    fixStringTieAnnotationProblem = FALSE
)

IsoseqsSwitchList <- addORFfromGTF(IsoseqsSwitchList, args$orfanage_gtf)

# Pre-filtering
# IsoseqsSwitchList <- preFilter(
#     switchAnalyzeRlist = IsoseqsSwitchList,
#     geneExpressionCutoff = 0,
#     isoformExpressionCutoff = 0,
#     removeSingleIsoformGenes = TRUE
# )

IsoseqsSwitchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = IsoseqsSwitchList,
    reduceToSwitchingGenes = FALSE,
    showProgress = TRUE
)

# Gene-level differential expression via DESeq2
geneCountMatrix <- extractGeneExpression(
    IsoseqsSwitchList,
    extractCounts = TRUE
)

geneDDS <- DESeqDataSetFromMatrix(
    countData = round(geneCountMatrix[, -(1:2)]),
    colData = IsoseqsSwitchList$designMatrix,
    design = ~ condition
)
geneDDS$condition <- relevel(geneDDS$condition, ref = "Unstim")
geneDDS <- DESeq(geneDDS)
geneRes <- results(geneDDS)
geneDE <- data.frame(
    gene_id = geneCountMatrix[, 1],
    gene_q_value = geneRes$padj
)
geneDE$gene_q_value[is.na(geneDE$gene_q_value)] <- 1

# Isoform-level differential expression via DESeq2
isoCountMatrix <- IsoseqsSwitchList$isoformCountMatrix
rownames(isoCountMatrix) <- isoCountMatrix$isoform_id
isoCountMatrix <- isoCountMatrix[, -1]

isoDDS <- DESeqDataSetFromMatrix(
    countData = round(isoCountMatrix),
    colData = IsoseqsSwitchList$designMatrix,
    design = ~ condition
)
isoDDS$condition <- relevel(isoDDS$condition, ref = "Unstim")
isoDDS <- DESeq(isoDDS)
isoRes <- results(isoDDS)
isoDE <- data.frame(
    isoform_id = rownames(isoCountMatrix),
    iso_q_value = isoRes$padj
)
isoDE$iso_q_value[is.na(isoDE$iso_q_value)] <- 1

# Import DE results into switchAnalyzeRlist
idx_gene <- match(
    IsoseqsSwitchList$isoformFeatures$gene_id,
    geneDE$gene_id
)
IsoseqsSwitchList$isoformFeatures$gene_q_value <- geneDE$gene_q_value[idx_gene]

idx_iso <- match(
    IsoseqsSwitchList$isoformFeatures$isoform_id,
    isoDE$isoform_id
)
IsoseqsSwitchList$isoformFeatures$iso_q_value <- isoDE$iso_q_value[idx_iso]

# Analyze switch consequences
IsoseqsSwitchList <- analyzeIntronRetention(IsoseqsSwitchList)
IsoseqsSwitchList <- analyzeCPAT(
    switchAnalyzeRlist   = IsoseqsSwitchList,
    pathToCPATresultFile = args$cpat_results,
    codingCutoff         = 0.725,
    removeNoncodinORFs   = TRUE
)
IsoseqsSwitchList <- analyzeAlternativeSplicing(IsoseqsSwitchList)

IsoseqsSwitchList <- extractSequence(IsoseqsSwitchList)
IsoseqsSwitchList <- analyzeSwitchConsequences(
    switchAnalyzeRlist = IsoseqsSwitchList, 
    consequencesToAnalyze = c('intron_retention', 'coding_potential', 'ORF_seq_similarity', 'NMD_status')
)

IsoseqsSwitchList <- analyzePFAM(
    switchAnalyzeRlist = IsoseqsSwitchList,
    pathToPFAMresultFile = args$pfam_results
)

saveRDS(IsoseqsSwitchList, "IsoformSwitchAnalyzeR.rds")

IsoseqsSwitchList$isoformFeatures %>% write.csv("isoformFeatures.csv", row.names = FALSE)