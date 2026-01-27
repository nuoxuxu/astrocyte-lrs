library(IsoformSwitchAnalyzeR)
library(dplyr)
library(arrow)
library(rtracklayer)
library(readr)

final_expression <- read_parquet("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_expression.parquet")

primer_to_sample <- read_csv("data/primer_to_sample.csv") %>% 
    rename(index = `Index primer`, sampleID = `Sample_ID`)

indices <- match(colnames(final_expression), primer_to_sample$index)
colnames(final_expression) <- primer_to_sample$sampleID[indices]
myDesign <- primer_to_sample %>% select(-index)

final_expression <- final_expression %>%
    rename(isoform_id = 1)

orfanage_tx <- rtracklayer::import("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/orfanage.gtf") %>%
    as_tibble() %>%
    filter(
        type == "transcript"
    ) %>% 
    pull(transcript_id)

final_expression <- final_expression %>%
    filter(
        isoform_id %in% orfanage_tx
    )

aSwitchList <- importRdata(
    isoformCountMatrix = final_expression,
    designMatrix = myDesign,
    isoformExonAnnoation = "nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/orfanage.gtf",
    isoformNtFasta = "nextflow_results/sqanti3/isoseq/sqanti3_qc/sqanti_qc_results_corrected.fasta",
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

gene_id_mapping <- rtracklayer::import("/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf") %>% 
  as_tibble() %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  dplyr::rename(associated_gene = gene_id)

classification <- read_parquet("nextflow_results/sqanti3/isoseq/sqanti3_filter/filter_by_expression/final_classification.parquet") %>% 
    left_join(
        gene_id_mapping,
        by = "associated_gene"
    )

IsoseqsSwitchList$isoformFeatures <- IsoseqsSwitchList$isoformFeatures %>% 
    left_join(
        dplyr::rename(dplyr::select(classification, c(isoform, gene_name)), isoform_id = isoform),
        by = "isoform_id",
    ) %>% 
    select(-gene_name.x) %>% 
    dplyr::rename(gene_name=gene_name.y)

IsoseqsSwitchList$isoformFeatures %>% 
    ggplot(aes(x=dIF, y=-log10(isoform_switch_q_value))) + 
    geom_point(
        aes( color=abs(dIF) > 0.5 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    ggrepel::geom_text_repel(data = IsoseqsSwitchList$isoformFeatures %>% filter(abs(dIF) > 0.5 & isoform_switch_q_value < 0.05), aes(label=gene_name),size=3, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
ggsave("figures/Isoform_Switching_DexSeq_Isoseqs.png", width=8, height=6)

AstroDEGs <- read_tsv("from_bharti/AstroDEGs_stim_vs_unstim.txt")

AstroDEGs %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(
        aes( color=abs(log2FoldChange) > 1 & padj < 0.05 ), # custom cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') +
    geom_vline(xintercept = c(-1, 1), linetype='dashed') +
    scale_color_manual('Signficant\nDEG', values = c('black','red')) +
    labs(x='Log2 Fold Change', y='-Log10 ( Adjusted P Value )') +
    coord_cartesian(ylim=c(0, 10), xlim=c(-6,6)) +
    theme_bw()
ggsave("figures/Astrocyte_DEGs_stim_vs_unstim.png", width=8, height=6)