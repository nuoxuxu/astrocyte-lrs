library(fgsea)
library(msigdbr)
library(dplyr)
library(readr)
library(ggplot2)

go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    split(x = .$gene_symbol, f = .$gs_name)
DEG <- read.csv("from_collaborator/AstroDEGs_stim_vs_unstim.txt", sep="\t") %>%
    filter(!is.na(gene_name))

ranked_genes <- setNames(-log10(DEG$padj) * sign(DEG$log2FoldChange), DEG$gene_name)
ranked_genes <- ranked_genes[is.finite(ranked_genes)]
ranked_genes

fgsea_res <- fgsea(
    pathways  = go_bp,
    stats     = ranked_genes,
    minSize   = 15,
    maxSize   = 500,
    nPermSimple = 10000
)

top_pathways <- fgsea_res %>%
    mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = 10) %>%
    ungroup() %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway), pathway = stringr::str_wrap(pathway, 40))

ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
    labs(x = NULL, y = "Normalized Enrichment Score", title = "GSEA — Hallmark gene sets") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

ggsave("figures/figure_1/gsea_barplot.pdf", width = 7, height = 6)