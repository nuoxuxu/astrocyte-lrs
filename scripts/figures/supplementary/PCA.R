library(Rsubread)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Define function to make PCA plot
make_pca_plot <- function(dds, output_file, title) {
  # Filter low-count genes (keeps ~60-70% of genes typically)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  # Variance-stabilizing transform (blind = TRUE for QC/exploration)
  vsd <- vst(dds, blind = TRUE)

  # DESeq2 uses top 500 most variable genes by default
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pct_var <- round(100 * attr(pca_data, "percentVar"), 1)

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 4, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    theme_bw(base_size = 13) +
    ggtitle(title)

  ggsave(output_file, plot = p, width = 6, height = 5)
  p
}

# For Ribo-seq samples
counts <- readRDS("proc/featureCounts_riboseq.rds")

colnames(counts) <- gsub("_Unmapped\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))

metadata <- list(
  "merged_astro_A" = "Unstim",
  "merged_astro_B" = "Unstim",
  "merged_astro_C" = "Stim",
  "Astro_D" = "Stim",
  "merged_astro_E" = "Stim",
  "merged_astro_F" = "Unstim",
  "Astro_J" = "Stim",
  "Astro_H" = "Unstim",
  "merged_astro_I" = "Unstim",
  "merged_astro_J1" = "Unstim",
  "Astro_K" = "Stim",
  "merged_astro_L" = "Stim"
) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  rename("sample" = "name", "condition" = "value")

metadata$condition <- factor(metadata$condition, levels = c("Unstim", "Stim"))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$sample
metadata <- metadata["condition"]
metadata <- metadata[match(colnames(counts), rownames(metadata)), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

make_pca_plot(
  dds,
  "figures/figure_1/riboseq_pca_plot.pdf",
  "PCA of Ribo-seq samples"
)

# For RNA-seq samples
counts <- readRDS("proc/featureCounts_rnaseq.rds")
colnames(counts) <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))

stimulation_label <- read.csv("data/stimulation_label.csv") %>%
  filter(Sample!="Astro_J") %>% 
  mutate(Sample = gsub("1", "", Sample))

samplessheet <- read.csv("data/samplesheet.csv") %>% 
  mutate(sample = gsub("_RNA", "", sample))%>% 
  filter(type=="rnaseq") %>% 
  select(sample, fastq_1) %>% 
  mutate(fastq_1 = gsub("/project/rrg-shreejoy/nxu/astrocytes/short_read/", "", fastq_1)) %>% 
  mutate(fastq_1 = gsub("_R1_001.fastq.gz", "", fastq_1)) %>% 
  left_join(stimulation_label, by = c("sample" = "Sample")) %>% 
  rename("condition" = "Condition")

name_map <- setNames(samplessheet$fastq_1, samplessheet$sample)
counts <- counts %>%
  as_tibble() %>% 
  rename(any_of(name_map))

metadata <- samplessheet %>% as.data.frame()
rownames(metadata) <- metadata$sample
metadata <- metadata %>% select(condition)
metadata$condition <- factor(metadata$condition, levels = c("Unstim", "Stim"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

make_pca_plot(
  dds,
  "figures/figure_1/rnaseq_pca_plot.pdf",
  "PCA of RNA-seq samples"
)