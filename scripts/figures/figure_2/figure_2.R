library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(arrow)
library(magick)
library(figpatch)

ribotie_cpm1_3sample <- read_csv("nextflow_results/quality/collaborator/ribotie_cpm1_3sample_with_lncRNA.csv")

orf_colors <- c(
  "uORF"      = "#117733",  # dark green
  "uoORF"     = "#3B3181",  # dark purple
  "intORF"    = "#B5408A",  # magenta/pink
  "doORF"     = "#2AADA8",  # teal
  "dORF"      = "#8B8B2A",  # olive yellow
  "lncRNA ORF" = "#A63228"  # dark red/brick
)

# Add proportion validated by peptide evidence from mass spec
peptide_mapping <- read_parquet("nextflow_results/proteomic/peptide_mapping.parquet")

ncORF_list <- ribotie_cpm1_3sample %>%
    filter(!ORF_type %in% c(
        "annotated CDS",
        "N-terminal truncation",
        "N-terminal extension"
    )) %>%
    distinct(ORF_id) %>%
    pull(ORF_id)

novel_peptides <- read_csv(
    "nextflow_results/proteomic/novel_peptides.csv",
    show_col_types = FALSE
) %>%
    filter(novel_peptide)

validated_ncORFs <- peptide_mapping %>%
    filter(
        transcript_id %in% ncORF_list,
    ) %>%
    pull(transcript_id)

df_prop <- ribotie_cpm1_3sample %>%
    filter(!ORF_type %in% c("annotated CDS", "N-terminal extension", "N-terminal truncation")) %>%
    mutate(validated = ORF_id %in% validated_ncORFs) %>%
    group_by(ORF_type) %>%
    summarise(
        n = n(),
        n_validated = sum(validated),
        prop_validated = n_validated / n,
        .groups = "drop"
    )

ncorf_type_labels <- setNames(
  paste0(df_prop$ORF_type, " (n = ", df_prop$n, ")"),
  df_prop$ORF_type
)

orf_type_prop <- df_prop %>%
    ggplot(aes(x = "", y = n, fill = ORF_type)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) +
      coord_polar(theta = "y") +
      scale_fill_manual(
        name = "ncORF type",
        values = orf_colors,
        labels = ncorf_type_labels
      ) +
      labs(x = NULL, y = NULL) +
      scale_y_continuous(expand = c(0, 0)) + 
      theme_void() +
      theme(
        text = element_text(size = 7),
        legend.position = "right",
        legend.key.size  = unit(4, "mm"), 
        legend.margin = margin(0, 0, 0, 0)
      )

orf_type_validation_prop <- df_prop %>%
    ggplot(aes(x = ORF_type, fill = ORF_type)) +
    geom_col(aes(y = 1), width = 0.8, alpha = 0.25) +
    geom_col(aes(y = prop_validated), width = 0.8) +
    scale_fill_manual(
      name = "ncORF type",
      values = orf_colors
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "% peptide validation") +
    theme_classic() +
    theme(
      text = element_text(size = 7),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

# Do this at the isoform level
classification <- read_parquet("nextflow_results/sqanti3/isoseq/sqanti3_filter/final_classification.parquet")
ribotie_cpm1_3sample <- ribotie_cpm1_3sample %>% 
  left_join(classification %>% select("isoform", "structural_category"), by = c("transcript_id" = "isoform"))

validated_ORFs <- peptide_mapping %>% 
  group_by(transcript_id) %>%
  summarise(n_peptides = n_distinct(pep), .groups = "drop") %>% 
  filter(n_peptides > 1) %>% 
  pull(transcript_id)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC",
    "genic"                    = "Genic",
    "fusion"                   = "Fusion"
)

df_prop_isoform <- ribotie_cpm1_3sample %>%
    filter(ORF_type %in% c("annotated CDS", "N-terminal extension", "N-terminal truncation")) %>%
    mutate(validated = ORF_id %in% validated_ORFs) %>%
    group_by(structural_category) %>%
    summarise(
        n = n(),
        n_validated = sum(validated),
        prop_validated = n_validated / n,
        .groups = "drop"
    ) %>%
    mutate(
      structural_category = structural_category_labels[structural_category]
    )

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00",
    "Genic" = "#000000",
    "Fusion" = "#CC79A7"
)

structural_category_labels_with_n <- setNames(
  paste0(df_prop_isoform$structural_category, " (n = ", df_prop_isoform$n, ")"),
  df_prop_isoform$structural_category
)

orf_type_prop_isoform <- df_prop_isoform %>%
    ggplot(aes(x = "", y = n, fill = structural_category)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) +
      coord_polar(theta = "y") +
      scale_fill_manual(
        name = "Structural category",
        values = colorVector,
        labels = structural_category_labels_with_n
      ) +
      labs(x = NULL, y = NULL) +
      theme_void() +
      theme(
        text = element_text(size = 7),
        legend.key.size = unit(4, "mm"),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0)
      )

orf_type_validation_prop_isoform <- df_prop_isoform %>%
    ggplot(aes(x = structural_category, fill = structural_category)) +
    geom_col(aes(y = 1), width = 0.8, alpha = 0.25) +
    geom_col(aes(y = prop_validated), width = 0.8) +
    scale_fill_manual(
      name = "Structural category",
      values = colorVector
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "% peptide validation") +
    theme_classic() +
    theme(
      text = element_text(size = 7),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

isoform_schematic <- image_ggplot(image_read_pdf("/scratch/nxu/astrocytes/figures/figure_2/isoform_schematic.pdf"))
ncORF_schematic <- image_ggplot(image_read_pdf("/scratch/nxu/astrocytes/figures/figure_2/ncORF_schematic.pdf"))
pipeline_schematic <- image_ggplot(image_read_pdf("/scratch/nxu/astrocytes/figures/figure_2/pipeline_schematic.pdf"))

p <- pipeline_schematic / (ncORF_schematic | isoform_schematic) / (orf_type_prop | orf_type_validation_prop | orf_type_prop_isoform | orf_type_validation_prop_isoform)
  # plot_layout(heights = c(1, 1, 1), widths = c(1, 1, 1)) +
  # plot_annotation(tag_levels = "A") & 
  # theme(plot.tag = element_text(size = 9))
# ggsave("figures/figure_2.pdf", width = 180, height = 120, units = "mm")

# orf_type_prop | orf_type_validation_prop
# codon_colors <- c(
#   "ACG" = "#573D91",  # dark purple
#   "ATG" = "#903B8E",  # orange
#   "CTG" = "#D7507B",  # pink/magenta
#   "GTG" = "#FA8458",  # yellow
#   "TTG" = "#FDCB49"   # light green
# )

# start_codon_prop <- ribotie_cpm1_3sample %>% 
#   filter(!ORF_type %in% c("annotated CDS", "N-terminal extension", "N-terminal truncation")) %>%
#   count(ORF_type, start_codon) %>%
#   group_by(ORF_type) %>%
#   mutate(prop = n / sum(n) * 100) %>%
#   ungroup() %>%
#   ggplot(aes(x = ORF_type, y = prop, fill = start_codon)) +
#   geom_col(width = 0.8) +
#   scale_fill_manual(name = "Start codon", values = codon_colors) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = NULL, y = "Percentage of ORFs") +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.title = element_text(size = 9),
#     legend.text  = element_text(size = 8),
#     aspect.ratio = 0.7
#   )

# orf_type_validation_count <- df_prop %>%
#     ggplot(aes(x = ORF_type, fill = ORF_type)) +
#     geom_col(aes(y = n), width = 0.8, alpha = 0.25) +
#     geom_col(aes(y = n_validated), width = 0.8) +
#     scale_fill_manual(values = orf_colors) +
#     scale_y_continuous(expand = c(0, 0)) +
#     labs(x = NULL, y = "Number of ncORFs") +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     )
