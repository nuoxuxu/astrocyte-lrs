library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(rtracklayer)

# Function to read GTF and extract transcript IDs
read_gtf_transcripts <- function(path, feature_filter = NULL) {
    gtf <- rtracklayer::import(path) %>%
        as.data.frame()

    if (!is.null(feature_filter)) {
        gtf <- gtf %>% filter(type == feature_filter)
    }

    gtf %>%
        pull(transcript_id) %>%
        unique()
}

# Function to compute percentages for a given param set
get_ribotie_percent <- function(param_set_name) {
    stim <- read_csv(
        paste0("nextflow_results/ribotie/", param_set_name, "/ribotie_res_Stim.csv"),
        show_col_types = FALSE
    )

    orfanage_transcripts <- if (param_set_name %in% c("high_stringency", "low_stringency")) {
        read_gtf_transcripts(
            paste0("nextflow_results/orfanage/", param_set_name, "/orfanage.gtf")
        )
    } else {
        read_gtf_transcripts("data/gencode.v47.annotation.gtf", feature_filter = "CDS")
    }

    n_total <- length(orfanage_transcripts)
    n_in_ribotie <- stim %>%
        pull(transcript_id) %>%
        unique() %>%
        length()
    n_not_in_ribotie <- n_total - n_in_ribotie

    tibble(
        param_set  = param_set_name,
        category   = c("In RiboTIE", "Not in RiboTIE"),
        count      = c(n_in_ribotie, n_not_in_ribotie),
        total      = n_total,
        percentage = count / total * 100
    )
}

# Collect data for all param sets
param_sets <- c("gencode", "high_stringency", "low_stringency")

plot_data <- map_dfr(param_sets, get_ribotie_percent) %>%
    mutate(
        param_set = factor(param_set, levels = param_sets),
        category  = factor(category, levels = c("Not in RiboTIE", "In RiboTIE")) # stack order
    )

# Print summary
plot_data %>%
    filter(category == "In RiboTIE") %>%
    select(param_set, percentage) %>%
    mutate(label = sprintf("%.1f%%", percentage)) %>%
    print()

# Plot
ggplot(plot_data, aes(x = param_set, y = count, fill = category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(
        data = plot_data %>% filter(category == "In RiboTIE"),
        aes(label = sprintf("%.1f%%", percentage)),
        position = position_stack(vjust = 0.5),
        color = "white", fontface = "bold", size = 4.5
    ) +
    scale_fill_manual(
        values = c("In RiboTIE" = "#2C7BB6", "Not in RiboTIE" = "#D7191C"),
        name   = NULL
    ) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
        title    = "ORFanage ORFs predicted by RiboTIE",
        subtitle = "Bar height = total ORF count; label = % predicted",
        x        = "Parameter Set",
        y        = "Number of ORFs"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        legend.position    = "top",
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold"),
        axis.text.x        = element_text(size = 11)
    )

ggsave("figures/percent_ORFS_predicted_by_RiboTIE.png")
