library(ggplot2)
library(dplyr)
library(readr)

unstim <- read_csv("ribotie_res/isoseq_ORFanage/custom_Unstim_with_gencode_ids.csv")
stim <- read_csv("ribotie_res/isoseq_ORFanage/custom_Stim_with_gencode_ids.csv")

unstim_n_ORF_type <- unstim %>% 
    mutate(GENCODE = !is.na(transcript_id_right)) %>%
    group_by(GENCODE, ORF_type) %>%
    summarise(count = n()) %>% 
    mutate(
        group = "Unstim",
        percentage = (count / sum(count)) * 100
    )

stim_n_ORF_type <- stim %>% 
    mutate(GENCODE = !is.na(transcript_id_right)) %>%
    group_by(GENCODE, ORF_type) %>%
    summarise(count = n()) %>% 
    mutate(
        group = "Stim",
        percentage = (count / sum(count)) * 100
    )

combined_n_ORF_type <- bind_rows(unstim_n_ORF_type, stim_n_ORF_type)

combined_n_ORF_type %>% 
    ggplot(aes(x = ORF_type, y = percentage, fill = GENCODE)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
        name = "Annotated in GENCODE",
        labels = c("No", "Yes"),
        values = c("FALSE" = "#fc8d62", "TRUE" = "#66c2a5")
    ) +
    facet_wrap(~group) +
    labs(
        title = "Distribution of ORF Types in Unstimulated and Stimulated Conditions",
        x = "ORF Type according to ORFanage annotation",
        y = "Percentage of ORFs"
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    )

ggsave("figures/RiboTIE_ORF_type_distribution.png", width = 8, height = 6)

unstim_null_counts <- unstim %>%
    mutate(is_null = is.na(transcript_id_right)) %>%
    group_by(is_null) %>%
    summarise(count = n()) %>%
    mutate(
        condition = "Unstim",
        percentage = (count / sum(count)) * 100
    )

stim_null_counts <- stim %>%
    mutate(is_null = is.na(transcript_id_right)) %>%
    group_by(is_null) %>%
    summarise(count = n()) %>%
    mutate(
        condition = "Stim",
        percentage = (count / sum(count)) * 100
    )

combined_null_counts <- bind_rows(unstim_null_counts, stim_null_counts)

combined_null_counts %>% 
    ggplot(aes(x = is_null, y = percentage, fill = is_null)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3.5) +
    facet_wrap(~condition) +
    labs(
        title = "Percentage of ORFs annotated in GENCODE",
        x = "Condition",
        y = "Percentage of ORFs"
    ) +
    scale_fill_manual(
        name = "Annotated in GENCODE",
        labels = c("No", "Yes"),
        values = c("FALSE" = "#fc8d62", "TRUE" = "#66c2a5")
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)
    )
ggsave("figures/RiboTIE_ORF_annotation_counts.png", width = 8, height = 6)
