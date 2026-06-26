library(ggplot2)
library(dplyr)

df <- read.delim(
  "nextflow_results/quality/minlen_quality_metrics.tsv",
  stringsAsFactors = FALSE
) |>
  filter(in_gencode_ribotie %in% c("True", "False")) |>
  mutate(
    ORBLv = as.numeric(ORBLv),
    in_gencode_ribotie = factor(
      in_gencode_ribotie,
      levels = c("True", "False"),
      labels = c("In GENCODE RiboTIE", "Not in GENCODE RiboTIE")
    )
  ) |>
  filter(!is.na(ORBLv))

counts <- df |> count(in_gencode_ribotie)
label_map <- setNames(
  sprintf("%s\n(n = %s)", counts$in_gencode_ribotie, format(counts$n, big.mark = ",")),
  as.character(counts$in_gencode_ribotie)
)
df$group_label <- label_map[as.character(df$in_gencode_ribotie)]
df$group_label <- factor(df$group_label, levels = label_map)

wt <- wilcox.test(
  ORBLv ~ in_gencode_ribotie,
  data = df,
  exact = FALSE
)
p_str <- if (wt$p.value < 1e-300) "p < 1e-300" else sprintf("p = %.2e", wt$p.value)

medians <- df |> group_by(in_gencode_ribotie) |> summarise(med = median(ORBLv))

p <- ggplot(df, aes(x = ORBLv, colour = group_label, fill = group_label)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(
    data = medians,
    aes(xintercept = med, colour = factor(
      in_gencode_ribotie,
      levels = c("In GENCODE RiboTIE", "Not in GENCODE RiboTIE"),
      labels = label_map
    )),
    linetype = "dashed", linewidth = 0.6
  ) +
  scale_colour_manual(values = c("In GENCODE RiboTIE\n(n = 22,626)" = "#2196F3",
                                  "Not in GENCODE RiboTIE\n(n = 28,399)" = "#E53935")) +
  scale_fill_manual(values   = c("In GENCODE RiboTIE\n(n = 22,626)" = "#2196F3",
                                  "Not in GENCODE RiboTIE\n(n = 28,399)" = "#E53935")) +
  annotate("text", x = 0.98, y = Inf, label = paste0("Wilcoxon rank-sum\n", p_str),
           hjust = 1, vjust = 1.4, size = 3.2, colour = "grey30") +
  labs(
    x = "ORBLv score",
    y = "Density",
    title = "ORBLv score distribution by GENCODE RiboTIE membership",
    colour = NULL, fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

out <- "nextflow_results/quality/orblv_by_gencode_ribotie.png"
ggsave(out, p, width = 7, height = 4.5, dpi = 150)
cat("Saved", out, "\n")

df |>
  group_by(in_gencode_ribotie) |>
  summarise(n = n(), median = median(ORBLv), mean = mean(ORBLv)) |>
  print()
cat("Wilcoxon W =", wt$statistic, ",", p_str, "\n")
