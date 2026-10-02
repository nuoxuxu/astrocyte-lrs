library(ggVennDiagram)
library(Biostrings)
library(arrow)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
# ------ INPUT FILES ------
# CSV files with protein_seq column
csv_file1 <- args[1]
csv_file2 <- args[2]
csv_file3 <- args[3]

# ------ READ DATA ------
# Read CSV files
df1 <- read.csv(csv_file1, stringsAsFactors = FALSE)
df2 <- read.csv(csv_file2, stringsAsFactors = FALSE)
df3 <- read.csv(csv_file3, stringsAsFactors = FALSE)

# Extract unique protein sequences from CSVs
protein_seq_1 <- unique(df1$protein_seq)
protein_seq_2 <- unique(df2$protein_seq)
protein_seq_3 <- unique(df3$protein_seq)

# ------ CREATE VENN DIAGRAM ------
# Create list of protein sets
protein_list <- list(
  "gencode" = protein_seq_1,
  "high" = protein_seq_2,
  "low" = protein_seq_3
)

# Generate Venn diagram
venn_plot <- ggVennDiagram(protein_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "right") +
  labs(title = "Protein Sequence Overlap")

ggsave(args[4], venn_plot, width = 8, height = 6)