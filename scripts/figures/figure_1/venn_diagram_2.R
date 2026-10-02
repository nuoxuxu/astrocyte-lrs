library(ggVennDiagram)
library(Biostrings)
library(arrow)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
# ------ INPUT FILES ------
# CSV files with protein_seq column
csv_file1 <- args[1]
csv_file2 <- args[2]

# FASTA files
fasta_file <- file.path(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE", "gencode.v47.pc_translations.fa")
orfanage_fasta_file <- args[3]

# ------ READ DATA ------
# Read CSV files
df1 <- read.csv(csv_file1, stringsAsFactors = FALSE)
df2 <- read.csv(csv_file2, stringsAsFactors = FALSE)

# Extract unique protein sequences from CSVs
protein_seq_1 <- unique(df1$protein_seq)
protein_seq_2 <- unique(df2$protein_seq)

# Read ORFanage FASTA file
orfanage_seqs <- readAAStringSet(orfanage_fasta_file)
proteins_orfanage <- unique(as.character(orfanage_seqs))

# ------ CREATE VENN DIAGRAM ------
# Create list of protein sets
protein_list <- list(
  "unstim" = protein_seq_1,
  "stim" = protein_seq_2
)

# Generate Venn diagram
venn_plot <- ggVennDiagram(protein_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "right") +
  labs(title = "Protein Sequence Overlap")

ggsave(args[4], venn_plot, width = 8, height = 6)