library(ggVennDiagram)
library(Biostrings)
library(arrow)
library(ggplot2)
# ------ INPUT FILES ------
# CSV files with protein_seq column
csv_file1 <- "ribotie_res/GRCh38v47/GRCh38v47_Unstim.csv"
csv_file2 <- "ribotie_res/isoseq_ORFanage/custom_Unstim.csv"

# FASTA files
fasta_file <- file.path(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE", "gencode.v47.pc_translations.fa")
orfanage_fasta_file <- "nextflow_results/orfanage/orfanage_proteins.fasta"

# ------ READ DATA ------
# Read CSV files
df1 <- read.csv(csv_file1, stringsAsFactors = FALSE)
df2 <- read.csv(csv_file2, stringsAsFactors = FALSE)

# Extract unique protein sequences from CSVs
GRCh38v47_unstim <- unique(df1$protein_seq)
custom_unstim <- unique(df2$protein_seq)

# Read FASTA file and extract sequences
fasta_seqs <- readAAStringSet(fasta_file)

# Extract transcript ID (second element separated by "|") as names
original_names <- names(fasta_seqs)
transcript_ids <- sapply(strsplit(original_names, "\\|"), `[`, 2)
names(fasta_seqs) <- transcript_ids

# Load parquet file and get associated transcripts to include
# classification_df <- read_parquet(parquet_file)
# transcripts_to_include <- unique(classification_df$associated_transcript)

proteins_fasta <- unique(as.character(fasta_seqs))

# Read ORFanage FASTA file
orfanage_seqs <- readAAStringSet(orfanage_fasta_file)
proteins_orfanage <- unique(as.character(orfanage_seqs))

# ------ CALCULATE OVERLAP PERCENTAGE ------
# Percentage of custom_unstim proteins found in proteins_fasta
overlap_count <- sum(custom_unstim %in% proteins_fasta)
overlap_pct <- (overlap_count / length(custom_unstim)) * 100
cat(sprintf("Proteins in custom_unstim found in GENCODE: %d / %d (%.2f%%)\n",
            overlap_count, length(custom_unstim), overlap_pct))

# ------ CREATE VENN DIAGRAM ------
# Create list of protein sets
protein_list <- list(
  "GRCh38v47_Unstim" = GRCh38v47_unstim,
  "GENCODE" = proteins_fasta
)

# Generate Venn diagram
venn_plot <- ggVennDiagram(protein_list,
                            label_alpha = 0,
                            category.names = c("GRCh38v47 + RiboTIE", "GENCODEv47")) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "right") +
  labs(title = "Protein Sequence Overlap")

# Display the plot
# print(venn_plot)

# Optionally save the plot
ggsave("figures/protein_venn_diagram.pdf", venn_plot, width = 8, height = 6)

# ------ CREATE 4-SET VENN DIAGRAM (with ORFanage) ------
protein_list_4 <- list(
  "custom_Unstim" = custom_unstim,
  "GRCh38v47_Unstim" = GRCh38v47_unstim,
  "GENCODE" = proteins_fasta,
  "ORFanage" = proteins_orfanage
)

venn_plot_4 <- ggVennDiagram(protein_list_4,
                              label_alpha = 0,
                              category.names = c("LR + RiboTIE", "GRCh38v47 + RiboTIE", "GENCODEv47", "ORFanage")) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "right") +
  labs(title = "Protein Sequence Overlap (with ORFanage)")

ggsave("figures/protein_venn_diagram_4set.pdf", venn_plot_4, width = 10, height = 8)
