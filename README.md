# Data preparation
## Short-read RNA-seq
- library prepared using takara stranded mRNA kit and sequenced on Novaseq X platform.

# Preprocessing
- Install SQANTI3
```bash
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.5.4/SQANTI3_v5.5.4.zip
mkdir sqanti3
unzip SQANTI3_v5.5.4.zip -d sqanti3
```

# Download data
```bash
# TransCODE phase 1 data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01369-0/MediaObjects/41587_2022_1369_MOESM2_ESM.xlsx -o data/41587_2022_1369_MOESM2_ESM.xlsx

# TransCODE phase 2 data
wget https://www.biorxiv.org/content/biorxiv/early/2025/07/07/2025.07.03.662928/DC1/embed/media-1.xlsx -o data/media-1.xlsx
```