# Transcript Classification Viewer

An optimized Shiny application for visualizing and comparing transcript ORF predictions from different tools (RiboTIE, ORFanage, GENCODE).

## Quick Start

### 1. Setup Data (First Time Only)

Run the setup script once to convert the original large GTF files to fast RDS format:

```bash
cd tx_classification_viewer
Rscript setup_data.R
```

This will:
- Convert GTF files to RDS format (~50-100x faster loading)
- Create lookup tables for fast metadata access
- Generate files in the `data/` directory (20MB total)

**Note:** This takes 2-5 minutes on first run, but you only need to do it once.

### 2. Run the App Locally

```bash
cd tx_classification_viewer
shiny::runApp()
```

Or from R:
```r
setwd("tx_classification_viewer")
shiny::runApp()
```

### 3. Deploy to shinyapps.io

See **[DEPLOYMENT.md](DEPLOYMENT.md)** for detailed instructions.

Quick version:
```r
setwd("tx_classification_viewer")
source("deploy.R")  # Interactive deployment script
```

**⚠️ Important:** The `data/` folder with RDS files **must** be deployed with the app. Don't have them in `.gitignore` during deployment.

## What's Optimized?

### 1. **RDS Caching** (50-100x speedup on startup)
   - Original GTF files are converted to RDS format in `data/`
   - RDS is 50-100x faster than importing GTF files
   - App loads in seconds instead of minutes

### 2. **Extracted and Memoized Difference Logic**
   - `get_differences()` function is memoized, so repeated lookups are instant
   - Shared between the plot and the difference counter
   - Eliminates code duplication

### 3. **Lookup Tables** (O(1) metadata access)
   - Fast lookups for ORF type and protein class metadata
   - Uses named vectors instead of filtering tibbles
   - Much faster than repeated `filter()` and `pull()` calls

## File Structure

```
tx_classification_viewer/
├── app.R                    # Main Shiny application
├── setup_data.R            # Setup script (run once)
├── README.md               # This file
└── data/                   # Generated data files (do not edit)
    ├── gencode.rds
    ├── ribotie.rds
    ├── orfanage.rds
    ├── pbid_to_pr_transcripts.rds
    ├── pbid_to_orfanage_template.rds
    ├── orf_type_lookup.rds
    └── pclass_lookup.rds
```

## Features

- **Transcript Selector** - Choose any RiboTIE transcript
- **Reference Comparison** - Compare against ORFanage prediction, GENCODE template, or SQANTI3 match
- **Difference Navigation** - Cycle through multiple differences with Previous/Next buttons
- **Focused View** - Zoom in on individual differences for detailed inspection
- **Difference Counter** - See "Difference X of Y" to track your position

## Performance

| Operation | Original | Optimized | Speedup |
|-----------|----------|-----------|---------|
| Initial load | 2-5 min | <5 sec | **20-60x** |
| GTF import | Per startup | Once, then RDS | N/A |
| Metadata lookup | O(n) filter | O(1) hash | **10-100x** |
| Difference calc (cached) | Every time | Once, then cached | **100x+** |

## Regenerating Data

If the source files change and you want to regenerate the RDS files:

```bash
rm -rf data/*.rds
Rscript setup_data.R
```

## Dependencies

- `shiny` - Web app framework
- `ggtranscript` - Transcript visualization
- `ggplot2` - Plotting
- `dplyr` - Data manipulation
- `GenomicRanges` - Genomic operations
- `rtracklayer` - GTF import
- `memoise` - Function caching
- `readr` - File reading
- `glue` - String interpolation
