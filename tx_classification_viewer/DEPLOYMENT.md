# Deployment Guide for tx_classification_viewer

## Local Development

```r
setwd("tx_classification_viewer")
shiny::runApp()
```

---

## Deployment to shinyapps.io

### Important: Data Files Must Be Deployed

The app requires the `data/` folder with RDS files. These files **must** be included in the deployment.

### Step 1: Prepare for Deployment

Make sure you're on the Trillium cluster with all data files generated:

```bash
cd /scratch/nxu/astrocytes/tx_classification_viewer

# Verify all RDS files exist
ls -lh data/*.rds
```

You should see:
- `gencode.rds` (12M)
- `ribotie.rds` (1.9M)
- `orfanage.rds` (4.3M)
- `pbid_to_pr_transcripts.rds` (487K)
- `pbid_to_orfanage_template.rds` (449K)
- `orf_type_lookup.rds` (179K)
- `pclass_lookup.rds` (354K)

### Step 2: Update .gitignore (Critical!)

**For shinyapps.io deployment**, you need to remove `data/*.rds` from `.gitignore` so the files are included.

Edit `.gitignore` to look like:

```
# Comment out data files for deployment
# data/*.rds

# Keep other entries
rsconnect/
.Rhistory
.RData
```

Or remove this line entirely if you don't mind committing RDS files to git.

### Step 3: Deploy to shinyapps.io

From your local machine (with R and the remote dev directory accessible):

```r
# Install rsconnect if needed
install.packages("rsconnect")

library(rsconnect)

# Set your shinyapps.io credentials (do this once)
setAccountInfo(account = "your-account",
               token = "your-token",
               secret = "your-secret")

# Deploy from the tx_classification_viewer directory
setwd("tx_classification_viewer")

# Deploy with explicit file list to include data/
rsconnect::deployApp(
  appDir = ".",
  appName = "transcript-viewer",
  appTitle = "Transcript Classification Viewer",
  account = "your-account"
  # appFiles parameter ensures data/ is included
)
```

### Step 4: Monitor the Deployment

Watch the deployment progress. It may take 1-2 minutes due to the RDS file sizes (20MB total).

---

## Troubleshooting

### Error: "Unable to connect to worker after 60.00 seconds"

**Cause:** Data files were not deployed.

**Solution:**
1. Check that `data/` folder is NOT in `.gitignore`
2. Redeploy with: `rsconnect::deployApp(appDir=".", appFiles=c("app.R", "data"))`
3. Check deployment logs on shinyapps.io dashboard

### Error: "data/ directory not found"

**Cause:** RDS files missing from server.

**Solution:**
1. Verify files exist locally: `ls data/*.rds`
2. Force redeploy: `rsconnect::forceDeployApp(appDir=".")`
3. Or delete the deployed app and redeploy from scratch

### Slow Deployment

The app is ~20MB due to the RDS files. Deployment typically takes 1-2 minutes on first upload.

---

## File Size Notes

| Component | Size | Purpose |
|-----------|------|---------|
| GENCODE GTF (RDS) | 12M | Reference transcripts |
| RiboTIE GTF (RDS) | 1.9M | Predicted ORFs |
| ORFanage GTF (RDS) | 4.3M | ORF predictions |
| Mappings + Lookups | 1.5M | Metadata tables |
| **Total** | **~20M** | All data needed |

This is much smaller than the original GTF files (~500MB+) and loads instantly.

---

## Alternative: Deploy Only Selected Transcripts

If file size is still an issue, you can pre-filter the data to include only transcripts of interest:

```r
# In setup_data.R, modify loading to filter
gencode <- import(gencode_gtf) %>%
    as_tibble() %>%
    select(c("seqnames", "start", "end", "strand", "transcript_id", "gene_id", "type")) %>%
    filter(transcript_id %in% selected_transcripts) %>%  # Add filter
    mutate(source = "GENCODE")
```

This would reduce the GENCODE RDS from 12M to a few MB.

---

## Verify Deployment Works

After deployment, check:

1. App loads without timeout
2. Transcript selector populates correctly
3. Plot renders when you select a transcript
4. Navigation buttons work smoothly
5. No "data/ directory not found" errors in logs
