# Fix: Out of Memory Error on shinyapps.io

Your app is failing with: `"Container event from container-14326642: oom (out of memory)"`

This is happening because the app loads ALL GTF data at startup, consuming too much RAM for shinyapps.io's containers.

## 🚀 Quick Fix (5 minutes)

### Step 1: Regenerate Optimized Data

```bash
cd /scratch/nxu/astrocytes/tx_classification_viewer

# Remove old data files
rm data/*.rds

# Regenerate with optimized filtering (50% smaller)
Rscript setup_data_optimized.R
```

This creates smaller RDS files (9.3MB instead of 20MB) by filtering to only relevant transcripts.

### Step 2: Redeploy

```r
setwd("tx_classification_viewer")

# Use the interactive deploy script
source("deploy.R")

# Or manual rsconnect:
rsconnect::forceDeployApp(appDir = ".", appName = "your-app-name")
```

That's it! The app should now deploy and run successfully.

---

## Why This Fixes It

| Issue | Root Cause | Solution |
|-------|------------|----------|
| OOM Error | Loading 20MB of GTF data at startup | Filter to only 9.3MB of used transcripts |
| 60s timeout | High memory usage during initialization | Smaller dataset loads faster |
| Container crash | RAM exceeded during app startup | Use optimized setup script |

---

## File Size Comparison

```
Original Setup
├── GENCODE GTF: 12MB (all ~500k transcripts)
├── ORFanage: 4.3MB
├── RiboTIE: 1.9MB
└── Metadata: 1.5MB
Total: 20MB ❌ OOM on shinyapps.io

Optimized Setup
├── GENCODE: 1.8MB (only ~118k referenced transcripts) ✅
├── ORFanage: 4.2MB
├── RiboTIE: 1.9MB
└── Metadata: 1.5MB
Total: 9.3MB ✅ Works on shinyapps.io free tier
```

---

## Verification

After redeployment, verify it works:

1. ✅ App loads without timeout (should take 5-10 seconds)
2. ✅ No "out of memory" in logs
3. ✅ Transcript dropdown populates
4. ✅ Can select a transcript and plot renders
5. ✅ Navigation buttons work

---

## If It Still Fails

### Check 1: Verify data was regenerated

```bash
ls -lh data/*.rds
# Should see 9 files, ~9.3MB total
```

### Check 2: Verify correct files were deployed

In shinyapps.io dashboard → "Deployment Details", confirm:
- ✅ `app.R` present
- ✅ `data/` folder included with RDS files
- ❌ No original GTF files (shouldn't be deployed)

### Check 3: Check deployment logs

View logs in shinyapps.io dashboard to see:
- Package loading times
- Data loading success/failure
- Memory usage

If you see "Loading data..." and then timeout, the RDS files weren't deployed. Make sure `data/` folder is included in deployment.

---

## Alternative: Lazy Loading (Advanced)

If optimized setup still doesn't work, try lazy loading:

```bash
cd tx_classification_viewer
cp app.R app_eager.R
cp app_lazy.R app.R
rsconnect::forceDeployApp()
```

This defers data loading until first use, reducing startup memory.

---

## Prevention for Future

The optimized setup is now the default:
- Use `setup_data_optimized.R` going forward
- Always deploy with `data/` folder included
- Monitor shinyapps.io logs if you make changes

For more details, see:
- `MEMORY_OPTIMIZATION.md` - Deep dive into both approaches
- `DEPLOYMENT.md` - Full deployment guide
