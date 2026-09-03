# Immediate Fix for OOM on Transcript Selection

## What's Happening

1. ✅ App loads (startup is fine)
2. ✅ Transcript dropdown works
3. ❌ When you select a transcript, server crashes with "oom"

This is because the **plot rendering** uses too much memory.

## Fix (3 Simple Steps)

### Step 1: Backup & Switch to Memory-Efficient Version

```bash
cd /scratch/nxu/astrocytes/tx_classification_viewer

# Backup current app.R
cp app.R app_original.R

# Use memory-efficient version (optimized for low RAM)
cp app_memory_efficient.R app.R
```

### Step 2: Verify Data (Make Sure You Used Optimized Setup)

```bash
# Check data file sizes
ls -lh data/*.rds | awk '{print $5, $9}'

# Should look like:
# 1.8M gencode.rds
# 4.2M orfanage.rds
# 1.9M ribotie.rds
# ... (others around 200K-500K)
# Total: ~9.3M
```

**If total is 20M instead of 9.3M:**

```bash
# Regenerate with optimized filtering
rm data/*.rds
Rscript setup_data_optimized.R
```

### Step 3: Redeploy

```r
setwd("tx_classification_viewer")

# Option A: Interactive deploy script
source("deploy.R")

# Option B: Manual rsconnect (if source() doesn't work)
library(rsconnect)
rsconnect::forceDeployApp(appDir = ".")
```

**Done!** 🚀 The app should now work without OOM crashes.

---

## What Changed

The `app_memory_efficient.R` version:

1. **Lazy loads annotation data** - Only loads full GTF when first transcript is selected
2. **Forces garbage collection** - Cleans up memory before and after plotting
3. **Removes temporary objects** - Doesn't keep intermediate data in memory
4. **Same functionality** - Everything works exactly the same, just uses less RAM

## Expected Behavior After Fix

```
Select transcript #1
  ↓
First plot takes 2-3 seconds (loading annotation data)
  ↓
Select transcript #2
  ↓
Plot renders instantly (data already loaded & cached)
  ↓
Change reference, toggle focused, navigate differences
  ↓
All instant (no more OOM)
```

## Verification

After redeploying, test it works:

1. Open the app URL on shinyapps.io
2. Select a transcript from dropdown
3. Plot should render without timeout
4. Select another transcript
5. Plot should be instant (data cached)
6. Use all features without crashes

## Still Having Issues?

### If deployment fails:

```bash
# Check that files exist
ls -la tx_classification_viewer/app.R
ls -la tx_classification_viewer/data/*.rds

# Check they're readable
wc -l app.R app_memory_efficient.R
```

### If app still crashes on select:

Try the ultra-lightweight lazy-loading version:

```bash
cp app_lazy.R app.R
rsconnect::forceDeployApp()

# This is slower on first select (8-10s) but more memory-efficient
```

### If you're unsure what version is running:

Check shinyapps.io logs for:

```
# Memory-efficient version shows:
"Core data loaded. Annotation data will load on demand."

# Lazy version shows:
"Initializing app (lazy loading enabled)..."

# Original shows neither message
```

## Files Changed

Only these need to be redeployed:
- `app.R` ← Changed from app.R to app_memory_efficient.R
- `data/*.rds` ← Should already be 9.3M (optimized)

Everything else stays the same.

---

## Summary

| Step | Action | Expected Result |
|------|--------|-----------------|
| 1 | `cp app_memory_efficient.R app.R` | Switched to memory-safe version |
| 2 | Verify `du -sh data/` is 9.3M | Using optimized data |
| 3 | `rsconnect::forceDeployApp()` | App redeployed |

After this, selecting transcripts should work without OOM! 

If not, the **one-line** fix is to try lazy loading:
```bash
cp app_lazy.R app.R && rsconnect::forceDeployApp()
```
