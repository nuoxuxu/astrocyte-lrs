# App Versions Guide

You have multiple app versions optimized for different scenarios. Choose the right one based on your needs.

## Quick Decision Tree

```
Are you getting OOM errors?
│
├─ YES (on startup)
│  └─ Use OPTIMIZED data + app_memory_efficient.R
│
├─ YES (when selecting transcript)
│  └─ Use app_memory_efficient.R (required)
│
└─ NO (works fine)
   ├─ Local use → Use app.R (fastest)
   └─ shinyapps.io → Use app_memory_efficient.R (safest)
```

## Available Versions

### 1. `app.R` (Original, Standard)

**Best for:** Local development, self-hosted with good RAM

```bash
shiny::runApp("app.R")
```

**Features:**
- ✅ Fastest startup (5s)
- ✅ Instant plot rendering
- ✅ Instant navigation
- ❌ High memory usage (400-500MB)
- ❌ May OOM on shinyapps.io

**Data requirements:**
- Works with both `setup_data.R` and `setup_data_optimized.R`
- Prefer optimized for smaller RAM

**When to use:**
- Running locally
- Server has 2GB+ RAM
- Want absolute fastest performance

---

### 2. `app_memory_efficient.R` (Recommended for shinyapps.io) ⭐

**Best for:** shinyapps.io, limited resources, reliability

```bash
cp app_memory_efficient.R app.R
shiny::runApp()
```

**Features:**
- ✅ Low memory usage (150-200MB peak)
- ✅ Works reliably on free tier
- ✅ Instant navigation (cached)
- ✅ Lazy loads annotation data
- ⚠️ Slightly slower first plot (2-3s vs instant)

**Data requirements:**
- Must use `setup_data_optimized.R`
- Data folder required: 9.3MB

**When to use:**
- Deploying to shinyapps.io ✅
- Limited server RAM (< 2GB)
- Need reliability over speed
- Getting OOM errors ✅

**How it works:**
1. Startup: Load only essential data (10MB)
2. First transcript: Lazy load full annotation (3s delay)
3. Subsequent: All cached and instant

---

### 3. `app_lazy.R` (Ultra-Lightweight) 

**Best for:** Very limited RAM, slow startup acceptable

```bash
cp app_lazy.R app.R
shiny::runApp()
```

**Features:**
- ✅ Minimal startup memory (5-10MB)
- ✅ Maximum RAM efficiency
- ⚠️ Slower first plot (8-10s)
- ✅ Very stable

**Data requirements:**
- Works with `setup_data.R` (20MB) or optimized
- All data kept on disk until needed

**When to use:**
- Startup timeout issues
- Extremely limited memory (<512MB)
- Memory-efficient is still too heavy
- Slow first load acceptable

**How it works:**
```
Start app (5s)
  ↓
Select transcript
  ↓
Load annotation from disk (8-10s)
  ↓
Render plot
  ↓
All cached afterward (instant)
```

---

### 4. `app_eager.R` (Backup)

Backup of original if you rename it during testing.

---

### 5. `app_original.R` (Backup)

Another backup created during troubleshooting.

## Version Comparison Table

| Feature | app.R | app_memory_efficient.R | app_lazy.R |
|---------|-------|----------------------|-----------|
| **Startup** | 5s | 10s | 2s |
| **First plot** | Instant | 2-3s | 8-10s |
| **Subsequent plots** | Instant | Instant | Instant |
| **Memory startup** | 100MB | 10MB | 5MB |
| **Memory peak** | 400-500MB | 150-200MB | 150-200MB |
| **shinyapps.io free** | ❌ OOM | ✅ Reliable | ⚠️ Risky |
| **Local dev** | ✅ Best | ⚠️ Slower | ⚠️ Slow |
| **Self-hosted** | ✅ Best | ✅ Good | ✅ Safe |

## Setup Data Versions

### `setup_data.R` (Original)
- Creates 20MB data files
- No filtering
- For self-hosted servers only

### `setup_data_optimized.R` (Recommended) ⭐
- Creates 9.3MB data files
- Filters to only used transcripts
- Required for shinyapps.io

## Deployment Scenarios

### Scenario 1: Local Development
```bash
# 1. Generate data (either version fine)
Rscript setup_data_optimized.R

# 2. Use original app for speed
shiny::runApp("app.R")

# 3. Edit and test locally
```

### Scenario 2: shinyapps.io (Recommended Path)
```bash
# 1. Generate optimized data
rm data/*.rds
Rscript setup_data_optimized.R

# 2. Use memory-efficient version
cp app_memory_efficient.R app.R

# 3. Deploy
source("deploy.R")
```

### Scenario 3: Getting OOM on shinyapps.io
```bash
# 1. If data is 20MB, regenerate optimized
rm data/*.rds
Rscript setup_data_optimized.R

# 2. Switch to memory-efficient app
cp app_memory_efficient.R app.R

# 3. Redeploy
rsconnect::forceDeployApp()
```

### Scenario 4: Still Getting OOM
```bash
# Try lazy loading
cp app_lazy.R app.R
rsconnect::forceDeployApp()

# First select will be slow, but should work
```

## Testing Checklist

Before deploying, test locally:

```r
setwd("tx_classification_viewer")

# Test 1: Startup
shiny::runApp()  # Should start in 2-10s

# Test 2: Selection
# → Select "PB.11906.164_253"
# → Plot should render without lag/crash

# Test 3: Navigation  
# → Click "Next →" and "← Previous"
# → Should be instant

# Test 4: Reference change
# → Change "Reference to Highlight" radio button
# → Plot should update without lag/crash

# Test 5: Toggle focused
# → Check "Focused View (zoom on differences)"
# → Plot should update instantly
```

If all tests pass locally, safe to deploy.

## Memory Optimization Timeline

**Evolution of the app:**

1. **Original** (OOM on shinyapps.io)
   - All data loaded at startup
   - Peak: 500MB
   - ❌ Crashes

2. **Optimized Setup** (Better, but still risky)
   - Smaller data (9.3MB vs 20MB)
   - Peak: 400MB
   - ⚠️ May still OOM

3. **Memory Efficient** (Recommended)
   - Lazy loads annotation
   - Explicit garbage collection
   - Peak: 150-200MB
   - ✅ Works reliably

4. **Lazy Loading** (Maximum safety)
   - Defers all data loading
   - Peak: 150-200MB
   - ✅ Works even with tiny RAM
   - ⚠️ Slower startup

## Recommendation

**For shinyapps.io:**
1. Use `setup_data_optimized.R` (9.3MB data)
2. Use `app_memory_efficient.R` (recommended)
3. Deploy with `source("deploy.R")`

**If you get OOM:**
1. Verify data is optimized: `du -sh data/` → should be 9.3MB
2. Switch to `app_memory_efficient.R` if not already
3. If still failing, try `app_lazy.R`

**For local development:**
- Use `app.R` (fastest)
- Can use either data version

## Switching Versions

To switch versions:

```bash
cd tx_classification_viewer

# Backup current
cp app.R app_backup.R

# Switch to different version
cp app_memory_efficient.R app.R

# Or switch to lazy
cp app_lazy.R app.R

# Test locally
shiny::runApp()

# If good, redeploy
rsconnect::forceDeployApp()
```

## Questions?

See detailed guides:
- `MEMORY_OPTIMIZATION.md` - Deep dive on memory issues
- `FIX_OOM_ON_SELECT.md` - Specific fix for selection crashes
- `DEPLOYMENT.md` - Full deployment guide
