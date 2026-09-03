#!/usr/bin/env Rscript
# Quick deployment script for shinyapps.io

library(rsconnect)

cat("=== Transcript Classification Viewer - shinyapps.io Deployment ===\n\n")

# Check if data files exist
data_files <- list.files("data", pattern = "\\.rds$", full.names = TRUE)
if (length(data_files) == 0) {
  stop("ERROR: No RDS files found in data/ directory.\n",
       "Run 'Rscript setup_data.R' first to generate data files.")
}

cat("✓ Found", length(data_files), "RDS files\n")
cat("✓ Total data size:",
    format(sum(file.size(data_files)), units = "Mb"), "\n\n")

# Get user input
cat("Enter your shinyapps.io account name: ")
account <- trimws(readline())

cat("Enter your app name (for the URL): ")
app_name <- trimws(readline())

cat("Enter app title (optional, press Enter for default): ")
app_title <- trimws(readline())
if (app_title == "") {
  app_title <- "Transcript Classification Viewer"
}

cat("\n=== Deployment Details ===\n")
cat("Account:", account, "\n")
cat("App name:", app_name, "\n")
cat("App title:", app_title, "\n")
cat("Files to deploy: app.R + data/\n\n")

cat("Ready to deploy? (yes/no): ")
confirm <- trimws(tolower(readline()))

if (confirm != "yes") {
  cat("Deployment cancelled.\n")
  quit(save = "no", status = 0)
}

cat("\nDeploying...\n")

tryCatch({
  rsconnect::deployApp(
    appDir = ".",
    appName = app_name,
    appTitle = app_title,
    account = account,
    launch.browser = FALSE
  )

  cat("\n✓ Deployment complete!\n")
  cat("Your app is available at: https://", account, ".shinyapps.io/", app_name, "/\n", sep = "")

}, error = function(e) {
  cat("\n✗ Deployment failed:\n")
  cat(conditionMessage(e), "\n\n")
  cat("Troubleshooting:\n")
  cat("1. Make sure you've set your shinyapps.io credentials:\n")
  cat("   rsconnect::setAccountInfo(account='...', token='...', secret='...')\n")
  cat("2. Check that data/*.rds files exist locally\n")
  cat("3. See DEPLOYMENT.md for detailed instructions\n")
})
