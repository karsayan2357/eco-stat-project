# =====================================================================
# 00_config.R  -  settings shared by all scripts
# Project: Rajasthan vs Sikkim, NSS 80th Round (2025), Schedule 25.0 Health
# Put rajasthan.csv and sikkim.csv in the folder "data/" next to these scripts,
# then run  source("run_all.R")  from this folder (setwd() here first).
# =====================================================================
DATA_FILES <- c(Rajasthan = "data/rajasthan.csv", Sikkim = "data/sikkim.csv")
OUT_DIR    <- "output"          # tables (csv/xlsx) and figures are written here
B_REPS     <- as.integer(Sys.getenv("B_REPS", "500"))   # bootstrap replicates (report used 500; set B_REPS=20 for a quick test)
SEED       <- 20261025
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figs"), showWarnings = FALSE)

# age bands used for direct standardisation
BINS_ALL <- list(lo = c(0, 5, 15, 25, 35, 45, 55, 65), hi = c(4, 14, 24, 34, 44, 54, 64, 200))
# bands excluding ages 15-49 (childbirth ages), used for non-reproductive hospitalisation
BINS_NR  <- list(lo = c(0, 5, 50, 60, 70),             hi = c(4, 14, 59, 69, 200))
SEC_OK   <- c(6, 7, 8, 10, 11, 12, 13)   # b3c7 codes: secondary or higher
