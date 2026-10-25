# =====================================================================
# 05_differences_tables.R - state differences (Sikkim minus Rajasthan), scorecard (Table 12), Excel workbook
# =====================================================================
source("00_config.R")
r <- readRDS(file.path(OUT_DIR, "results.rds")); res <- r$res
rj <- res[res$state == "Rajasthan", ]; sk <- res[res$state == "Sikkim", ]
m <- merge(rj, sk, by = "key", suffixes = c("_RJ", "_SK"), sort = FALSE)
diffs <- data.frame(key = m$key, est_RJ = m$est_RJ, est_SK = m$est_SK, diff = m$est_SK - m$est_RJ)
diffs$se <- sqrt(m$se_RJ^2 + m$se_SK^2)                       # the two samples are independent
diffs$lo <- diffs$diff - 1.96 * diffs$se; diffs$hi <- diffs$diff + 1.96 * diffs$se
diffs$clear <- abs(diffs$diff / diffs$se) > 1.96
write.csv(diffs, file.path(OUT_DIR, "state_differences.csv"), row.names = FALSE)

sc_keys <- c(median_age = "Median age (years)", dep_total = "Total dependency ratio", tfi = "Fertility index (TFI)",
             smam_F = "Women's age at marriage (SMAM, years)", gpi_sec25 = "Female/male secondary-plus parity",
             chronic_std = "Chronic ailment, age-std.", hosp_std = "Hospitalisation, age-std.", cov_insured = "Any health scheme",
             mpce_median = "Median MPCE (Rs)", gini = "Gini coefficient")
score <- diffs[match(names(sc_keys), diffs$key), c("key", "diff", "lo", "hi", "clear")]
score$indicator <- unname(sc_keys)
score$verdict <- ifelse(!score$clear, "No clear difference", ifelse(score$diff > 0, "Higher in Sikkim", "Lower in Sikkim"))
write.csv(score, file.path(OUT_DIR, "scorecard_table12.csv"), row.names = FALSE)
print(score[, c("indicator", "diff", "lo", "hi", "verdict")], digits = 3)

# one Excel workbook with everything (needs openxlsx: install.packages("openxlsx"))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  sheets <- list(estimates = res, differences = diffs, scorecard = score)
  rf <- file.path(OUT_DIR, "regression_or.csv"); if (file.exists(rf)) sheets$regression_OR <- read.csv(rf)
  openxlsx::write.xlsx(sheets, file.path(OUT_DIR, "NSS80_Rajasthan_Sikkim_tables.xlsx"))
  cat("Excel workbook written\n")
}
