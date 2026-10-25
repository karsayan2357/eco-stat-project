# Run the whole analysis in order. setwd() to this folder first.
# Total run time with 500 bootstrap replicates is a few minutes.
for (f in c("02_estimates.R", "03_regression.R", "04_figures.R", "05_differences_tables.R")) {
  cat("\n=====", f, "=====\n"); source(f, echo = FALSE)
}
