# =====================================================================
# 02_estimates.R  -  survey-weighted estimates + Rao-Wu bootstrap CIs  -> output/results_long.csv
# =====================================================================
source("00_config.R"); source("01_functions.R")
set.seed(SEED)
Ds <- lapply(DATA_FILES, prep)
stdw_all <- pooled_stdw(Ds, "bi_all", length(BINS_ALL$lo))
stdw_nr  <- pooled_stdw(Ds, "bi_nr",  length(BINS_NR$lo))
cat("Standard population (8 bands):", round(stdw_all, 4), "\n")

results <- list(); extras <- list()
for (nm in names(Ds)) {
  D <- add_quintiles(Ds[[nm]]); w0 <- D$w0
  pt  <- stats(D, w0, stdw_all, stdw_nr, D$quint)
  fac <- make_reps(D, B_REPS)
  reps <- t(vapply(seq_len(B_REPS), function(b) stats(D, w0 * fac[D$psu, b], stdw_all, stdw_nr, D$quint), numeric(length(pt))))
  se <- apply(reps, 2, sd)
  results[[nm]] <- data.frame(state = nm, key = names(pt), est = unname(pt), se = unname(se),
                              lo = unname(pt - 1.96 * se), hi = unname(pt + 1.96 * se), row.names = NULL)
  # ---- point-only extras: pyramid, Lorenz curve, sample sizes
  d <- D$df; o <- D$o_m; xs <- D$mpce[o]; ws <- w0[o]
  sexed <- D$b3c4 %in% 1:2
  pyr <- sapply(1:2, function(code) sapply(1:17, function(g) sum(w0[D$ag == g & D$b3c4 == code]) / sum(w0[sexed]) * 100))
  colnames(pyr) <- c("M", "F")
  cwp <- cumsum(ws) / sum(ws); cwx <- cumsum(ws * xs) / sum(ws * xs); grid <- seq(0, 1, length.out = 101)
  lorenz <- data.frame(x = grid, y = approx(c(0, cwp), c(0, cwx), xout = grid)$y)
  extras[[nm]] <- list(
    pyramid = pyr, lorenz = lorenz, quintile_cuts = D$q_cuts,
    n = c(persons = nrow(d), households = sum(d$hh1), fsu = length(unique(d$fsu)),
          persons_rural = sum(d$sec == 1), persons_urban = sum(d$sec == 2),
          women15_49 = sum(d$b3c4 == 2 & d$b3c5 >= 15 & d$b3c5 <= 49), children0_4 = sum(d$b3c5 <= 4),
          girls0_4 = sum(d$b3c5 <= 4 & d$b3c4 == 2), transgender = sum(d$b3c4 == 3), missing_health = sum(is.na(d$b3c14)),
          hosp_unweighted = sum(D$hosp), hosp_unweighted_share_pct = 100 * mean(D$hosp)),
    weighted_total_persons = sum(w0), second_stage_stratum_hh = table(d$sss[d$hh1]))
  cat(nm, "done\n")
}
res <- do.call(rbind, results)
write.csv(res, file.path(OUT_DIR, "results_long.csv"), row.names = FALSE)
saveRDS(list(res = res, extras = extras, stdw_all = stdw_all, stdw_nr = stdw_nr), file.path(OUT_DIR, "results.rds"))
print(subset(res, key %in% c("median_age", "dep_total", "sex_ratio", "tfi", "smam_F", "gini", "cov_insured", "ci_insured", "mpce_median")), digits = 4)
