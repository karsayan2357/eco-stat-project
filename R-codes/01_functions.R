# =====================================================================
# 01_functions.R  -  data preparation and estimators (base R only)
# =====================================================================
eq <- function(x, v) !is.na(x) & x == v

bin_index <- function(age, bins) {              # 1..K, NA if outside all bands
  idx <- rep(NA_integer_, length(age))
  for (i in seq_along(bins$lo)) idx[age >= bins$lo[i] & age <= bins$hi[i]] <- i
  idx
}
bsum <- function(x, g, nb) {                     # sum of x by integer group g (1..nb)
  r <- rowsum(x, g); out <- numeric(nb); out[as.integer(rownames(r))] <- r[, 1]; out
}

# ---------------------------------------------------------------- data preparation
prep <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$hh1  <- !duplicated(d[c("fsu", "samp", "sss", "hhd", "sec")])    # first row of each household
  d$mpce <- d$umce / d$hhsz                                          # monthly per-capita expenditure
  D <- as.list(d)
  D$df <- d; D$n <- nrow(d)
  D$w0 <- as.numeric(d$final_weight)                                 # final_weight = mult / 100
  D$age <- d$b3c5
  D$female <- eq(d$b3c4, 2); D$male <- eq(d$b3c4, 1)
  D$rural <- d$sec == 1;     D$urban <- d$sec == 2
  D$valid   <- !is.na(d$b3c14)                                       # ailment/coverage fields present
  D$chronic <- eq(d$b3c14, 1); D$ail15 <- eq(d$b3c15, 1)
  D$hosp    <- eq(d$b3c9, 1);  D$vacc  <- eq(d$b3c8, 1)
  D$insured  <- !is.na(d$b3c17) & d$b3c17 != 19
  D$pmjay    <- eq(d$b3c17, 1); D$stscheme <- eq(d$b3c17, 2)
  D$emp      <- !is.na(d$b3c17) & d$b3c17 %in% 3:7                   # ESIS, CGHS/ECHS, state emp., PSU, other employer
  D$stemp    <- eq(d$b3c17, 5); D$private <- eq(d$b3c17, 10)
  D$comm     <- !is.na(d$b3c13) & d$b3c13 != 19
  D$preg     <- eq(d$b3c11, 1); D$nevermar <- eq(d$b3c6, 1)
  D$edu      <- d$b3c7
  D$ag       <- d$age_group                                          # 1 = 0-4, 2 = 5-9, ..., 17 = 80+
  D$psu      <- match(d$fsu, unique(d$fsu))                          # design: PSU = FSU
  D$strat    <- d$strat_var                                          # design: stratum
  D$o_m      <- order(d$mpce)                                        # ranking by MPCE
  ms         <- d$mpce[D$o_m]
  D$gid_m    <- cumsum(c(TRUE, ms[-1] != ms[-length(ms)]))           # tie groups in sorted order
  D$o_a      <- order(D$age)
  D$bi_all   <- bin_index(D$age, BINS_ALL)
  D$bi_nr    <- bin_index(D$age, BINS_NR)
  D
}

# ---------------------------------------------------------------- weighted helpers
wmean <- function(y, w, m = NULL) {
  if (is.null(m)) { den <- sum(w); num <- sum(w * y) } else { den <- sum(w * m); num <- sum(w * y * m) }
  if (den > 0) num / den else NA_real_
}
wq_sorted <- function(x_s, w_s, probs) {        # weighted quantile (first value whose cum. weight reaches p * total)
  cw <- cumsum(w_s); tot <- cw[length(cw)]
  idx <- pmin(findInterval(probs * tot, cw, left.open = TRUE) + 1L, length(x_s))
  x_s[idx]
}
std_rate <- function(y, w, mask, bidx, stdw) {   # direct age standardisation
  nb <- length(stdw); sel <- mask & !is.na(bidx)
  num <- bsum((w * y)[sel], bidx[sel], nb); den <- bsum(w[sel], bidx[sel], nb)
  ok <- den > 0
  if (!any(ok)) return(NA_real_)
  sum(stdw[ok] * num[ok] / den[ok]) / sum(stdw[ok])
}
gini_sorted <- function(x_s, w_s) {
  W <- sum(w_s); if (W <= 0) return(NA_real_)
  F <- (cumsum(w_s) - 0.5 * w_s) / W
  2 * sum(w_s * x_s * F) / sum(w_s * x_s) - 1
}
frac_rank <- function(ws, gid) {                 # fractional rank, ties share the mid-rank
  gtot <- as.numeric(rowsum(ws, gid)[, 1])
  before <- cumsum(gtot) - gtot
  ((before + 0.5 * gtot) / sum(ws))[gid]
}
conc_index <- function(h, w, D, standardise = FALSE) {  # concentration index over MPCE rank
  o <- D$o_m; h <- as.numeric(h)
  if (standardise) {                              # indirect standardisation by age-group x sex cell
    cell <- D$ag * 4 + D$b3c4
    num <- bsum(w * h, cell, max(cell)); den <- bsum(w, cell, max(cell))
    cmean <- ifelse(den > 0, num / pmax(den, 1e-12), 0)
    mu <- sum(w * h) / sum(w)
    h <- h - cmean[cell] + mu
  }
  hs <- h[o]; ws <- w[o]
  r <- frac_rank(ws, D$gid_m)
  mu <- sum(ws * hs) / sum(ws); rbar <- sum(ws * r) / sum(ws)
  cov <- sum(ws * (hs - mu) * (r - rbar)) / sum(ws)
  2 * cov / mu
}
pooled_stdw <- function(Ds, key, nb) {          # standard population = mean of the two states' age shares
  rowMeans(sapply(Ds, function(D) {
    m <- !is.na(D[[key]]); b <- bsum(D$w0[m], D[[key]][m], nb); b / sum(b)
  }))
}
add_quintiles <- function(D) {                   # person-weighted MPCE quintiles within the state
  o <- D$o_m; cuts <- wq_sorted(D$mpce[o], D$w0[o], c(.2, .4, .6, .8))
  D$quint <- findInterval(D$mpce, cuts, left.open = TRUE) + 1L
  D$q_cuts <- cuts
  D
}

# ---------------------------------------------------------------- all point estimates for one weight vector
stats <- function(D, w, stdw_all, stdw_nr, quint) {
  o <- D$o_m; s <- list()
  F <- D$female; M <- D$male; age <- D$age; ag <- D$ag
  wv <- w * D$valid
  # ---- demography
  s$sex_ratio <- 1000 * sum(w[F]) / sum(w[M])
  for (lab in c("rural", "urban")) { m <- D[[lab]]; s[[paste0("sex_ratio_", lab)]] <- 1000 * sum(w[F & m]) / sum(w[M & m]) }
  for (b in list(c("0_4", 0, 4), c("0_14", 0, 14), c("15_59", 15, 59), c("60p", 60, 200))) {
    m <- age >= as.numeric(b[2]) & age <= as.numeric(b[3])
    s[[paste0("sex_ratio_", b[1])]] <- 1000 * sum(w[F & m]) / sum(w[M & m])
  }
  oa <- D$o_a
  s$median_age <- wq_sorted(age[oa], w[oa], 0.5)
  s$mean_age   <- wmean(age, w)
  sh <- list()
  for (b in list(c("0_14", 0, 14), c("15_59", 15, 59), c("60p", 60, 200))) {
    sh[[b[1]]] <- wmean(as.numeric(age >= as.numeric(b[2]) & age <= as.numeric(b[3])), w)
    s[[paste0("share_", b[1])]] <- sh[[b[1]]]
  }
  s$dep_young <- sh[["0_14"]] / sh[["15_59"]]; s$dep_old <- sh[["60p"]] / sh[["15_59"]]
  s$dep_total <- (sh[["0_14"]] + sh[["60p"]]) / sh[["15_59"]]
  s$rural_share <- wmean(as.numeric(D$rural), w)
  hw <- w * D$hh1                                                      # household-level weight
  s$hhsize <- wmean(as.numeric(D$hhsz), hw)
  s$female_head <- wmean(as.numeric(F), w, D$b3c3 == 1)
  for (r in list(c("hindu", 1), c("muslim", 2), c("christian", 3), c("sikh", 4), c("jain", 5), c("buddhist", 6)))
    s[[paste0("rel_", r[1])]] <- wmean(as.numeric(D$b5i2 == as.numeric(r[2])), w)
  for (r in list(c("st", 1), c("sc", 2), c("obc", 3), c("gen", 9)))
    s[[paste0("soc_", r[1])]] <- wmean(as.numeric(D$b5i3 == as.numeric(r[2])), w)
  # ---- education
  edu <- D$edu
  for (sx in c("F", "M")) {
    m <- if (sx == "F") F else M; a25 <- m & age >= 25; a15 <- m & age >= 15
    s[[paste0("notlit25_", sx)]] <- wmean(as.numeric(eq(edu, 1)), w, a25)
    s[[paste0("sec25_", sx)]]    <- wmean(as.numeric(edu %in% SEC_OK), w, a25)
    s[[paste0("lit15_", sx)]]    <- wmean(as.numeric(!eq(edu, 1)), w, a15)
    s[[paste0("grad25_", sx)]]   <- wmean(as.numeric(edu %in% c(12, 13)), w, a25)
  }
  s$gpi_sec25 <- s$sec25_F / s$sec25_M; s$gpi_lit15 <- s$lit15_F / s$lit15_M
  # ---- singulate mean age at marriage (Hajnal)
  for (sx in c("F", "M")) {
    m <- if (sx == "F") F else M
    sa <- sapply(4:11, function(g) wmean(as.numeric(D$nevermar), w, m & ag == g)); names(sa) <- 4:11
    s50 <- 0.5 * (sa["10"] + sa["11"])
    s[[paste0("smam_", sx)]] <- unname((15 + 5 * sum(sa[as.character(4:10)]) - 50 * s50) / (1 - s50))
    s[[paste0("evermar15_19_", sx)]] <- unname(1 - sa["4"]); s[[paste0("evermar20_24_", sx)]] <- unname(1 - sa["5"])
    s[[paste0("never50_", sx)]] <- unname(s50)
  }
  s$smam_gap <- s$smam_M - s$smam_F
  # ---- pregnancy-based fertility index (proxy for TFR)
  asfr <- sapply(4:10, function(g) wmean(as.numeric(D$preg), w, F & ag == g))
  for (i in seq_along(asfr)) s[[paste0("asfr_", i + 3)]] <- asfr[i]
  s$tfi <- 5 * sum(asfr); s$mac <- sum((17.5 + 5 * (0:6)) * asfr) / sum(asfr)
  # ---- health indicators: age-standardised rates and female/male ratios
  for (k in c("chronic", "ail15", "hosp", "comm")) {
    msk <- if (k == "hosp") rep(TRUE, D$n) else D$valid
    y <- as.numeric(D[[k]])
    s[[paste0(k, "_std")]]   <- std_rate(y, w, msk & (F | M), D$bi_all, stdw_all)
    s[[paste0(k, "_crude")]] <- wmean(y, w, msk)
    for (sx in c("F", "M")) {
      m <- if (sx == "F") F else M
      s[[paste0(k, "_std_", sx)]]   <- std_rate(y, w, msk & m, D$bi_all, stdw_all)
      s[[paste0(k, "_crude_", sx)]] <- wmean(y, w, msk & m)
    }
    s[[paste0(k, "_ratioFM")]] <- s[[paste0(k, "_std_F")]] / s[[paste0(k, "_std_M")]]
  }
  y <- as.numeric(D$hosp)
  for (sx in c("F", "M")) { m <- if (sx == "F") F else M; s[[paste0("hospnr_std_", sx)]] <- std_rate(y, w, m, D$bi_nr, stdw_nr) }
  s$hospnr_ratioFM <- s$hospnr_std_F / s$hospnr_std_M
  s$hosp_w15_49 <- wmean(y, w, F & age >= 15 & age <= 49); s$hosp_m15_49 <- wmean(y, w, M & age >= 15 & age <= 49)
  for (sx in c("F", "M")) { m <- if (sx == "F") F else M; s[[paste0("vacc04_", sx)]] <- wmean(as.numeric(D$vacc), w, m & age <= 4) }
  s$vacc04_ratioFM <- s$vacc04_F / s$vacc04_M
  s$vacc04 <- wmean(as.numeric(D$vacc), w, age <= 4); s$vacc1_4 <- wmean(as.numeric(D$vacc), w, age >= 1 & age <= 4)
  # ---- health-scheme coverage
  for (k in c("insured", "pmjay", "stscheme", "emp", "stemp", "private")) {
    y <- as.numeric(D[[k]])
    s[[paste0("cov_", k)]] <- wmean(y, w, D$valid)
    for (sx in c("F", "M")) { m <- if (sx == "F") F else M; s[[paste0("cov_", k, "_", sx)]] <- wmean(y, w, D$valid & m) }
    for (lab in c("rural", "urban")) s[[paste0("cov_", k, "_", lab)]] <- wmean(y, w, D$valid & D[[lab]])
    for (r in list(c("st", 1), c("sc", 2), c("obc", 3), c("gen", 9)))
      s[[paste0("cov_", k, "_soc_", r[1])]] <- wmean(y, w, D$valid & D$b5i3 == as.numeric(r[2]))
  }
  s$cov_insured_ratioFM <- s$cov_insured_F / s$cov_insured_M
  s$cov_none <- 1 - s$cov_insured
  for (lab in c("ins", "unins")) {
    m <- if (lab == "ins") D$insured else !D$insured
    s[[paste0("hosp_std_", lab)]] <- std_rate(as.numeric(D$hosp), w, D$valid & m & (F | M), D$bi_all, stdw_all)
  }
  s$hosp_ratio_ins <- s$hosp_std_ins / s$hosp_std_unins
  for (lab in c("rural", "urban")) for (k in c("chronic", "hosp")) {
    msk <- if (k == "hosp") rep(TRUE, D$n) else D$valid
    s[[paste0(k, "_std_", lab)]] <- std_rate(as.numeric(D[[k]]), w, msk & D[[lab]] & (F | M), D$bi_all, stdw_all)
  }
  # ---- by MPCE quintile
  for (q in 1:5) {
    mq <- quint == q
    s[[sprintf("q%d_share_pop", q)]] <- wmean(as.numeric(mq), w)
    for (k in c("chronic", "hosp")) {
      msk <- if (k == "hosp") rep(TRUE, D$n) else D$valid
      s[[sprintf("q%d_%s_std", q, k)]] <- std_rate(as.numeric(D[[k]]), w, msk & mq & (F | M), D$bi_all, stdw_all)
    }
    for (k in c("insured", "pmjay", "stscheme", "emp", "private"))
      s[[sprintf("q%d_cov_%s", q, k)]] <- wmean(as.numeric(D[[k]]), w, D$valid & mq)
  }
  # ---- consumption and inequality (person-weighted MPCE)
  x <- D$mpce; xs <- x[o]; ws <- w[o]
  s$mpce_mean <- wmean(x, w)
  qs <- wq_sorted(xs, ws, c(.1, .2, .5, .8, .9))
  s$mpce_p10 <- qs[1]; s$mpce_p20 <- qs[2]; s$mpce_median <- qs[3]; s$mpce_p80 <- qs[4]; s$mpce_p90 <- qs[5]
  s$p90_p10 <- qs[5] / qs[1]; s$p80_p20 <- qs[4] / qs[2]
  s$gini <- gini_sorted(xs, ws)
  tot <- sum(ws * xs); cuts <- wq_sorted(xs, ws, c(.2, .4, .6, .8))
  qi <- findInterval(xs, cuts, left.open = TRUE)                       # 0..4, ties at a cut go lower
  for (q in 0:4) s[[sprintf("qshare_%d", q + 1)]] <- sum(ws * xs * (qi == q)) / tot
  s$top10_share    <- sum(ws * xs * (xs > wq_sorted(xs, ws, 0.9))) / tot
  s$bottom40_share <- sum(ws * xs * (xs <= wq_sorted(xs, ws, 0.4))) / tot
  for (lab in c("rural", "urban")) {
    wm <- ws * D[[lab]][o]
    s[[paste0("mpce_median_", lab)]] <- wq_sorted(xs, wm, 0.5)
    s[[paste0("mpce_mean_", lab)]]   <- sum(wm * xs) / sum(wm)
    s[[paste0("gini_", lab)]]        <- gini_sorted(xs, wm)
  }
  s$urb_rur_ratio <- s$mpce_median_urban / s$mpce_median_rural
  for (r in list(c("st", 1), c("sc", 2), c("obc", 3), c("gen", 9))) {
    wm <- ws * (D$b5i3[o] == as.numeric(r[2]))
    s[[paste0("mpce_median_soc_", r[1])]] <- if (sum(wm) > 0) wq_sorted(xs, wm, 0.5) else NA_real_
  }
  # ---- household-level expenditure composition
  U <- sum(hw * D$umce)
  s$hh_umce_mean <- U / sum(hw)
  s$comp_purchased <- sum(hw * D$b5i7) / U; s$comp_homegrown <- sum(hw * D$b5i8) / U; s$comp_inkind <- sum(hw * D$b5i9) / U
  s$comp_annual <- sum(hw * (D$b5i10 + D$b5i11) / 12) / U
  s$comp_clothing <- sum(hw * D$b5i10 / 12) / U; s$comp_durables <- sum(hw * D$b5i11 / 12) / U
  for (lab in c("rural", "urban")) {
    m <- D[[lab]]; Um <- sum(hw * m * D$umce)
    s[[paste0("comp_homegrown_", lab)]] <- sum(hw * m * D$b5i8) / Um
    s[[paste0("comp_annual_", lab)]]    <- sum(hw * m * (D$b5i10 + D$b5i11) / 12) / Um
    s[[paste0("comp_inkind_", lab)]]    <- sum(hw * m * D$b5i9) / Um
  }
  s$hh_premium_share <- wmean(as.numeric(D$b5i6 > 0), hw); s$hh_outbreak <- wmean(as.numeric(D$b5i5 == 1), hw)
  # ---- concentration indices (rank by MPCE)
  for (k in c("insured", "pmjay", "stscheme", "emp")) s[[paste0("ci_", k)]] <- conc_index(as.numeric(D[[k]]), wv, D)
  for (k in c("chronic", "hosp")) s[[paste0("ci_", k, "_std")]] <- conc_index(as.numeric(D[[k]]), if (k == "hosp") w else wv, D, standardise = TRUE)
  s$ci_mpce_health_ail15 <- conc_index(as.numeric(D$ail15), wv, D, standardise = TRUE)
  unlist(s)
}

# ---------------------------------------------------------------- Rao-Wu rescaled bootstrap factors
make_reps <- function(D, B) {
  npsu <- max(D$psu); psu_strat <- integer(npsu); psu_strat[D$psu] <- D$strat
  fac <- matrix(1, npsu, B)
  for (h in unique(psu_strat)) {
    idx <- which(psu_strat == h); n <- length(idx)
    if (n < 2) next
    m <- n - 1
    fac[idx, ] <- rmultinom(B, m, rep(1 / n, n)) * n / m           # resample n-1 PSUs with replacement, rescale
  }
  fac
}
