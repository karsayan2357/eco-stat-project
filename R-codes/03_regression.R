# =====================================================================
# 03_regression.R  -  weighted logistic regressions (odds ratios), FSU-clustered SEs  -> output/regression_or.csv
#   chronic      : chronic ailment, all ages, men and women
#   hosp_nonrep  : hospitalised, excluding ages 15-49 (childbirth ages)
#   hosp_all     : hospitalised, all ages
# Needs the 'sandwich' package: install.packages("sandwich")
# =====================================================================
source("00_config.R"); source("01_functions.R")
library(sandwich)
Ds <- lapply(DATA_FILES, prep)
age_bands <- BINS_ALL

design <- function(D) {
  X <- data.frame(female = as.numeric(D$female), rural = as.numeric(D$rural))
  bi <- bin_index(D$age, age_bands)
  for (i in 2:length(age_bands$lo)) {
    hi <- age_bands$hi[i]
    X[[sprintf("age_%d_%s", age_bands$lo[i], if (hi < 200) hi else "plus")]] <- as.numeric(!is.na(bi) & bi == i)
  }
  for (q in 2:5) X[[paste0("quintile_", q)]] <- as.numeric(D$quint == q)
  X$insured <- as.numeric(D$insured)
  X$soc_ST <- as.numeric(D$b5i3 == 1); X$soc_SC <- as.numeric(D$b5i3 == 2); X$soc_OBC <- as.numeric(D$b5i3 == 3)
  X
}

out <- list()
for (nm in names(Ds)) {
  D <- add_quintiles(Ds[[nm]]); X <- design(D); sexed <- D$female | D$male
  specs <- list(
    chronic     = list(y = D$chronic, mask = D$valid & sexed),
    hosp_nonrep = list(y = D$hosp,    mask = sexed & D$valid & (D$age < 15 | D$age >= 50)),
    hosp_all    = list(y = D$hosp,    mask = sexed & D$valid))
  for (lab in names(specs)) {
    sp <- specs[[lab]]; m <- sp$mask; Xm <- X[m, , drop = FALSE]
    if (lab == "hosp_nonrep") {   # sample has only ages <15 and >=50, so use bands 0-4 (ref), 5-14, 50-64, 65+
      a <- D$age[m]
      Xm <- Xm[, c("female", "rural", "age_5_14", paste0("quintile_", 2:5), "insured", "soc_ST", "soc_SC", "soc_OBC")]
      Xm$age_65_plus <- as.numeric(a >= 65); Xm$age_50_64 <- as.numeric(a >= 50 & a <= 64)
    }
    wm <- D$w0[m]; ws <- wm * length(wm) / sum(wm)              # weights rescaled to sum to the sample size
    dat <- cbind(y = as.numeric(sp$y[m]), Xm)
    fit <- glm(y ~ ., data = dat, family = quasibinomial(), weights = ws)
    V <- vcovCL(fit, cluster = D$psu[m], type = "HC1")          # cluster-robust by FSU
    b <- coef(fit); se <- sqrt(diag(V))
    out[[paste(nm, lab)]] <- data.frame(state = nm, model = lab, term = names(b), n = sum(m),
      OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), row.names = NULL)
  }
}
reg <- do.call(rbind, out); reg <- reg[reg$term != "(Intercept)", ]
write.csv(reg, file.path(OUT_DIR, "regression_or.csv"), row.names = FALSE)
print(subset(reg, term %in% c("female", "quintile_5", "insured")), digits = 3)
