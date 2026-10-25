# =====================================================================
# 04_figures.R - the seven report figures (ggplot2 + patchwork) -> output/figs/*.png
# install.packages(c("ggplot2", "patchwork"))
# =====================================================================
source("00_config.R")
library(ggplot2); library(patchwork)
r <- readRDS(file.path(OUT_DIR, "results.rds")); res <- r$res; ex <- r$extras
COL <- c(Rajasthan = "#2a78d6", Sikkim = "#eb6834")
est <- function(st, k) res[res$state == st & res$key == k, ]
get <- function(keys, labels, scale = 1, floor0 = FALSE) {       # tidy table of est/lo/hi for both states
  do.call(rbind, lapply(names(COL), function(st) {
    d <- do.call(rbind, lapply(keys, function(k) est(st, k)))
    d$label <- factor(labels, levels = rev(labels)); d$state <- st
    d[c("est", "lo", "hi")] <- d[c("est", "lo", "hi")] * scale
    if (floor0) d$lo <- pmax(d$lo, 0)
    d
  }))
}
th <- theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e3e2de", linewidth = 0.3), legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 9.5), axis.title = element_text(colour = "#52514e"))
sv <- function(p, f, w, h) ggsave(file.path(OUT_DIR, "figs", f), p, width = w, height = h, dpi = 200, bg = "white")
pd <- position_dodge(width = 0.5)

# ---- Figure 1: age-sex pyramids
ages <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+")
py <- do.call(rbind, lapply(names(COL), function(st) {
  p <- ex[[st]]$pyramid
  rbind(data.frame(state = st, sex = "Male", age = factor(ages, levels = ages), pct = -p[, "M"]),
        data.frame(state = st, sex = "Female", age = factor(ages, levels = ages), pct = p[, "F"]))
}))
f1 <- ggplot(py, aes(pct, age, fill = state, alpha = sex)) + geom_col(width = 0.78) +
  facet_wrap(~state) + scale_fill_manual(values = COL, guide = "none") + scale_alpha_manual(values = c(Male = 0.95, Female = 0.55)) +
  scale_x_continuous(labels = abs) + geom_vline(xintercept = 0, colour = "#8a8982") +
  labs(x = "% of state population (males left, females right)", y = NULL) + th + theme(legend.position = "bottom")
sv(f1, "f1_pyramid.png", 6.6, 3.7)

# ---- Figure 2: age pattern of reported pregnancies
xl <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")
f2d <- do.call(rbind, lapply(names(COL), function(st) { d <- do.call(rbind, lapply(paste0("asfr_", 4:10), function(k) est(st, k)))
  data.frame(state = st, age = factor(xl, levels = xl), est = d$est * 1000, lo = pmax(d$lo, 0) * 1000, hi = d$hi * 1000) }))
f2 <- ggplot(f2d, aes(age, est, colour = state, fill = state, group = state)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) + scale_colour_manual(values = COL) + scale_fill_manual(values = COL, guide = "none") +
  labs(x = "Woman's age (years)", y = "Reported pregnancies per 1,000 women",
       title = sprintf("TFI (sum of age-specific rates x 5):  Rajasthan %.2f   Sikkim %.2f", est("Rajasthan", "tfi")$est, est("Sikkim", "tfi")$est)) +
  coord_cartesian(ylim = c(0, 200)) + th
sv(f2, "f2_fertility.png", 6, 2.9)

# ---- Figure 3: female/male ratios
k3 <- c("chronic_ratioFM", "ail15_ratioFM", "hosp_ratioFM", "hospnr_ratioFM", "cov_insured_ratioFM", "vacc04_ratioFM", "gpi_lit15", "gpi_sec25")
l3 <- c("Chronic ailment (age-standardised)", "Any ailment, last 15 days (age-std.)", "Hospitalised, all ages (age-std.)", "Hospitalised, excluding ages 15-49",
        "Covered by a health scheme", "Any vaccine, ages 0-4", "Literacy, ages 15+ (parity index)", "Secondary+ education, ages 25+ (parity index)")
d3 <- get(k3, l3); d3$lo <- pmax(d3$lo, 0); d3$hi <- pmin(d3$hi, 2.9)
f3 <- ggplot(d3, aes(est, label, colour = state)) + geom_vline(xintercept = 1, linetype = "dashed", colour = "#6b6a64") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = pd, linewidth = 0.7) + geom_point(position = pd, size = 2) +
  scale_colour_manual(values = COL) + coord_cartesian(xlim = c(0, 3)) +
  labs(x = "Female / male ratio (1.0 = parity); bars = 95% CI", y = NULL) + th
sv(f3, "f3_gender.png", 6.6, 3.6)

# ---- Figure 4: Lorenz curves and MPCE levels
lz <- do.call(rbind, lapply(names(COL), function(st) cbind(state = st, ex[[st]]$lorenz)))
g <- sprintf("Gini: Rajasthan %.3f, Sikkim %.3f", est("Rajasthan", "gini")$est, est("Sikkim", "gini")$est)
f4a <- ggplot(lz, aes(x, y, colour = state)) + geom_abline(linetype = "dashed", colour = "#8a8982") + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = COL) + coord_equal() + labs(x = "Cumulative share of persons (poorest first)", y = "Cumulative share of consumption", title = "Lorenz curves of MPCE", subtitle = g) + th + theme(legend.position = "none")
d4 <- get(c("mpce_p10", "mpce_median", "mpce_mean", "mpce_p90"), c("P10", "Median", "Mean", "P90"))
f4b <- ggplot(d4, aes(est, label, colour = state)) + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = pd) + geom_point(position = pd, size = 2) +
  geom_text(aes(x = hi, label = format(round(est), big.mark = ",")), position = pd, hjust = -0.2, size = 2.6, colour = "#52514e") +
  scale_colour_manual(values = COL) + coord_cartesian(xlim = c(0, 12000)) + labs(x = "MPCE (Rs per person per month); bars = 95% CI", y = NULL, title = "Level of MPCE") + th
sv(f4a + f4b + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "f4_inequality.png", 6.6, 3.1)

# ---- Figure 5: gradients by MPCE quintile
pan <- list(c("Chronic ailment, age-std.", "q%d_chronic_std"), c("Hospitalised, age-std.", "q%d_hosp_std"), c("Any health scheme", "q%d_cov_insured"))
d5 <- do.call(rbind, lapply(pan, function(p) do.call(rbind, lapply(names(COL), function(st) {
  d <- do.call(rbind, lapply(1:5, function(q) est(st, sprintf(p[2], q))))
  data.frame(panel = p[1], state = st, q = 1:5, est = d$est * 100, lo = pmax(d$lo, 0) * 100, hi = d$hi * 100) }))))
d5$panel <- factor(d5$panel, levels = sapply(pan, `[`, 1))
f5 <- ggplot(d5, aes(q, est, colour = state)) + geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(0.25), alpha = 0.8) +
  geom_line(position = position_dodge(0.25)) + geom_point(position = position_dodge(0.25), size = 1.8) +
  facet_wrap(~panel, scales = "free_y") + scale_x_continuous(breaks = 1:5, labels = c("Q1\npoorest", "Q2", "Q3", "Q4", "Q5\nrichest")) +
  expand_limits(y = 0) + scale_colour_manual(values = COL) + labs(x = NULL, y = "% of persons") + th
sv(f5, "f5_quintile.png", 6.8, 2.7)

# ---- Figure 6: coverage by scheme type
k6 <- c("cov_stscheme", "cov_pmjay", "cov_emp", "cov_private", "cov_none")
l6 <- c("State health insurance scheme", "AB-PMJAY", "Employer / govt. employee schemes", "Private insurance only", "Not covered")
d6 <- get(k6, l6, 100)
f6 <- ggplot(d6, aes(est, label, fill = state)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", est)), position = position_dodge(0.75), hjust = -0.15, size = 2.7) +
  scale_fill_manual(values = COL) + scale_x_continuous(limits = c(0, 100)) + labs(x = "% of persons (one scheme recorded per person)", y = NULL) + th
sv(f6, "f6_schemes.png", 6.4, 2.9)

# ---- Figure 7: concentration indices
k7 <- c("ci_insured", "ci_pmjay", "ci_stscheme", "ci_emp", "ci_chronic_std", "ci_hosp_std")
l7 <- c("Any health scheme", "AB-PMJAY", "State health insurance scheme", "Employer / govt. employee schemes", "Chronic ailment (age-sex std.)", "Hospitalisation (age-sex std.)")
d7 <- get(k7, l7)
f7 <- ggplot(d7, aes(est, label, colour = state)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "#6b6a64") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = pd) + geom_point(position = pd, size = 2) +
  scale_colour_manual(values = COL) + coord_cartesian(xlim = c(-0.7, 0.8)) +
  labs(x = "Concentration index (> 0 richer, < 0 poorer); bars = 95% CI", y = NULL) + th
sv(f7, "f7_concentration.png", 6.6, 3.0)
cat("Seven figures written to", file.path(OUT_DIR, "figs"), "\n")
