# =============================================================================
# Project 2, Stage 2: position sizing
# Part 1: simulate many possible futures of the 50/50 portfolio
# =============================================================================
# The model from Stage 1 (AR(0)-GARCH(1,1) skew-t marginals + survival BB1
# copula, estimated on 2000-2025) is used as the "world". Each path starts from
# the real market state at the end of 2025 and runs forward week by week.
#
# Run from the repository root:  Rscript code/stage2_position_sizing.R
# Outputs: printed tables, figures/*.png and results/stage2_results.rds
# (the report reports/stage2_position_sizing_report.Rmd reads the .rds file)
DATA_DIR <- "data"
FIG_DIR  <- "figures"
RES_DIR  <- "results"
if (!file.exists(file.path(DATA_DIR, "FTSE100.csv")))
  stop("Data not found. Run this script from the repository root folder.")
dir.create(FIG_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)
show_and_save <- function(file, draw, width = 9, height = 6) {
  draw()
  png(file.path(FIG_DIR, file), width = width, height = height, units = "in",
      res = 150)
  draw()
  invisible(dev.off())
}

suppressPackageStartupMessages({
  library(fGarch)
  library(VineCopula)
})

# -----------------------------------------------------------------------------
# 1. Fit the Stage 1 model (same data and choices as the report)
# -----------------------------------------------------------------------------
load_index <- function(file, close_col) {
  df <- na.omit(unique(read.csv(file, stringsAsFactors = FALSE)))
  df$Date <- as.Date(trimws(df$Date)); df <- df[order(df$Date), ]
  data.frame(Date = df$Date[-1], Week = format(df$Date[-1], "%G-%V"),
             Return = diff(log(df[[close_col]])))
}
dat <- merge(load_index(file.path(DATA_DIR, "FTSE100.csv"), "FTSE.Close"),
             load_index(file.path(DATA_DIR, "SNP500.csv"),  "GSPC.Close"),
             by = "Week", suffixes = c("_ftse", "_snp"))
dat <- dat[order(dat$Date_ftse), ]
dat <- dat[is.finite(dat$Return_ftse) & is.finite(dat$Return_snp), ]  # as Stage 1
ret_ftse <- dat$Return_ftse; ret_snp <- dat$Return_snp

model_ftse <- garchFit(~ garch(1, 1), data = ret_ftse, cond.dist = "sstd", trace = FALSE)
model_snp  <- garchFit(~ garch(1, 1), data = ret_snp,  cond.dist = "sstd", trace = FALSE)
cf_f <- coef(model_ftse); cf_s <- coef(model_snp)
u_ftse <- as.numeric(psstd(as.numeric(residuals(model_ftse, standardize = TRUE)), 0, 1, cf_f["shape"], cf_f["skew"]))
u_snp  <- as.numeric(psstd(as.numeric(residuals(model_snp,  standardize = TRUE)), 0, 1, cf_s["shape"], cf_s["skew"]))
cop <- BiCopSelect(u_ftse, u_snp, familyset = NA, selectioncrit = "AIC")
cat("Copula:", cop$familyname, "\n")

# The world's average return (drift). GARCH maximum likelihood gives each week a
# weight that falls as its volatility rises, so the volatile weeks, when most
# of the big losses happened, count less, and the estimated mean comes out much
# higher than the simple average of the returns. A world with the GARCH mean
# grows about twice as fast as the real indices did in 2000-2025, so the world
# uses the sample mean of the weekly log returns instead.
mu_garch  <- c(FTSE = unname(cf_f["mu"]), SNP = unname(cf_s["mu"]))
mu_sample <- c(FTSE = mean(ret_ftse), SNP = mean(ret_snp))
cat("Weekly mean log return, GARCH vs sample:\n")
print(signif(rbind(GARCH = mu_garch, sample = mu_sample), 3))
cf_f["mu"] <- mu_sample["FTSE"]; cf_s["mu"] <- mu_sample["SNP"]

# -----------------------------------------------------------------------------
# 2. Path simulator
# -----------------------------------------------------------------------------
# n_paths possible futures, each n_weeks long. Inside each path the GARCH
# variance is chained week by week (volatility clusters form), and the two
# indices are linked by the copula. All paths share the same starting point:
# the real last shock and variance at the end of 2025.
# Returns a matrix of weekly portfolio net returns (rows = paths, cols = weeks),
# with the portfolio rebalanced to 50/50 every week.
simulate_paths <- function(n_paths, n_weeks, seed = 2026) {
  set.seed(seed)
  U <- as.matrix(BiCopSim(n_paths * n_weeks, cop$family, cop$par, cop$par2))
  eps_f <- matrix(as.numeric(qsstd(U[, 1], 0, 1, cf_f["shape"], cf_f["skew"])), nrow = n_paths)
  eps_s <- matrix(as.numeric(qsstd(U[, 2], 0, 1, cf_s["shape"], cf_s["skew"])), nrow = n_paths)

  # starting state: real last shock u_T and variance sigma^2_T
  u_f  <- rep(tail(as.numeric(residuals(model_ftse)), 1), n_paths)
  u_s  <- rep(tail(as.numeric(residuals(model_snp)),  1), n_paths)
  s2_f <- rep(tail(as.numeric(volatility(model_ftse)), 1)^2, n_paths)
  s2_s <- rep(tail(as.numeric(volatility(model_snp)),  1)^2, n_paths)

  R <- matrix(NA_real_, n_paths, n_weeks)
  for (w in 1:n_weeks) {               # loop over weeks, all paths at once
    s2_f <- cf_f["omega"] + cf_f["alpha1"] * u_f^2 + cf_f["beta1"] * s2_f
    s2_s <- cf_s["omega"] + cf_s["alpha1"] * u_s^2 + cf_s["beta1"] * s2_s
    u_f <- sqrt(s2_f) * eps_f[, w]
    u_s <- sqrt(s2_s) * eps_s[, w]
    R[, w] <- 0.5 * (exp(cf_f["mu"] + u_f) - 1) + 0.5 * (exp(cf_s["mu"] + u_s) - 1)
  }
  R
}

R1 <- simulate_paths(n_paths = 10000, n_weeks = 52)

# -----------------------------------------------------------------------------
# 3. Checks: does the simulated world behave like the real one?
# -----------------------------------------------------------------------------
# (a) Week 1 must reproduce the Stage 1 conditional VaR (about 4.0% at 99%)
cat("\nWeek 1 VaR 95% / 99% (%):",
    round(-100 * quantile(R1[, 1], c(0.05, 0.01)), 2), "\n")

# (b) Volatility clustering inside the paths: big weeks follow big weeks
abs_now  <- abs(R1[, -52]); abs_next <- abs(R1[, -1])
cat("Correlation of |return| this week vs next week:",
    round(cor(as.vector(abs_now), as.vector(abs_next)), 3), "\n")

# (c) One-year outcomes at 1x (no leverage) vs the real 2000-2025 history
wealth   <- t(apply(1 + R1, 1, cumprod))              # wealth path, start = 1
max_dd   <- function(w) { w <- c(1, w); max(1 - w / cummax(w)) }
sim_year <- wealth[, 52] - 1
sim_dd   <- apply(wealth, 1, max_dd)

port_hist <- 0.5 * (exp(ret_ftse) - 1) + 0.5 * (exp(ret_snp) - 1)
hist_year <- sapply(1:(length(port_hist) - 51), function(i)
  prod(1 + port_hist[i:(i + 51)]) - 1)               # every 52-week window
hist_dd   <- sapply(1:(length(port_hist) - 51), function(i)
  max_dd(cumprod(1 + port_hist[i:(i + 51)])))

summ <- function(x) round(100 * c(q05 = quantile(x, 0.05, names = FALSE),
                                  median = median(x),
                                  q95 = quantile(x, 0.95, names = FALSE)), 1)
cat("\nOne-year return (%)           5%  median  95%\n")
cat("  Simulated (from end-2025):", summ(sim_year), "\n")
cat("  Real 52-week windows:     ", summ(hist_year), "\n")
cat("\nWorst drawdown within a year (%)\n")
cat("  Simulated:", summ(sim_dd), "\n")
cat("  Real:     ", summ(hist_dd), "\n")
cat("\nP(lose more than 20% within the year at 1x):",
    round(100 * mean(sim_dd > 0.20), 1), "% simulated vs",
    round(100 * mean(hist_dd > 0.20), 1), "% real windows\n")

# The same comparison by calendar year (non-overlapping, so a binomial
# confidence interval is valid), and for simulated years that have forgotten
# the calm end-2025 start (100-week burn-in).
yr <- format(dat$Date_ftse, "%Y")
cal_dd <- tapply(port_hist, yr, function(x) max_dd(cumprod(1 + x)))
k_cal  <- sum(cal_dd > 0.20)
ci_cal <- binom.test(k_cal, length(cal_dd))$conf.int
R_burn <- simulate_paths(n_paths = 10000, n_weeks = 152, seed = 7)[, 101:152]
dd_burn <- apply(t(apply(1 + R_burn, 1, cumprod)), 1, max_dd)
cat("Calendar years with a fall of 20%+ from the peak:", k_cal, "of",
    length(cal_dd), sprintf("(%.1f%%, 95%% CI %.1f%% to %.1f%%)\n",
                            100 * k_cal / length(cal_dd), 100 * ci_cal[1], 100 * ci_cal[2]))
cat("Simulated years after a 100-week burn-in:", round(100 * mean(dd_burn > 0.20), 1), "%\n")
growth_fold <- exp(c(sample = sum(ret_snp), garch = unname(mu_garch["SNP"]) * length(ret_snp)))
cat("100 put into the S&P 500 in 2000 became: actually", round(100 * growth_fold[1]),
    "; with the GARCH mean it would have become", round(100 * growth_fold[2]), "\n")

# =============================================================================
# Part 2: the Kelly leverage in this world
# =============================================================================
# An investor holds L times their wealth in the 50/50 portfolio (L > 1 means
# borrowing), rebalanced every week:  W(t+1) = W(t) * (1 + L * R(t)).
# If a week's loss is bigger than 1/L, wealth hits zero: the investor is wiped
# out. Simplification: price returns exclude dividends and borrowing costs are
# set to zero; over 2000-2025 the two roughly offset each other.
#
# Long paths (200 weeks after a 100-week burn-in, so the calm end-2025 start is
# forgotten) describe the long run.
R_long <- simulate_paths(n_paths = 10000, n_weeks = 300)[, 101:300]
n_years <- ncol(R_long) / 52

# How often does the model produce very bad weeks?
cat("\nHow often the model produces a very bad week:\n")
for (th in c(0.15, 0.196, 0.25, 0.30)) {
  k <- sum(R_long < -th)
  cat(sprintf("  worse than -%4.1f%%: once every %5.0f years\n",
              100 * th, (length(R_long) / 52) / max(k, 1)))
}
cat("  (real data: worst week -19.6%, October 2008, once in 26 years)\n")

# Kelly leverage with the standard formula: mean / variance of weekly returns.
# This is the usual approximation to the leverage that maximises expected log
# growth. Strictly, the model's shocks are unbounded, so for any L > 1 there is
# a tiny chance of a week that loses more than 1/L (more than everything).
# log(0) = -infinity, so the exact expected log growth is -infinity for all
# L > 1. The approximation ignores this very rare event (frequency below).
L_kelly <- mean(R_long) / var(as.vector(R_long))
cat("\nKelly leverage (mean / variance):", round(L_kelly, 2), "\n")
cat("Chance per week of losing more than 1/L at L =", round(L_kelly, 2), ":",
    signif(mean(R_long < -1 / L_kelly), 2), "\n")

# Sensitivity to the two omitted items. Dividends add to the return and the
# borrowing rate is a cost, so what matters is (dividend yield - borrowing
# rate). Kelly leverage if that net carry were +1% or -1% a year:
carry <- c(-0.01, 0, 0.01)
L_carry <- (mean(R_long) + carry / 52) / var(as.vector(R_long))
names(L_carry) <- paste0("carry ", 100 * carry, "%/yr")
cat("Kelly leverage under different net carry:\n"); print(round(L_carry, 2))

# Check it directly: typical (median) growth of wealth for different leverages
leverage_outcome <- function(L, R) {
  x <- pmax(1 + L * R, 0)                         # wealth cannot go below zero
  wiped <- apply(x == 0, 1, any)
  log_w <- rowSums(log(pmax(x, 1e-300)))          # log of final wealth
  c(L = L,
    median_growth_pct = 100 * median(log_w) / n_years,
    wiped_out_pct     = 100 * mean(wiped))
}
L_values <- c(0.5, 1, 1.5, 2, round(L_kelly, 2), 3, 4, 5)
lev_tab <- t(sapply(L_values, leverage_outcome, R = R_long))
cat("\nLeverage vs typical growth per year and chance of being wiped out (",
    round(n_years, 1), " years):\n", sep = "")
print(round(lev_tab, 2))

# An investor who trusts the GARCH means instead of the sample means
cf_f_keep <- cf_f; cf_s_keep <- cf_s
cf_f["mu"] <- mu_garch["FTSE"]; cf_s["mu"] <- mu_garch["SNP"]
R_garch <- simulate_paths(n_paths = 10000, n_weeks = 300)[, 101:300]
cf_f <- cf_f_keep; cf_s <- cf_s_keep
L_garch <- mean(R_garch) / var(as.vector(R_garch))
cat("\nKelly leverage computed with the GARCH means:", round(L_garch, 2),
    "(", round(L_garch / L_kelly, 1), "x the true Kelly )\n")
print(round(leverage_outcome(L_garch, R_long), 2))

# Growth versus leverage
L_grid <- seq(0, 6, by = 0.1)
g_grid <- sapply(L_grid, function(L) leverage_outcome(L, R_long)["median_growth_pct"])
cat("Leverage with the highest typical (median) growth on the grid:",
    L_grid[which.max(g_grid)], "\n")
plot_growth <- function() {
  op <- par(mar = c(4.5, 4.5, 2.5, 1))
  plot(L_grid, g_grid, type = "l", lwd = 2,
       xlab = "Leverage L (position per 1 of own money)",
       ylab = "Typical (median) growth of wealth (% per year)",
       main = "Growth versus leverage in the simulated world")
  abline(h = 0, lty = 2)
  abline(v = L_kelly, col = "red", lty = 2); abline(v = L_garch, col = "blue", lty = 2)
  legend("bottomleft", c("Kelly (true average)", "Kelly computed with GARCH means"),
         col = c("red", "blue"), lty = 2, bty = "n")
  par(op)
}
show_and_save("growth_vs_leverage.png", plot_growth, 8, 5)

# =============================================================================
# Part 3: investors who KNOW the odds
# =============================================================================
# Each investor knows the true average and variance of the weekly returns
# (so they can compute the true Kelly leverage), but not the future returns.
# This is the best case; in Part 4 the investors have to estimate the odds.
investors <- c("No leverage" = 1,
               "Half Kelly"  = 0.5 * L_kelly,
               "Full Kelly"  = L_kelly,
               "2x Kelly"    = 2 * L_kelly)
cat("\nInvestors (position per 1 of own money):\n"); print(round(investors, 2))

# ---- (a) Weekly risk, from the simulated world (the validated part) ---------
# With weekly rebalancing, the investor's equity return in a week is L * R.
weekly_risk <- function(L, R) {
  eq <- pmax(L * as.vector(R), -1)                  # cannot lose more than 100%
  v99 <- -quantile(eq, 0.01, names = FALSE)
  every <- function(th) {               # Inf = never happened in the simulation
    k <- sum(eq <= -th); if (k == 0) Inf else (length(eq) / 52) / k
  }
  c("VaR 99% (%)"              = 100 * v99,
    "ES 99% (%)"               = 100 * -mean(eq[eq <= -v99]),
    "Week losing 25%+: every (yrs)" = every(0.25),
    "Week losing 50%+: every (yrs)" = every(0.50),
    "Wiped out in a week: every (yrs)" = every(1))
}
tab_week <- t(sapply(investors, weekly_risk, R = R_long))
cat("\n=== (a) One-week risk to the investor's own money (simulated long run) ===\n")
print(round(tab_week, 1))

# ---- (b) Real history: replay 2000-2025 with each leverage -------------------
dates_all <- dat$Date_ftse
port_real <- 0.5 * (exp(ret_ftse) - 1) + 0.5 * (exp(ret_snp) - 1)
replay <- function(L) {
  eq <- pmax(1 + L * port_real, 0)
  w  <- cumprod(eq)
  # Yahoo dates each weekly bar on the Friday/Saturday BEFORE its trading week,
  # so the labels 2020-02-21 to 2020-03-15 are the four crash weeks of
  # 24 February - 20 March 2020. (FIX: the earlier window, 15 Feb - 21 Mar,
  # also included the quiet week before and the rebound week after.)
  covid <- dates_all >= as.Date("2020-02-21") & dates_all <= as.Date("2020-03-15")
  c("Worst real week (%)"        = 100 * (min(eq) - 1),
    "Covid, 24 Feb-20 Mar 2020 (%)" = 100 * (prod(eq[covid]) - 1),
    "Worst drawdown (%)"         = 100 * max(1 - w / cummax(w)),
    "1 became"                   = tail(w, 1),
    "Growth per year (%)"        = 100 * (tail(w, 1)^(52 / length(w)) - 1))
}
tab_hist <- t(sapply(investors, replay))
cat("\n=== (b) Replaying the real weeks of 2000-2025 ===\n")
print(round(tab_hist, 2))

# ---- (c) The year ahead, simulated from end-2025 (a lower bound) -------------
year_risk <- function(L, R) {
  w <- t(apply(pmax(1 + L * R, 0), 1, cumprod))
  dd <- apply(w, 1, function(x) max(1 - c(1, x) / cummax(c(1, x))))
  c("P(drawdown 50%+ in year) (%)" = 100 * mean(dd >= 0.5),
    "P(drawdown 25%+ in year) (%)" = 100 * mean(dd >= 0.25),
    "Median 1-year return (%)"     = 100 * (median(w[, ncol(w)]) - 1))
}
tab_year <- t(sapply(investors, year_risk, R = R1))
cat("\n=== (c) The year ahead from end-2025, simulated (optimistic) ===\n")
print(round(tab_year, 2))

# ---- (d) Half Kelly vs full Kelly ---------------------------------------------
# Long-run growth from the simulated world, risk from all three views
half_vs_full <- sapply(c("Half Kelly" = 0.5 * L_kelly, "Full Kelly" = L_kelly), function(L) {
  c("Growth per year, simulated (%)" = unname(leverage_outcome(L, R_long)["median_growth_pct"]),
    "Growth per year, real 2000-2025 (%)" = unname(replay(L)["Growth per year (%)"]),
    "99% one-week loss (%)"          = unname(weekly_risk(L, R_long)["VaR 99% (%)"]),
    "Worst real drawdown (%)"        = unname(replay(L)["Worst drawdown (%)"]),
    "P(drawdown 50%+ next year) (%)" = unname(year_risk(L, R1)["P(drawdown 50%+ in year) (%)"]))
})
half_vs_full <- cbind(half_vs_full,
                      "Half as % of full" = 100 * half_vs_full[, 1] / half_vs_full[, 2])
cat("\n=== (d) Half Kelly vs full Kelly ===\n")
print(round(half_vs_full, 2))

# =============================================================================
# Part 4: investors who ESTIMATE the odds
# =============================================================================
# Real investors do not know the true average return; they estimate it from
# past data. With 26 years of weekly data, the uncertainty of that estimate is
# its standard error, sd / sqrt(n). Each Kelly investor is given two guesses:
#   overestimate  = true average + 1 standard error
#   underestimate = true average - 1 standard error
# About 1 investor in 6 guesses too high by at least this much, and 1 in 6 too
# low. The variance is assumed known: it is estimated far more precisely.
m_true <- mean(R_long)
v_true <- var(as.vector(R_long))
se_mean <- sd(port_real) / sqrt(length(port_real))
guesses <- c(Over = m_true + se_mean, True = m_true, Under = m_true - se_mean)
cat("\nAverage weekly return: true", signif(m_true, 3), " standard error",
    signif(se_mean, 3), "\n")
cat("Real 2000-2025 average", signif(mean(port_real), 3), " t-statistic",
    round(mean(port_real) / se_mean, 2), " 95% CI", signif(mean(port_real) - 1.96 * se_mean, 2),
    "to", signif(mean(port_real) + 1.96 * se_mean, 2),
    "(includes zero: the data cannot rule out a zero average return)\n")
cat("Guesses:\n"); print(signif(guesses, 3))

# Position each investor takes with each guess
part4 <- expand.grid(Odds = names(guesses), Rule = c("Full Kelly", "Half Kelly"),
                     stringsAsFactors = FALSE)
part4$Leverage <- ifelse(part4$Rule == "Full Kelly", 1, 0.5) *
                  guesses[part4$Odds] / v_true
part4 <- rbind(data.frame(Odds = "-", Rule = "No leverage", Leverage = 1), part4)

res4 <- t(sapply(part4$Leverage, function(L) c(
  "Position per 100"          = 100 * L,
  "Growth/yr, simulated (%)"  = unname(leverage_outcome(L, R_long)["median_growth_pct"]),
  "VaR 99% 1 week (%)"        = unname(weekly_risk(L, R_long)["VaR 99% (%)"]),
  "25%+ week: every (yrs)"    = unname(weekly_risk(L, R_long)["Week losing 25%+: every (yrs)"]),
  "P(fall 50%+ next yr) (%)"  = unname(year_risk(L, R1)["P(drawdown 50%+ in year) (%)"]),
  "Real max drawdown (%)"     = -unname(replay(L)["Worst drawdown (%)"]),
  "Real: 1 became"            = unname(replay(L)["1 became"]))))
rownames(res4) <- ifelse(part4$Odds == "-", part4$Rule,
                         paste0(part4$Rule, ", ", tolower(part4$Odds), "estimate"))
rownames(res4) <- sub("trueestimate", "true odds", rownames(res4))
cat("\n=== Part 4: Kelly investors with over- and underestimated odds ===\n")
print(round(res4, 2))

# Margin calls: with a 25% maintenance margin and weekly rebalancing to L, a
# call is triggered by a single-week fall of x = (1 - 0.25 L) / (0.75 L).
# (The replays above assume no forced selling; this shows how often that
# assumption would have been broken in the real data.)
margin_tab <- t(sapply(c("Full Kelly, true" = L_kelly,
                         "Full Kelly, over" = part4$Leverage[part4$Rule == "Full Kelly" & part4$Odds == "Over"],
                         "Half Kelly, over" = part4$Leverage[part4$Rule == "Half Kelly" & part4$Odds == "Over"],
                         "Half Kelly, true" = 0.5 * L_kelly), function(L) {
  x <- (1 - 0.25 * L) / (0.75 * L)
  c(L = L, fall_triggering_call_pct = 100 * x, real_weeks_with_call = sum(port_real < -x))
}))
cat("\n=== Margin calls in the real weeks (25% maintenance margin) ===\n")
print(round(margin_tab, 2))

# =============================================================================
# Part 5: estimation error, averaged over all the guesses an investor could make
# =============================================================================
# Part 4 looks at two particular errors (+/- 1 standard error). Here each of
# 4,000 investors gets their own estimate of the average, drawn from its
# sampling distribution N(true average, se^2). Each one sizes with a fraction c
# of their estimated Kelly leverage (no short selling: a negative estimate means
# no position). We then average their outcomes.
# Simple theory (quadratic approximation of growth, variance known): the
# fraction of estimated Kelly that maximises average growth is
#   c* = m^2 / (m^2 + se^2) = 1 / (1 + (se/m)^2).
c_star <- 1 / (1 + (se_mean / m_true)^2)
# Why the variance can be treated as known: its relative standard error,
# sqrt((kurtosis - 1) / n), is far smaller than that of the average (se / mean)
kurt_real <- mean((port_real - mean(port_real))^4) / var(port_real)^2
rel_se_var  <- sqrt((kurt_real - 1) / length(port_real))
rel_se_mean <- se_mean / mean(port_real)
cat("Relative standard error: average", round(100 * rel_se_mean), "%, variance",
    round(100 * rel_se_var), "%\n")
cat("\nGrowth-maximising fraction of Kelly under estimation error (theory):",
    round(c_star, 2), "\n")

# Outcomes on a grid of leverages, then interpolated for each investor
p_half_year <- function(L, R) {                 # fast: loop over 52 weeks
  w <- rep(1, nrow(R)); peak <- w; worst <- rep(0, nrow(R))
  for (j in seq_len(ncol(R))) {
    w <- pmax(w * (1 + L * R[, j]), 0); peak <- pmax(peak, w)
    worst <- pmax(worst, 1 - w / peak)
  }
  mean(worst >= 0.5)
}
L_fine  <- seq(0, 10, by = 0.1)
g_fine  <- sapply(L_fine, function(L) leverage_outcome(L, R_long)["median_growth_pct"])
ph_fine <- sapply(L_fine, p_half_year, R = R1)
set.seed(5)
mu_hat <- rnorm(4000, m_true, se_mean)
part5 <- t(sapply(c(0.25, 0.5, 0.75, 1), function(cf) {
  L <- pmin(cf * pmax(mu_hat, 0) / v_true, max(L_fine))
  g <- approx(L_fine, g_fine, L)$y; ph <- approx(L_fine, ph_fine, L)$y
  c(fraction_of_kelly       = cf,
    avg_growth_pct          = mean(g),
    avg_P_fall50_next_yr    = 100 * mean(ph),
    share_more_than_2xKelly = 100 * mean(L > 2 * L_kelly),
    share_P_fall50_over_10  = 100 * mean(ph > 0.10))
}))
cat("\n=== Part 5: averaged over 4,000 investors' estimation errors ===\n")
print(round(part5, 2))
# The same on a finer grid of fractions, to locate the best one
c_grid <- seq(0.3, 1.3, by = 0.05)
g_c <- sapply(c_grid, function(cf) {
  L <- pmin(cf * pmax(mu_hat, 0) / v_true, max(L_fine))
  mean(approx(L_fine, g_fine, L)$y)
})
cat("Fraction of Kelly with the highest average growth (simulated):",
    c_grid[which.max(g_c)], " theory:", round(c_star, 2), "\n")
cat("Share of investors whose estimate is negative (no position):",
    round(100 * mean(mu_hat <= 0), 1), "%\n")

saveRDS(list(
  mu_garch = mu_garch, mu_sample = mu_sample, week1_var = -quantile(R1[, 1], c(0.05, 0.01)),
  sim_year = sim_year, hist_year = hist_year, sim_dd = sim_dd, hist_dd = hist_dd,
  k_cal = k_cal, n_cal = length(cal_dd), ci_cal = ci_cal, dd_burn = dd_burn,
  growth_fold = growth_fold, bad_week_freq = sapply(c(0.15, 0.196, 0.25, 0.30), function(th)
    (length(R_long) / 52) / max(sum(R_long < -th), 1)),
  L_kelly = L_kelly, L_garch = L_garch, L_carry = L_carry, lev_tab = lev_tab,
  L_grid = L_grid, g_grid = g_grid, p_wipe_week = mean(R_long < -1 / L_kelly),
  investors = investors, tab_week = tab_week, tab_hist = tab_hist,
  tab_year = tab_year, half_vs_full = half_vs_full, m_true = m_true,
  real_mean = mean(port_real), se_mean = se_mean, res4 = res4,
  margin_tab = margin_tab, c_star = c_star, rel_se_var = rel_se_var, rel_se_mean = rel_se_mean, part5 = part5, c_grid = c_grid, g_c = g_c,
  share_negative = mean(mu_hat <= 0), n_weeks = length(port_real)),
  file.path(RES_DIR, "stage2_results.rds"))
cat("\nSaved results/stage2_results.rds\n")
sessionInfo()
