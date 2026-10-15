# =============================================================================
# Project 2, Stage 1: risk measurement
# FTSE 100 + S&P 500, weekly log returns 2000-2025
# AR(0)-GARCH(1,1) skew-t marginals + copula -> one-week portfolio VaR and ES,
# backtest 2005-2025 and out-of-sample check on 2026
#
# Run from the repository root:  Rscript code/stage1_analysis.R
# Outputs: printed tables, figures/*.png and results/stage1_results.rds
# (the report reports/stage1_risk_measurement_report.Rmd reads the .rds file)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages and settings
# -----------------------------------------------------------------------------
# Run from the repository root; the data files are read from the data/ folder.
DATA_DIR <- "data"
FIG_DIR  <- "figures"
RES_DIR  <- "results"
if (!file.exists(file.path(DATA_DIR, "FTSE100.csv")))
  stop("Data not found. Run this script from the repository root folder.")
dir.create(FIG_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(fGarch)      # GARCH fitting, psstd / qsstd
  library(VineCopula)  # copula selection, simulation, tail dependence
  library(ADGofTest)   # Anderson-Darling test
})

# Draw a plot on the screen and save the same plot as a PNG in figures/
show_and_save <- function(file, draw, width = 9, height = 6) {
  draw()
  png(file.path(FIG_DIR, file), width = width, height = height, units = "in",
      res = 150)
  draw()
  invisible(dev.off())
}

N_SIM <- 100000   # number of simulated scenarios for the new VaR figures
SEED  <- 333

# -----------------------------------------------------------------------------
# 1. Data: read, clean, compute log returns, align the two indices by week
# -----------------------------------------------------------------------------
# (Download code from the original script, run once to create the csv files)
# library(quantmod)
# getSymbols("^FTSE", src = "yahoo", from = "2000-01-01", to = "2025-12-31",
#            periodicity = "weekly")
# getSymbols("^GSPC", src = "yahoo", from = "2000-01-01", to = "2025-12-31",
#            periodicity = "weekly")
# write.csv(cbind(Date = as.Date(index(FTSE)), data.frame(coredata(FTSE))),
#           file.path(DATA_DIR, "FTSE100.csv"), row.names = FALSE)
# write.csv(cbind(Date = as.Date(index(GSPC)), data.frame(coredata(GSPC))),
#           file.path(DATA_DIR, "SNP500.csv"), row.names = FALSE)

load_index <- function(file, close_col) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df <- na.omit(unique(df))
  df$Date <- as.Date(trimws(df$Date))
  df <- df[order(df$Date), ]
  data.frame(Date   = df$Date[-1],
             Week   = format(df$Date[-1], "%G-%V"),   # ISO year-week
             Return = diff(log(df[[close_col]])))
}

ftse <- load_index(file.path(DATA_DIR, "FTSE100.csv"), "FTSE.Close")
snp  <- load_index(file.path(DATA_DIR, "SNP500.csv"),  "GSPC.Close")

# FIX: the original paired the two series by row position. Yahoo stamps some
# weeks on Friday and others on Saturday, so we merge on the ISO week instead.
dat <- merge(ftse, snp, by = "Week", suffixes = c("_ftse", "_snp"))
dat <- dat[order(dat$Date_ftse), ]
dat <- dat[is.finite(dat$Return_ftse) & is.finite(dat$Return_snp), ]

ret_ftse <- dat$Return_ftse
ret_snp  <- dat$Return_snp
dates    <- dat$Date_ftse
cat("Weeks in the aligned sample:", nrow(dat),
    "(", format(min(dates)), "to", format(max(dates)), ")\n")

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(dates, ret_ftse, type = "l", main = "FTSE 100 weekly log returns",
     xlab = "", ylab = "Log return")
plot(dates, ret_snp, type = "l", main = "S&P 500 weekly log returns",
     xlab = "", ylab = "Log return")
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# 2. Helper functions for the diagnostics
# -----------------------------------------------------------------------------
# Ljung-Box p-value; fitdf removes the degrees of freedom used by the AR terms
lb_p <- function(x, lag, fitdf = 0) {
  Box.test(x, lag = lag, type = "Ljung-Box", fitdf = fitdf)$p.value
}

# Engle's ARCH-LM test: regress z^2 on its own lags, T * R^2 ~ chi^2(lags)
arch_lm_p <- function(z, lags = 10) {
  X <- embed(z^2, lags + 1)
  r2 <- summary(lm(X[, 1] ~ X[, -1]))$r.squared
  pchisq(nrow(X) * r2, df = lags, lower.tail = FALSE)
}

# One row of the diagnostics table for a fitted fGarch model
diag_row <- function(model, p) {
  z <- residuals(model, standardize = TRUE)
  if (p > 0) z <- z[-(1:p)]   # FIX: fGarch sets the first p residuals to 0
  ic <- model@fit$ics
  c(AIC      = unname(ic["AIC"]),
    BIC      = unname(ic["BIC"]),
    LB_z_10  = lb_p(z, 10, p),
    LB_z_20  = lb_p(z, 20, p),
    LB_z2_10 = lb_p(z^2, 10),
    LB_z2_20 = lb_p(z^2, 20),
    ARCH_LM  = arch_lm_p(z, 10))
}

# -----------------------------------------------------------------------------
# 3. Mean equation: how many AR lags?
# -----------------------------------------------------------------------------
# One table per index replaces the 15 pages of repeated output in the report.
#  - LB_z   : Ljung-Box on standardised residuals   -> mean left unexplained?
#  - LB_z2  : Ljung-Box on squared std. residuals   -> volatility clusters left?
#  - ARCH_LM: same question as LB_z2, different test
p_max <- 5

select_ar <- function(ret, name) {
  rows <- lapply(0:p_max, function(p) {
    fml <- as.formula(paste0("~ arma(", p, ",0) + garch(1,1)"))
    m <- garchFit(fml, data = ret, cond.dist = "sstd", trace = FALSE)
    diag_row(m, p)
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- paste0("AR(", 0:p_max, ")")
  cat("\n=== ", name, ": AR(p)-GARCH(1,1) with skew-t, model selection ===\n")
  print(round(tab, 4))
  invisible(tab)
}

tab_ar_ftse <- select_ar(ret_ftse, "FTSE 100")
tab_ar_snp  <- select_ar(ret_snp,  "S&P 500")

# -----------------------------------------------------------------------------
# 4. Shock distribution for AR(0)-GARCH(1,1)
# -----------------------------------------------------------------------------
dist_vec    <- c("norm", "std", "sstd", "snorm", "ged", "sged")
dist_labels <- c("Normal", "Student-t", "Skew Student-t",
                 "Skew normal", "GED", "Skew GED")

select_dist <- function(ret, name) {
  fits <- lapply(dist_vec, function(d)
    garchFit(~ arma(0, 0) + garch(1, 1), data = ret, cond.dist = d,
             trace = FALSE))
  ic <- do.call(rbind, lapply(fits, function(m) m@fit$ics[c("AIC", "BIC")]))
  rownames(ic) <- dist_labels
  best <- which.min(ic[, "BIC"])
  cat("\n=== ", name, ": shock distribution ===\n")
  print(round(ic, 4))
  cat("Chosen (min BIC):", dist_labels[best], "\n")
  dist_ic_tabs[[name]] <<- ic
  fits[[best]]
}
dist_ic_tabs <- list()

model_ftse <- select_dist(ret_ftse, "FTSE 100")
model_snp  <- select_dist(ret_snp,  "S&P 500")

# -----------------------------------------------------------------------------
# 5. Final marginal models
# -----------------------------------------------------------------------------
summarise_model <- function(m) {
  cf <- coef(m)
  persistence <- cf["alpha1"] + cf["beta1"]
  c(cf,
    persistence    = unname(persistence),
    half_life_wks  = unname(log(0.5) / log(persistence)),
    long_run_vol   = unname(sqrt(cf["omega"] / (1 - persistence))))
}
final_tab <- rbind(FTSE = summarise_model(model_ftse),
                   SNP  = summarise_model(model_snp))
cat("\n=== Final AR(0)-GARCH(1,1) skew-t models ===\n")
print(signif(final_tab, 4))
# skew < 1: large falls more common than equally large rises
# shape   : fat tails (lower = fatter; very large = normal)
# half-life: weeks for a volatility shock to halve

# -----------------------------------------------------------------------------
# 6. PIT: residuals -> uniforms, and check they really are uniform
# -----------------------------------------------------------------------------
# The PIT gives exact uniforms only if the FITTED distribution is the true one,
# so these tests check the whole marginal model (GARCH filter + skew-t).
#  - KS weights the whole distribution, so mostly the centre
#  - AD weights the tails most heavily, which is where VaR lives
# Both tests treat the GARCH parameters as known. Because they were estimated
# from the same data, the p-values are approximate (usually somewhat too high).
pit <- function(m) {
  z <- residuals(m, standardize = TRUE)
  psstd(z, mean = 0, sd = 1, nu = coef(m)["shape"], xi = coef(m)["skew"])
}
u_ftse <- pit(model_ftse)
u_snp  <- pit(model_snp)

uniform_tests <- function(u) {
  # FIX: the earlier version used KScorrect::LcKS(u, "punif"), which corrects
  # for estimating the min and max of a uniform from the data. That is not the
  # parameter uncertainty we have here, and it is simulation-based, so its
  # p-value changed from run to run. The standard KS test against U(0,1) is
  # used instead (deterministic; approximate for the reason given above).
  c(KS_p           = ks.test(u, "punif")$p.value,
    AD_p           = unname(ad.test(u, punif)$p.value),
    share_below_5  = mean(u < 0.05),   # should be about 0.05
    share_below_1  = mean(u < 0.01))   # should be about 0.01
}
cat("\n=== Uniformity of PIT values ===\n")
print(round(rbind(FTSE = uniform_tests(u_ftse), SNP = uniform_tests(u_snp)), 4))

par(mfrow = c(1, 2))
hist(u_ftse, breaks = 10, main = "FTSE 100 PIT values", xlab = "u",
     col = "lightblue")
abline(h = length(u_ftse) / 10, lty = 2)
hist(u_snp, breaks = 10, main = "S&P 500 PIT values", xlab = "u",
     col = "lightblue")
abline(h = length(u_snp) / 10, lty = 2)
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# 7. Copula: selection, Kendall's tau, tail dependence
# -----------------------------------------------------------------------------
plot(u_ftse, u_snp, pch = 16, cex = 0.5, col = "blue",
     main = "FTSE 100 vs S&P 500 PIT values", xlab = "FTSE 100", ylab = "S&P 500")

cop_best <- BiCopSelect(u_ftse, u_snp, familyset = NA, selectioncrit = "AIC",
                        indeptest = TRUE, level = 0.05, se = TRUE)
cat("\nSelected copula:\n"); print(cop_best)

# For comparison: copulas with no / symmetric crash dependence, fitted to the
# same data
cop_gauss <- BiCopEst(u_ftse, u_snp, family = 1)
cop_t     <- BiCopEst(u_ftse, u_snp, family = 2)
cop_indep <- BiCop(family = 0)

copula_row <- function(cop, aic = cop$AIC) {
  td <- if (cop$family == 0) list(lower = 0, upper = 0)
        else BiCopPar2TailDep(cop$family, cop$par, cop$par2)
  c(tau          = if (cop$family == 0) 0 else
                     BiCopPar2Tau(cop$family, cop$par, cop$par2),
    lambda_lower = td$lower,   # P(both in worst weeks)  -> crash dependence
    lambda_upper = td$upper,   # P(both in best weeks)   -> boom dependence
    AIC          = if (is.null(aic)) 0 else aic)
}
cop_tab <- rbind(Independence  = copula_row(cop_indep, 0),
                 Gaussian      = copula_row(cop_gauss),
                 Student_t     = copula_row(cop_t),
                 Selected      = copula_row(cop_best))
rownames(cop_tab)[4] <- cop_best$familyname
cat("\n=== Copula comparison (same PIT data) ===\n")
print(round(cop_tab, 3))
cat("Empirical Kendall's tau of raw returns:",
    round(cor(ret_ftse, ret_snp, method = "kendall"), 3), "\n")

# Model-free check of crash dependence: among the q-share of worst weeks of the
# FTSE 100, how many were also among the q-share of worst weeks of the S&P 500?
# Uses the same GARCH-standardised residuals the copula is fitted to (so
# "worst weeks" = worst shocks relative to the volatility at the time), but no
# copula. The model value is P(V < q | U < q) = C(q, q) / q for each copula.
co_exceed <- function(q) {
  rf_ <- rank(u_ftse) / (length(u_ftse) + 1)
  rs_ <- rank(u_snp)  / (length(u_snp)  + 1)
  n_q <- sum(rf_ < q)
  c(q = q, worst_weeks = n_q,
    shared_worst = sum(rf_ < q & rs_ < q),
    share_worst  = sum(rf_ < q & rs_ < q) / n_q,
    model_BB1    = BiCopCDF(q, q, cop_best$family, cop_best$par, cop_best$par2) / q,
    model_t      = BiCopCDF(q, q, cop_t$family, cop_t$par, cop_t$par2) / q,
    model_gauss  = BiCopCDF(q, q, cop_gauss$family, cop_gauss$par) / q,
    shared_best  = sum(rf_ > 1 - q & rs_ > 1 - q))
}
tail_tab <- t(sapply(c(0.05, 0.02, 0.01), co_exceed))
cat("\n=== Shared extreme weeks: data vs copulas ===\n")
print(round(tail_tab, 3))

# -----------------------------------------------------------------------------
# 8. FIX: one-week-ahead VaR, as in the lecturer's hint
# -----------------------------------------------------------------------------
# sigma^2(T+1) = omega + alpha * u(T)^2 + beta * sigma^2(T), from REAL data.
# Every scenario shares this volatility; only the copula shocks differ.
next_week_sigma2 <- function(m) {
  cf <- coef(m)
  u_T      <- tail(residuals(m, standardize = FALSE), 1)
  sigma2_T <- tail(volatility(m), 1)^2
  unname(cf["omega"] + cf["alpha1"] * u_T^2 + cf["beta1"] * sigma2_T)
}
s2_next <- c(FTSE = next_week_sigma2(model_ftse),
             SNP  = next_week_sigma2(model_snp))
cat("\nNext-week volatility vs long-run volatility:\n")
print(round(rbind(next_week = sqrt(s2_next),
                  long_run  = final_tab[, "long_run_vol"]), 4))

# Simulate n versions of next week's 50/50 portfolio return
simulate_next_week <- function(cop, n = N_SIM, seed = SEED) {
  set.seed(seed)
  U <- BiCopSim(n, cop$family, cop$par, cop$par2)
  cf_f <- coef(model_ftse); cf_s <- coef(model_snp)
  eps_f <- qsstd(U[, 1], 0, 1, nu = cf_f["shape"], xi = cf_f["skew"])
  eps_s <- qsstd(U[, 2], 0, 1, nu = cf_s["shape"], xi = cf_s["skew"])
  r_f <- cf_f["mu"] + sqrt(s2_next["FTSE"]) * eps_f
  r_s <- cf_s["mu"] + sqrt(s2_next["SNP"])  * eps_s
  0.5 * (exp(r_f) - 1) + 0.5 * (exp(r_s) - 1)   # net returns, then weight
}

# VaR = loss exceeded with probability 1 - level
# ES  = average loss in the scenarios worse than the VaR
risk_measures <- function(port) {
  v95 <- -quantile(port, 0.05, names = FALSE)
  v99 <- -quantile(port, 0.01, names = FALSE)
  c(VaR_95 = v95, ES_95 = -mean(port[port <= -v95]),
    VaR_99 = v99, ES_99 = -mean(port[port <= -v99]))
}

# -----------------------------------------------------------------------------
# 9. Original method, kept and relabelled: unconditional VaR ("a typical week")
# -----------------------------------------------------------------------------
# One GARCH path of 10,000 weeks. After a few dozen steps it forgets its
# starting point, so it describes weeks in general, not next week.
# Reproduces the coursework numbers exactly (same seed, same order of calls).
simulate_long_path <- function(cop, N = 10000, seed = 333) {
  set.seed(seed)
  U <- BiCopSim(N, family = cop$family, cop$par, cop$par2)
  path <- function(m, eps) {
    cf <- coef(m)
    s2 <- numeric(N); u <- numeric(N)
    s2[1] <- tail(volatility(m), 1)^2
    u[1]  <- sqrt(s2[1]) * eps[1]
    for (t in 2:N) {
      s2[t] <- cf["omega"] + cf["alpha1"] * u[t - 1]^2 + cf["beta1"] * s2[t - 1]
      u[t]  <- sqrt(s2[t]) * eps[t]
    }
    cf["mu"] + u
  }
  cf_f <- coef(model_ftse); cf_s <- coef(model_snp)
  r_f <- path(model_ftse, qsstd(U[, 1], nu = cf_f["shape"], xi = cf_f["skew"]))
  r_s <- path(model_snp,  qsstd(U[, 2], nu = cf_s["shape"], xi = cf_s["skew"]))
  list(port = 0.5 * (exp(r_f) - 1) + 0.5 * (exp(r_s) - 1), r_f = r_f, r_s = r_s)
}
long <- simulate_long_path(cop_best)
port_uncond <- long$port

# In-sample check against the real weeks
port_hist <- 0.5 * (exp(ret_ftse) - 1) + 0.5 * (exp(ret_snp) - 1)

# -----------------------------------------------------------------------------
# 10. Results
# -----------------------------------------------------------------------------
port_next <- simulate_next_week(cop_best)
res_main <- rbind(
  "Next week (conditional, corrected)" = risk_measures(port_next),
  "Typical week (unconditional, original method)" = risk_measures(port_uncond),
  "Real weeks 2000-2025 (historical)" = risk_measures(port_hist))
cat("\n=== One-week risk of the 50/50 portfolio (loss as % of value) ===\n")
print(round(100 * res_main, 2))

uncond_var <- risk_measures(port_uncond)
cat("\nShare of real weeks with a loss beyond the unconditional VaR:\n",
    " 95% VaR:", round(100 * mean(port_hist < -uncond_var["VaR_95"]), 2),
    "% of weeks (target 5%)\n",
    " 99% VaR:", round(100 * mean(port_hist < -uncond_var["VaR_99"]), 2),
    "% of weeks (target 1%)\n")

# Same marginals, same next-week volatility, only the copula changes
res_cop <- rbind(
  Independence = risk_measures(simulate_next_week(cop_indep)),
  Gaussian     = risk_measures(simulate_next_week(cop_gauss)),
  Student_t    = risk_measures(simulate_next_week(cop_t)),
  Selected     = risk_measures(port_next))
rownames(res_cop)[4] <- cop_best$familyname
cat("\n=== Next-week risk by copula (loss as % of value) ===\n")
print(round(100 * res_cop, 2))

# Probability that both indices fall into their worst 1% of weeks together
p_joint_1pct <- c(Selected = BiCopCDF(0.01, 0.01, cop_best$family, cop_best$par,
                                      cop_best$par2),
                  Gaussian = BiCopCDF(0.01, 0.01, 1, cop_gauss$par))
cat("P(both in worst 1% of weeks): selected", round(100 * p_joint_1pct[1], 2),
    "%, Gaussian", round(100 * p_joint_1pct[2], 2), "%\n")

# Sensitivity: does the AR(0) choice matter? BIC prefers AR(4) for the FTSE 100
# and AR(3) for the S&P 500 (section 3). Recompute next week's VaR with those
# AR orders (each index's forecast mean includes its AR terms).
var_with_ar <- function(p) {
  fml <- as.formula(paste0("~ arma(", p, ",0) + garch(1,1)"))
  fit_one <- function(r) {
    m  <- garchFit(fml, data = r, cond.dist = "sstd", trace = FALSE)
    cf <- coef(m)
    mean_next <- cf["mu"] + if (p > 0) sum(cf[paste0("ar", 1:p)] * rev(tail(r, p))) else 0
    s2_next   <- cf["omega"] + cf["alpha1"] * tail(residuals(m), 1)^2 +
                 cf["beta1"] * tail(volatility(m), 1)^2
    list(cf = cf, m = mean_next, s = sqrt(s2_next),
         u = psstd(residuals(m, standardize = TRUE)[-(seq_len(p))], 0, 1,
                   cf["shape"], cf["skew"]))
  }
  f <- fit_one(ret_ftse); s <- fit_one(ret_snp)
  cop_p <- BiCopSelect(f$u, s$u, familyset = NA, selectioncrit = "AIC")
  set.seed(SEED)
  U <- BiCopSim(N_SIM, cop_p$family, cop_p$par, cop_p$par2)
  port <- 0.5 * (exp(f$m + f$s * qsstd(U[, 1], 0, 1, f$cf["shape"], f$cf["skew"])) - 1) +
          0.5 * (exp(s$m + s$s * qsstd(U[, 2], 0, 1, s$cf["shape"], s$cf["skew"])) - 1)
  risk_measures(port)
}
ar_sens <- rbind("AR(0) (used)"        = risk_measures(port_next),
                 "AR(1)"               = var_with_ar(1),
                 "AR(3)"               = var_with_ar(3),
                 "AR(4)"               = var_with_ar(4))
cat("\n=== Next-week VaR/ES under different AR orders (%) ===\n")
print(round(100 * ar_sens, 2))

# Simulated vs historical densities (check of the marginal models)
plot_densities <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  for (k in 1:2) {
    hist_r <- list(ret_ftse, ret_snp)[[k]]
    sim_r  <- list(long$r_f, long$r_s)[[k]]
    plot(density(hist_r), col = "blue", lwd = 2,
         main = c("FTSE 100", "S&P 500")[k], xlab = "Weekly log return")
    lines(density(sim_r), col = "red", lwd = 2)
    legend("topleft", c("Historical", "Simulated"), col = c("blue", "red"),
           lwd = 2, bty = "n")
  }
  par(op)
}
show_and_save("simulated_vs_historical_density.png", plot_densities, 9, 4.5)

# =============================================================================
# 11. Backtest of the one-week VaR, 2005-2025
# =============================================================================
# At each week t, using ONLY data up to t: re-estimate the models every 13
# weeks, update the volatility every week, forecast next week's VaR, and
# compare it with the real portfolio return in week t+1.
START_TEST  <- as.Date("2005-01-01")  # first week that is forecast
REFIT_EVERY <- 13                     # full re-estimation every quarter
BT_N_SIM    <- 10000                  # simulated scenarios per forecast
set.seed(333)

ret_f <- ret_ftse; ret_s <- ret_snp   # aligned weekly log returns from section 1
port  <- port_hist                    # real 50/50 portfolio net returns

# Fit marginals + copula on data 1..t, simulate copula shocks, and compute the
# unconditional VaR with the original long-path method
fit_all <- function(t) {
  mf <- garchFit(~ garch(1, 1), data = ret_f[1:t], cond.dist = "sstd",
                 trace = FALSE)
  ms <- garchFit(~ garch(1, 1), data = ret_s[1:t], cond.dist = "sstd",
                 trace = FALSE)
  cf <- coef(mf); cs <- coef(ms)
  uf <- psstd(residuals(mf, standardize = TRUE), 0, 1, cf["shape"], cf["skew"])
  us <- psstd(residuals(ms, standardize = TRUE), 0, 1, cs["shape"], cs["skew"])
  cop <- BiCopSelect(uf, us, familyset = NA, selectioncrit = "AIC")
  U <- BiCopSim(BT_N_SIM, cop$family, cop$par, cop$par2)
  eps_f <- qsstd(U[, 1], 0, 1, cf["shape"], cf["skew"])
  eps_s <- qsstd(U[, 2], 0, 1, cs["shape"], cs["skew"])
  path <- function(c, eps, s2_start) {
    s2 <- s2_start; u <- numeric(BT_N_SIM)
    for (k in 1:BT_N_SIM) {
      if (k > 1) s2 <- c["omega"] + c["alpha1"] * u[k - 1]^2 + c["beta1"] * s2
      u[k] <- sqrt(s2) * eps[k]
    }
    c["mu"] + u
  }
  rp <- 0.5 * (exp(path(cf, eps_f, tail(volatility(mf), 1)^2)) - 1) +
        0.5 * (exp(path(cs, eps_s, tail(volatility(ms), 1)^2)) - 1)
  list(cf = cf, cs = cs, eps_f = eps_f, eps_s = eps_s, family = cop$familyname,
       s2_f = tail(volatility(mf), 1)^2, s2_s = tail(volatility(ms), 1)^2,
       uncond = -quantile(rp, c(0.05, 0.01), names = FALSE))
}

# Walk forward one week at a time
test_idx <- which(dates >= START_TEST)
bt <- data.frame(Date = dates[test_idx], actual = port[test_idx],
                 cVaR95 = NA, cVaR99 = NA, uVaR95 = NA, uVaR99 = NA, copula = NA)
for (k in seq_along(test_idx)) {
  t <- test_idx[k] - 1                         # last week we are allowed to see
  if ((k - 1) %% REFIT_EVERY == 0) {
    m <- fit_all(t)
    s2f <- m$s2_f; s2s <- m$s2_s               # variance of week t from the fit
  } else {                                     # update with last week's return
    s2f <- m$cf["omega"] + m$cf["alpha1"] * (ret_f[t - 1] - m$cf["mu"])^2 +
      m$cf["beta1"] * s2f
    s2s <- m$cs["omega"] + m$cs["alpha1"] * (ret_s[t - 1] - m$cs["mu"])^2 +
      m$cs["beta1"] * s2s
  }
  # forecast variance for week t+1, as in the hint
  s2f_next <- m$cf["omega"] + m$cf["alpha1"] * (ret_f[t] - m$cf["mu"])^2 +
    m$cf["beta1"] * s2f
  s2s_next <- m$cs["omega"] + m$cs["alpha1"] * (ret_s[t] - m$cs["mu"])^2 +
    m$cs["beta1"] * s2s
  rp <- 0.5 * (exp(m$cf["mu"] + sqrt(s2f_next) * m$eps_f) - 1) +
        0.5 * (exp(m$cs["mu"] + sqrt(s2s_next) * m$eps_s) - 1)
  bt[k, c("cVaR95", "cVaR99")] <- -quantile(rp, c(0.05, 0.01), names = FALSE)
  bt[k, c("uVaR95", "uVaR99")] <- m$uncond
  bt$copula[k] <- m$family
}

# =============================================================================
# 12. Backtest evaluation
# =============================================================================
# Kupiec test: is the NUMBER of hits consistent with 5% / 1%?
kupiec_test <- function(hit, p) {
  n <- length(hit); x <- sum(hit); ph <- x / n
  ll0 <- (n - x) * log(1 - p) + x * log(p)
  ll1 <- (n - x) * log(1 - ph) + ifelse(x > 0, x * log(ph), 0)
  c(weeks = n, hits = x, expected = n * p, hit_rate = ph,
    kupiec_p = pchisq(-2 * (ll0 - ll1), 1, lower.tail = FALSE))
}
# Christoffersen test: are the hits SPREAD OUT, or does a hit make another hit
# next week more likely? Compares P(hit | hit last week) with P(hit | no hit).
christoffersen_test <- function(hit) {
  h <- as.integer(hit); a <- h[-length(h)]; b <- h[-1]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1)
  n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  L <- function(x, k) ifelse(k > 0, k * log(x), 0)
  p01 <- n01 / (n00 + n01); p11 <- n11 / (n10 + n11)
  p <- (n01 + n11) / (n00 + n01 + n10 + n11)
  ll0 <- L(1 - p, n00 + n10) + L(p, n01 + n11)
  ll1 <- L(1 - p01, n00) + L(p01, n01) + L(1 - p11, n10) + L(p11, n11)
  c(hit_after_no_hit = p01, hit_after_hit = p11, back_to_back = n11,
    christoffersen_p = pchisq(-2 * (ll0 - ll1), 1, lower.tail = FALSE))
}
bt_tab <- NULL
for (v in c("cVaR95", "cVaR99", "uVaR95", "uVaR99")) {
  p <- ifelse(grepl("95", v), 0.05, 0.01)
  hit <- bt$actual < -bt[[v]]
  bt_tab <- rbind(bt_tab, c(kupiec_test(hit, p), christoffersen_test(hit)))
}
rownames(bt_tab) <- c("Conditional 95%", "Conditional 99%",
                      "Unconditional 95%", "Unconditional 99%")
cat("\n=== Backtest 2005-2025 ===\n")
print(round(bt_tab, 4))

# Hits per year: shows WHERE the hits happen
bt$year <- format(bt$Date, "%Y")
hits_year <- aggregate(cbind(cond95 = actual < -cVaR95, uncond95 = actual < -uVaR95,
                             cond99 = actual < -cVaR99, uncond99 = actual < -uVaR99)
                       ~ year, data = bt, FUN = sum)
# A conservative rule: use the larger of the two VaRs each week, so the VaR
# never drops below its long-run level in calm times but still rises in crises.
max_tab <- t(sapply(c("95", "99"), function(lv) {
  p <- ifelse(lv == "95", 0.05, 0.01)
  m <- pmax(bt[[paste0("cVaR", lv)]], bt[[paste0("uVaR", lv)]])
  hit <- bt$actual < -m
  c(kupiec_test(hit, p), christoffersen_test(hit), avg_VaR = mean(m),
    hits_2008 = sum(hit[bt$year == "2008"]))
}))
rownames(max_tab) <- c("max(cond, uncond) 95%", "max(cond, uncond) 99%")
cat("\n=== Conservative rule: larger of the two VaRs ===\n")
print(round(max_tab, 4))

cat("\n=== Hits per year ===\n")
print(hits_year)

# Copula chosen at each re-fit, using only past data
cat("\n=== Copula selected in the backtest (weeks) ===\n")
print(table(bt$copula))
cop_switch <- rle(bt$copula)
cat("Copula choice over time (first forecast week of each run):\n")
print(data.frame(copula = cop_switch$values,
                 from   = bt$Date[cumsum(cop_switch$lengths) - cop_switch$lengths + 1],
                 weeks  = cop_switch$lengths))

# Plot: real returns, both VaR lines, hits marked
plot_backtest <- function() {
op <- par(mfrow = c(2, 1), mar = c(3, 4.5, 2.5, 1))
for (lv in c("95", "99")) {
  cv <- bt[[paste0("cVaR", lv)]]; uv <- bt[[paste0("uVaR", lv)]]
  plot(bt$Date, 100 * bt$actual, type = "h", col = "grey60",
       ylim = c(min(100 * bt$actual, -100 * cv) - 1, 13), xlab = "",
       ylab = "Weekly return (%)",
       main = paste0(lv, "% VaR: conditional (red) vs unconditional (blue)"))
  lines(bt$Date, -100 * cv, col = "red", lwd = 1.5)
  lines(bt$Date, -100 * uv, col = "blue", lwd = 1.5)
  hc <- bt$actual < -cv; hu <- bt$actual < -uv
  points(bt$Date[hc], 100 * bt$actual[hc], pch = 16, col = "red", cex = 0.8)
  points(bt$Date[hu], 100 * bt$actual[hu], pch = 1, col = "blue", cex = 1.3)
}
par(op)
}
show_and_save("backtest_var_2005_2025.png", plot_backtest, 10, 7)

# =============================================================================
# 13. Out-of-sample test: January - September 2026
# =============================================================================
# The model is FROZEN as estimated on 2000-2025 (sections 4-7): no parameter
# is re-estimated. Each week of 2026, next week's volatility is updated with
# the newly observed return (exactly as in the hint), the VaR is forecast, and
# then compared with what really happened.
# Data: FTSE100_2026.csv and SNP500_2026.csv (same download, end date 2026-10-01)
new <- merge(load_index(file.path(DATA_DIR, "FTSE100_2026.csv"), "FTSE.Close"),
             load_index(file.path(DATA_DIR, "SNP500_2026.csv"),  "GSPC.Close"),
             by = "Week", suffixes = c("_ftse", "_snp"))
new <- new[order(new$Date_ftse), ]
new <- new[is.finite(new$Return_ftse) & is.finite(new$Return_snp), ]

# Note: the original csv files were downloaded before 31 December 2025, so
# their last week was incomplete. The new files contain the full week, so the
# last 2025 return is taken from the new files.
last_2025 <- max(which(new$Date_ftse < as.Date("2026-01-01")))
oos_idx   <- which(new$Date_ftse >= as.Date("2026-01-01"))

cf_f <- coef(model_ftse); cf_s <- coef(model_snp)
# Copula shocks from the frozen survival BB1 copula
set.seed(SEED)
U_oos  <- BiCopSim(N_SIM, cop_best$family, cop_best$par, cop_best$par2)
eps_f  <- qsstd(U_oos[, 1], 0, 1, cf_f["shape"], cf_f["skew"])
eps_s  <- qsstd(U_oos[, 2], 0, 1, cf_s["shape"], cf_s["skew"])

# Variance of the last 2025 week (uses data up to the week before, unchanged)
s2f <- tail(volatility(model_ftse), 1)^2
s2s <- tail(volatility(model_snp),  1)^2
uncond <- risk_measures(port_uncond)          # unconditional VaR, section 9

oos <- data.frame(Date = new$Date_ftse[oos_idx],
                  actual = 0.5 * (exp(new$Return_ftse[oos_idx]) - 1) +
                           0.5 * (exp(new$Return_snp[oos_idx])  - 1),
                  cVaR95 = NA, cVaR99 = NA)
t <- last_2025
for (k in seq_along(oos_idx)) {
  # forecast variance for the coming week from last week's REAL return
  s2f <- cf_f["omega"] + cf_f["alpha1"] * (new$Return_ftse[t] - cf_f["mu"])^2 +
         cf_f["beta1"] * s2f
  s2s <- cf_s["omega"] + cf_s["alpha1"] * (new$Return_snp[t] - cf_s["mu"])^2 +
         cf_s["beta1"] * s2s
  rp <- 0.5 * (exp(cf_f["mu"] + sqrt(s2f) * eps_f) - 1) +
        0.5 * (exp(cf_s["mu"] + sqrt(s2s) * eps_s) - 1)
  oos$cVaR95[k] <- -quantile(rp, 0.05, names = FALSE)
  oos$cVaR99[k] <- -quantile(rp, 0.01, names = FALSE)
  t <- t + 1
}
oos$uVaR95 <- uncond["VaR_95"]; oos$uVaR99 <- uncond["VaR_99"]

oos_tab <- NULL
for (v in c("cVaR95", "cVaR99", "uVaR95", "uVaR99")) {
  p <- ifelse(grepl("95", v), 0.05, 0.01)
  oos_tab <- rbind(oos_tab, kupiec_test(oos$actual < -oos[[v]], p))
}
rownames(oos_tab) <- rownames(bt_tab)
cat("\n=== Out-of-sample 2026 (model frozen at end of 2025) ===\n")
print(round(oos_tab, 4))

cat("\nWorst weeks of 2026 and what the model forecast for them (%):\n")
worst <- oos[order(oos$actual), ][1:5, ]
print(data.frame(Date = worst$Date, actual = round(100 * worst$actual, 2),
                 cVaR95 = round(100 * worst$cVaR95, 2),
                 cVaR99 = round(100 * worst$cVaR99, 2)))

plot_oos <- function() {
op <- par(mar = c(3, 4.5, 2.5, 1))
plot(oos$Date, 100 * oos$actual, type = "h", lwd = 3, col = "grey50",
     ylim = c(min(100 * oos$actual, -100 * oos$cVaR99) - 3,
              max(100 * oos$actual) + 1),
     xlab = "", ylab = "Weekly return (%)",
     main = "2026 out-of-sample: real returns vs forecast VaR")
lines(oos$Date, -100 * oos$cVaR95, col = "red", lwd = 2)
lines(oos$Date, -100 * oos$cVaR99, col = "red", lwd = 2, lty = 2)
abline(h = -100 * uncond[c("VaR_95", "VaR_99")], col = "blue", lty = c(1, 2))
legend("bottom", c("Conditional 95%", "Conditional 99%",
                   "Unconditional 95%", "Unconditional 99%"),
       col = c("red", "red", "blue", "blue"), lty = c(1, 2, 1, 2), bty = "n",
       cex = 0.8, ncol = 4)
par(op)
}
show_and_save("out_of_sample_2026.png", plot_oos, 9, 4.5)

# Power of the 2026 check: chance of seeing this few hits if the VaR were right
cat("\nIf the VaR were exactly right, P(0 hits of the 99% VaR in", nrow(oos),
    "weeks) =", round(0.99^nrow(oos), 2), "and P(at most 1 hit of the 95% VaR) =",
    round(pbinom(1, nrow(oos), 0.05), 2), "\n")

# =============================================================================
# 14. Save everything the report needs
# =============================================================================
saveRDS(list(
  n_weeks = nrow(dat), first_date = min(dates), last_date = max(dates),
  tab_ar_ftse = tab_ar_ftse, tab_ar_snp = tab_ar_snp, dist_ic = dist_ic_tabs,
  model_ftse = model_ftse, model_snp = model_snp, final_tab = final_tab,
  pit = rbind(FTSE = uniform_tests(u_ftse), SNP = uniform_tests(u_snp)),
  u_ftse = u_ftse, u_snp = u_snp, cop_best = cop_best, cop_tab = cop_tab,
  tail_tab = tail_tab, s2_next = s2_next, res_main = res_main, res_cop = res_cop,
  p_joint_1pct = p_joint_1pct, ar_sens = ar_sens, uncond_var = uncond_var,
  hist_share = c(p95 = mean(port_hist < -uncond_var["VaR_95"]),
                 p99 = mean(port_hist < -uncond_var["VaR_99"])),
  bt = bt, bt_tab = bt_tab, hits_year = hits_year, max_tab = max_tab, cop_switch = cop_switch,
  oos = oos, oos_tab = oos_tab, ret_ftse = ret_ftse, ret_snp = ret_snp,
  dates = dates, long_r_f = long$r_f, long_r_s = long$r_s),
  file.path(RES_DIR, "stage1_results.rds"))
cat("\nSaved results/stage1_results.rds\n")
sessionInfo()
