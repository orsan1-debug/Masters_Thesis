# sbw_wc.R -----------------------------------------------------------------
# Reproduce the "Variance" rows of Table 4 in Wang & Zubizarreta (2020) on
# the Wong & Chan design: stable balancing weights (Zubizarreta 2015; sbw)
# with exact balance and approximate balance tuned by their Algorithm 1
# using the paper's Wong & Chan settings (B = 10 bootstrap samples of size
# n / 10, first and second moments, K = 20).
# sbw::sbw(bal_alg = TRUE) bootstraps the full arm size bal_sam times and
# the subsample size cannot be changed, so the tuning loop is written out
# here around fixed-tolerance sbw fits. The weights themselves are sbw's.
# Solver is quadprog (exact active set). osqp inside sbw is too loose: a
# 4e-5 per-weight error moved model A estimates by 0.2. quadprog rejects
# delta = 0, so the exact arm uses the smallest tolerance it accepts
# (delta_exact, in SD units).

library(sbw)
source(paste0(dir, "dgp_wc.R"))

#' sbw weights for one tolerance delta (in SD units of the target group).
#' Returns the n-vector of weights in data order. For "ate" each arm sums to
#' one; for "att" the controls sum to one and each treated unit gets 1 / n1.
fit_sbw <- function(dat, bal_cov, delta, estimand, solver = "quadprog") {
  fit <- sbw(dat, ind = "W",
             bal = list(bal_cov = bal_cov, bal_alg = FALSE, bal_tol = delta,
                        bal_std = "target"),
             sol = list(sol_nam = solver),
             par = list(par_est = estimand), mes = FALSE)
  fit$dat_weights$sbw_weights
}

#' Algorithm 1 criterion: mean over B subsamples (size frac * arm size) of
#' the mean standardised absolute imbalance of the weighted arm(s) against
#' the target means, summed over arms as sbw does for the ATE.
cstat <- function(w, B_mat, W, arms, target, s, B = 10, frac = 0.1) {
  total <- 0
  for (a in arms) {
    idx <- which(W == a)
    m   <- round(frac * length(idx))
    for (b in seq_len(B)) {
      i  <- sample(idx, m, replace = TRUE)
      wm <- colSums(w[i] * B_mat[i, , drop = FALSE]) / sum(w[i])
      total <- total + mean(abs(wm - target) / s)
    }
  }
  total / B
}

#' One replication: exact and Algorithm-1-tuned sbw estimates of the ATE
#' and ATT for every outcome column of a dgp_wc() draw d. grid is in SD
#' units and stays below K^(-1/2) = 0.22 (Figure 4 points).
estimate_wc_sbw <- function(d, grid = c(0.001, 0.002, 0.005, 0.01, 0.02,
                                        0.05, 0.1, 0.2),
                            delta_exact = 1e-4, solver = "quadprog") {
  W <- d$W
  X <- cbind(d$X, d$X^2)
  colnames(X) <- c(paste0("x", 1:10), paste0("x", 1:10, "sq"))
  dat <- data.frame(W = W, X)
  s   <- apply(X, 2, sd)
  
  tau_row <- function(w, estimand, balance, delta) {
    if (is.character(w)) {
      return(data.frame(estimand, balance, delta, outcome = colnames(d$Y),
                        tau_hat = NA_real_, error = w))
    }
    stopifnot(abs(sum(w * W) - 1) < 1e-6, abs(sum(w * (1 - W)) - 1) < 1e-6)
    data.frame(estimand, balance, delta, outcome = colnames(d$Y),
               tau_hat = colSums(w * W * d$Y) - colSums(w * (1 - W) * d$Y),
               error = NA_character_)
  }
  
  rows <- lapply(c("ate", "att"), \(estimand) {
    arms   <- if (estimand == "ate") c(0, 1) else 0
    target <- if (estimand == "ate") colMeans(X) else colMeans(X[W == 1, ])
    fit    <- \(delta) tryCatch(
      fit_sbw(dat, colnames(X), delta, estimand, solver),
      error = \(e) conditionMessage(e))
    w_grid <- lapply(grid, fit)
    cs <- vapply(w_grid, \(w) if (is.character(w)) Inf else
      cstat(w, X, W, arms, target, s), numeric(1))
    k  <- which.min(cs)
    rbind(tau_row(fit(delta_exact), estimand, "exact", delta_exact),
          tau_row(w_grid[[k]], estimand, "approx", grid[k]))
  })
  do.call(rbind, rows)
}