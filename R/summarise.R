# ============================================================================
# analysis_setup.R -- shared setup for the E1-E6 analysis scripts.
# ============================================================================

library(dplyr)
library(tidyr)      
library(ggplot2)
library(patchwork)  # (l | r) panels
library(gt)
# pin dplyr verbs against masking (MASS, stats, ...)
select <- dplyr::select
filter <- dplyr::filter

drop0 <- function(x) sub("^(-?)0\\.", "\\1.", x)   # ".05" not "0.05" on axes

# ---- 1. Load ---------------------------------------------------------------

design_cols <- c("n", "p", "s", "decay", "misspec", "covcor", "outcome", "overlap",
                 "signs", "strength", "decay_ps", "decay_out", "outcome_set",
                 "treat_prop")
outcome_cols <- c("linear", "quad1", "quad2", "exp")

load_sim <- function(path) {
  d <- tibble::as_tibble(data.table::fread(path))
  d <- mutate(d, across(any_of(outcome_cols), \(x) suppressWarnings(as.numeric(x))))
  if ("err" %in% names(d)) d$err <- coalesce(as.character(d$err), "")
  pivot_longer(d, any_of(outcome_cols), names_to = "outcome", values_to = "tau_hat")
}

# diagnostics are logged as rows of `estimator`; identify them by name
is_diag <- function(x)
  grepl("^(lam_|nnz_|smd[01]_|cvloss_|logcvloss_|trunc05|prev|emin|emax|nout)", x) |
  x %in% c("trunc05", "prev", "emin", "emax", "nout05", "nout01")


# ---- 2. Cell summaries -----------------------------------------------------

# Screen n_na / n_extreme before reading anything else; bias/empse/rmse
# propagate NA loudly by design.
summarise_sim <- function(data, by = c(design_cols, "estimator")) {
  by <- intersect(by, names(data))
  data |>
    filter(!is_diag(estimator)) |>
    group_by(across(all_of(by))) |>
    summarise(
      m          = dplyr::n(),
      bias       = mean(tau_hat),
      empse      = sd(tau_hat),
      rmse       = sqrt(mean(tau_hat^2)),
      mcse_bias  = empse / sqrt(m),
      mcse_empse = empse / sqrt(2 * (m - 1)),
      mcse_rmse  = sd(tau_hat^2) / (2 * rmse * sqrt(m)),
      n_extreme  = sum(abs(tau_hat) > 10, na.rm = TRUE),
      n_na       = sum(is.na(tau_hat)),
      .groups    = "drop"
    )
}

# ---- 9. Misspecification helpers -------------------------------------------

# delta = misspecified - correct at matching (n, outcome, overlap, estimator)
delta_ms <- function(summ, ests = est_main) {
  summ |>
    filter(estimator %in% ests) |>
    mutate(absbias = abs(bias)) |>
    select(n, outcome, overlap, estimator, misspec, rmse, absbias) |>
    pivot_wider(names_from = misspec, values_from = c(rmse, absbias)) |>
    mutate(d_rmse = rmse_TRUE - rmse_FALSE,
           d_bias = absbias_TRUE - absbias_FALSE)
}

# ---- 11. Diagnostics infrastructure ----------------------------------------

# Diagnostics arrive as estimator rows; diag_wide() lifts them to one row per
# rep (per outcome -- values repeat across outcomes). as_diag() accepts raw
# or already-wide input, so every helper below can be fed replicate data
# directly.
diag_wide <- function(data) {
  id <- intersect(c(design_cols, "sim"), names(data))
  d  <- filter(data, is_diag(estimator))
  if (!nrow(d)) stop("no diagnostic rows found")
  dw <- pivot_wider(d, id_cols = all_of(id), names_from = estimator,
                    values_from = tau_hat)
  err_tab <- data |> distinct(across(all_of(id)), err) |>
    mutate(err = as.integer(err != ""))
  left_join(dw, err_tab, by = id)
}

as_diag <- function(data)
  if ("lam_glmcv" %in% names(data)) data else diag_wide(data)

# med / q10 / q90 columns -> long, one row per stat x quantile set
pivot_q <- function(dg, prefix)
  pivot_longer(dg, matches(paste0("^", prefix, "_(med|q10|q90)$")),
               names_to = c("stat", ".value"),
               names_pattern = paste0("^(", prefix, ")_(med|q10|q90)$"))

# cell-level med / q10 / q90 of every logged quantity + failure counts
summarise_diag <- function(data, by = design_cols) {
  dw <- as_diag(data)
  by <- intersect(by, names(dw))
  dw |>
    group_by(across(all_of(by))) |>
    summarise(
      m         = dplyr::n(),
      n_err     = sum(err != 0, na.rm = TRUE),
      p_trunc05 = mean(trunc05, na.rm = TRUE),
      across(matches("^(lam|nnz|smd|cvloss|logcvloss|prev|emin|emax|nout)"),
             list(med = \(x) median(x, na.rm = TRUE),
                  q10 = \(x) quantile(x, .10, na.rm = TRUE, names = FALSE),
                  q90 = \(x) quantile(x, .90, na.rm = TRUE, names = FALSE))),
      .groups = "drop"
    )
}

# ---- 13. KKT / attained-imbalance check (balnet eq 10) ---------------------

# Bound: attained max|SMD| <= lambda per arm. Reports worst gap and share of
# reps exceeding it by > tol for whichever of {endpoint, .05, .10, CV} the
# file logged. Gaps <= 0 with slack shrinking in n are the fixed-lambda
# plateau mechanism. NB the driver logs a self-normalised SMD; small positive
# gaps may be that definitional mismatch, not solver failure.
kkt_check <- function(data, tol = 1e-6, by = design_cols) {
  dw  <- as_diag(data)
  by  <- intersect(by, names(dw))
  has <- \(...) all(c(...) %in% names(dw))
  if (has("smd1_0", "lam_end1"))
    dw$gap0  <- pmax(dw$smd1_0 - dw$lam_end1, dw$smd0_0 - dw$lam_end0)
  if (has("smd1_05")) dw$gap05 <- pmax(dw$smd1_05, dw$smd0_05) - 0.05
  if (has("smd1_10")) dw$gap10 <- pmax(dw$smd1_10, dw$smd0_10) - 0.10
  if (has("smd1_cv", "lam_balcv1"))
    dw$gapcv <- pmax(dw$smd1_cv - dw$lam_balcv1, dw$smd0_cv - dw$lam_balcv0)
  dw |>
    group_by(across(all_of(by))) |>
    summarise(m = dplyr::n(),
              across(starts_with("gap"),
                     list(max      = \(x) max(x, na.rm = TRUE),
                          pct_viol = \(x) 100 * mean(x > tol, na.rm = TRUE))),
              .groups = "drop")
}

# ---- 15. CV vs fixed lambda (E8: where and why CV loses) -------------------

# Per-rep join of selected estimates with that rep's lambdas. lam_cv averages
# the two arms because the bias mechanism sums both (the 2-lambda KKT term).
rep_cv <- function(data, oc, ests) {
  id <- intersect(c(setdiff(design_cols, "outcome"), "sim"), names(data))
  est <- data |>
    filter(estimator %in% ests, outcome == oc) |>
    select(all_of(id), estimator, tau_hat) |>
    pivot_wider(names_from = estimator, values_from = tau_hat)
  lam <- as_diag(data) |>
    filter(outcome == oc) |>
    transmute(across(all_of(id)),
              lam_cv  = (lam_balcv1 + lam_balcv0) / 2,
              lam_end = (lam_end1 + lam_end0) / 2)
  left_join(est, lam, by = id)
}

# ---- 16. CV criterion curves (what the criterion sees) ---------------------

# simulate_grid stores the full CV loss curve for the first cv_curve_reps
# reps of every cell in <out>_cvcurves.rds. Flat curves near the minimum are
# the candidate mechanism for the slow, dispersed lambda_cv (E8).
load_cv_curves <- function(out_file) {
  d      <- tibble::as_tibble(data.table::fread(out_file))
  cells  <- distinct(d[, seq_len(match("sim", names(d)) - 1L)])  # grid order
  curves <- readRDS(sub("\\.csv\\.gz$", "_cvcurves.rds", out_file))
  stopifnot(nrow(cells) == length(curves))
  bind_rows(lapply(seq_along(curves), function(ci) {
    cc <- bind_rows(lapply(seq_along(curves[[ci]]), function(r) {
      x <- curves[[ci]][[r]]
      if (is.null(x)) return(NULL)                  # failed / unlogged rep
      bind_rows(lapply(c("treated", "control"), function(a)
        tibble(rep = r, arm = a, lambda = x$lambda[[a]], cv = x$cv.mean[[a]])))
    }))
    if (!nrow(cc)) return(NULL)
    bind_cols(cells[rep(ci, nrow(cc)), ], cc)
  }))
}

# ---- 17. Criterion insensitivity (E8: what the CV loss can't see) ----------

# Per-rep, per-arm shape of the CV loss surface from the sidecar curves.
# rel excess r(lambda) = cv(lambda)/min(cv) - 1: how much worse the criterion
# rates lambda than its own optimum. Flat region = {lambda : r < eps}; if the
# RMSE-best rung sits inside it, no selection rule on this loss can find it.
# Reference losses interpolated linearly on log-lambda; NA outside the grid.
curve_stats <- function(cur, ref_fixed = c(.05, .10), eps = .05) {
  grp <- setdiff(names(cur), c("lambda", "cv"))
  cur |>
    filter(is.finite(cv)) |>
    group_by(across(all_of(grp))) |>
    group_modify(\(d, key) {
      fmin <- min(d$cv)
      r_at <- \(l) approx(log(d$lambda), d$cv, xout = log(l), rule = 1)$y / fmin - 1
      flat <- d$lambda[d$cv / fmin - 1 < eps]
      tibble(lam_min   = d$lambda[which.min(d$cv)],
             r_end     = d$cv[which.min(d$lambda)] / fmin - 1,  # attained endpoint
             r_05      = r_at(ref_fixed[1]),
             r_10      = r_at(ref_fixed[2]),
             flat_lo   = min(flat), flat_hi = max(flat),
             flat_span = log10(max(flat) / min(flat)))          # decades wide
    }) |> ungroup()
}

# Cell summary: how flat, and which reference lambdas the criterion can't
# distinguish from its optimum. pct_* = share of stored reps with r < eps.
flatness_table <- function(cs, eps = .05, by = c(intersect(design_cols, names(cs)), "arm")) {
  cs |>
    group_by(across(all_of(by))) |>
    summarise(m_stored  = dplyr::n(),
              med_span  = median(flat_span),
              med_r_end = median(r_end),
              med_r_05  = median(r_05, na.rm = TRUE),
              med_r_10  = median(r_10, na.rm = TRUE),
              pct_end   = 100 * mean(r_end < eps),
              pct_05    = 100 * mean(r_05 < eps, na.rm = TRUE),
              pct_10    = 100 * mean(r_10 < eps, na.rm = TRUE),
              .groups = "drop")
}

# All-rep endpoint-vs-selected criterion gap from the logged two-point loss.
# Endpoint fits degenerate at small n (loss O(1e40+)), so log10 and tail
# shares only; the blowup IS the small-n endpoint-divergence evidence.
loss_gap <- function(data, by = design_cols) {
  as_diag(data) |>
    mutate(lr = pmax(log10(cvloss_end1 / cvloss_cv1),
                     log10(cvloss_end0 / cvloss_cv0))) |>
    group_by(across(all_of(intersect(by, names(data))))) |>
    summarise(m        = dplyr::n(),
              n_na     = sum(is.na(lr)),        # overflowed endpoint fits
              med_lr   = median(lr, na.rm = TRUE),
              q90_lr   = quantile(lr, .90, names = FALSE, na.rm = TRUE),
              pct_10x  = 100 * mean(lr > 1, na.rm = TRUE),
              pct_1e3x = 100 * mean(lr > 3, na.rm = TRUE),
              .groups = "drop")
}
