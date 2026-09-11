# ============================================================================
# analysis_setup.R -- shared setup for the E1-E6 analysis scripts.
# ============================================================================

# Packages and the select/filter aliases live in R/packages.R (source it first).

#' Drop the leading zero of decimals in axis labels
#'
#' @param x Character vector.
#' @return x with a leading "0." replaced by ".".
drop0 <- function(x) sub("^(-?)0\\.", "\\1.", x)   # ".05" not "0.05" on axes

# ---- 1. Load ---------------------------------------------------------------

design_cols <- c("n", "p", "s", "decay", "misspec", "covcor", "outcome", "overlap",
                 "signs", "strength", "decay_ps", "decay_out", "outcome_set",
                 "treat_prop")
outcome_cols <- c("linear", "quad1", "quad2", "exp")

#' Load one batch results file in long format
#'
#' Reads the csv.gz with data.table::fread, coerces the outcome columns to
#' numeric, replaces NA in err by "", and pivots the outcome columns to
#' (outcome, tau_hat).
#'
#' @param path Path of a results csv.gz written by simulate_grid().
#' @return A tibble with one row per replication x estimator x outcome.
load_sim <- function(path) {
  d <- tibble::as_tibble(data.table::fread(path))
  d <- mutate(d, across(any_of(outcome_cols), \(x) suppressWarnings(as.numeric(x))))
  if ("err" %in% names(d)) d$err <- coalesce(as.character(d$err), "")
  pivot_longer(d, any_of(outcome_cols), names_to = "outcome", values_to = "tau_hat")
}

# diagnostics are logged as rows of `estimator`; identify them by name
#' Identify diagnostic rows by their estimator label
#'
#' @param x Character vector of estimator labels.
#' @return Logical vector, TRUE for logged diagnostics (lam_, nnz_, smd,
#'   cvloss, trunc05, prev, emin, emax, nout).
is_diag <- function(x)
  grepl("^(lam_|nnz_|smd[01]_|cvloss_|logcvloss_|trunc05|prev|emin|emax|nout)", x) |
  x %in% c("trunc05", "prev", "emin", "emax", "nout05", "nout01")


# ---- 2. Cell summaries -----------------------------------------------------

# Screen n_na / n_extreme before reading anything else; bias/empse/rmse
# propagate NA loudly by design.
#' Cell-level performance measures with Monte Carlo standard errors
#'
#' Drops the diagnostic rows and, per group, computes bias (mean of tau_hat,
#' the true tau being 0), empirical SE, RMSE, their MCSEs, and the counts of
#' extreme (|tau_hat| > 10) and NA estimates. NA values propagate rather than
#' being dropped, by design.
#'
#' @param data Long results from load_sim().
#' @param by Grouping columns; intersected with the columns present.
#' @return A tibble with m, bias, empse, rmse, mcse_bias, mcse_empse,
#'   mcse_rmse, n_extreme and n_na per group.
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
#' Misspecified minus correct RMSE and |bias| per cell
#'
#' @param summ Output of summarise_sim() with a logical misspec column.
#' @param ests Estimators to keep.
#' @return A tibble with d_rmse and d_bias per (n, outcome, overlap, estimator).
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
#' Diagnostics as one row per replication
#'
#' Pivots the diagnostic estimator rows to columns, keyed by the design columns
#' and sim, and adds err (1 if the replication failed).
#'
#' @param data Long results from load_sim().
#' @return A wide tibble; stops if no diagnostic rows are present.
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

#' Coerce to the wide diagnostics layout
#'
#' @param data Long results or an already wide diagnostics table.
#' @return data unchanged when it has a lam_glmcv column, otherwise
#'   diag_wide(data).
as_diag <- function(data)
  if ("lam_glmcv" %in% names(data)) data else diag_wide(data)

# med / q10 / q90 columns -> long, one row per stat x quantile set
#' Pivot med / q10 / q90 summary columns to long
#'
#' @param dg Output of summarise_diag().
#' @param prefix Regular expression for the statistic name, for example
#'   "lam_[a-z0-9]+".
#' @return A tibble with stat, med, q10, q90 and the grouping columns.
pivot_q <- function(dg, prefix)
  pivot_longer(dg, matches(paste0("^", prefix, "_(med|q10|q90)$")),
               names_to = c("stat", ".value"),
               names_pattern = paste0("^(", prefix, ")_(med|q10|q90)$"))

# cell-level med / q10 / q90 of every logged quantity + failure counts
#' Cell-level median and 10/90 quantiles of every logged diagnostic
#'
#' @param data Long results or wide diagnostics.
#' @param by Grouping columns; intersected with the columns present.
#' @return A tibble with m, n_err, p_trunc05 and <stat>_med / _q10 / _q90 columns.
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
#' Attained imbalance against the lambda ceiling (balnet eq. 10)
#'
#' For whichever of the endpoint, fixed .05 / .10 and CV lambdas were logged,
#' computes the per-replication gap attained max|SMD| minus lambda and reports
#' its maximum and the share of replications exceeding tol per cell. The
#' logged SMD is self-normalised, so small positive gaps may reflect that
#' definitional mismatch rather than a solver failure.
#'
#' @param data Long results or wide diagnostics.
#' @param tol Tolerance for counting a violation.
#' @param by Grouping columns.
#' @return A tibble with m and gap*_max / gap*_pct_viol columns.
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
#' Per-replication estimates joined with that replication's lambdas
#'
#' @param data Long results from load_sim().
#' @param oc Outcome name to keep.
#' @param ests Estimators to keep, one column each.
#' @return A tibble with the design columns, sim, one column per estimator,
#'   lam_cv (mean of the two arms' selected lambdas) and lam_end (mean of the
#'   two path endpoints).
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
#' Load the stored per-replication CV loss curves of a batch
#'
#' Reads the cells (grid order) from the results csv.gz and the matching
#' <stem>_cvcurves.rds sidecar, and expands them to one row per replication,
#' arm and lambda.
#'
#' @param out_file Path of the results csv.gz.
#' @return A tibble with the grid columns, rep, arm, lambda and cv.
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
#' Shape of each CV loss curve relative to its own minimum
#'
#' Per replication and arm: the relative excess r(lambda) = cv / min(cv) - 1
#' at the path endpoint and at the reference lambdas (interpolated on log
#' lambda; NA outside the grid), and the extent of the flat region r < eps.
#'
#' @param cur Curves from load_cv_curves().
#' @param ref_fixed Two reference lambdas.
#' @param eps Flatness threshold.
#' @return A tibble with lam_min, r_end, r_05, r_10, flat_lo, flat_hi and
#'   flat_span (in decades).
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
#' Cell summary of CV curve flatness
#'
#' @param cs Output of curve_stats().
#' @param eps Flatness threshold.
#' @param by Grouping columns.
#' @return A tibble with m_stored, the median span and excesses, and the share
#'   of stored replications for which the endpoint / .05 / .10 lie within eps
#'   of the optimum.
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
#' Endpoint versus selected CV loss from the logged two-point losses
#'
#' Uses the larger (over arms) log10 ratio of the endpoint loss to the loss at
#' the selected lambda; overflowed endpoint fits give NA and are counted.
#'
#' @param data Long results or wide diagnostics.
#' @param by Grouping columns.
#' @return A tibble with m, n_na, med_lr, q90_lr, pct_10x and pct_1e3x.
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
