# what each thesis results section averages over. read by build.R (store) and thesis.R (labels)
# every runner uses base_seed 202 and stream (cell - 1) * m + r so same grid order = same draws

library(tibble)

# arms

arms_main <- c("balnet0", "abw_balnet0", "glmnetcv_hajek", "glmnetcv_ht",
               "aipw_glmnetcv", "oracle_hajek", "oracle_ht")
arms_true <- c("abw_balnet0_ortrue", "aipw_glmnetcv_ortrue")  # misspecified
arms_add <- c("ols", "oracle_aipw", "unadjusted")             # add-on files
arms_nc <- c("balnet0", "glmnetcv_hajek", "abw_balnet0", "aipw_glmnetcv")
# estimators drawn along n (fig-n-misspec). balnet and the modelling ones
arms_npath <- c("balnet0", "glmnetcv_hajek", "glmnetcv_ht", "aipw_glmnetcv")
arms_lam <- c("lam_end1", "lam_end0")                         # path ends
spec_invariant <- c("oracle_hajek", "oracle_ht")              # true e and true Y

#' balnet arm at endpoint 10^-k (estimate_tol() in R/estimators_ipw.R). 10^-4 = main balnet0
tol_arm <- function(k) ifelse(k == 4, "balnet0", paste0("balnet0_tol", k))
tol_exponents <- c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7)
arms_tol <- tol_arm(setdiff(tol_exponents, 4))               # tolerance run

# pairs. estimator A vs comparator B

# correct spec table (tbl-q1-adv)
pairs_main <- tribble(
  ~A,                   ~B,
  "balnet0",            "glmnetcv_hajek",
  "abw_balnet0",        "aipw_glmnetcv",
  "aipw_glmnetcv",      "glmnetcv_hajek",
  "balnet0",            "oracle_hajek",
  "balnet0",            "ols",
  "abw_balnet0_ortrue", "aipw_glmnetcv_ortrue"
)

# ABW-balnet vs each comparator (tbl-n-abw-m)
pairs_abw <- tribble(
  ~A,            ~B,
  "abw_balnet0", "balnet0",
  "abw_balnet0", "aipw_glmnetcv",
  "abw_balnet0", "glmnetcv_hajek",
  "abw_balnet0", "glmnetcv_ht"
)

# overlap and the other study 1 grids (tbl-ov-adv on)
pairs_s1 <- tribble(
  ~A,            ~B,
  "balnet0",     "aipw_glmnetcv",
  "balnet0",     "glmnetcv_hajek",
  "balnet0",     "glmnetcv_ht",
  "abw_balnet0", "balnet0",
  "abw_balnet0", "aipw_glmnetcv"
)

# sparsity regime table (tbl-sp-regimes) in its column order
pairs_sp <- tribble(
  ~A,            ~B,
  "balnet0",     "glmnetcv_hajek",
  "balnet0",     "aipw_glmnetcv",
  "abw_balnet0", "aipw_glmnetcv"
)

# balnet at every path endpoint vs OLS (tbl-app-tol)
pairs_tol <- tibble(A = tol_arm(tol_exponents), B = "ols")

# section helpers

# sparsity grid cov structures in display order. covcor value in results = name in tables
sparsity_corr <- c(iid = "independent", ar1 = "AR(1)",
                   equi = "equicorrelated")

#' sparsity grid. cov structure as a factor (corr) + a panel label per structure and decay
sparsity_labels <- function(d) {
  panels <- c(outer(unname(sparsity_corr), c(1, 0.5),
                    \(corr, nu) sprintf("%s, ν = %g", corr, nu)))
  d |>
    dplyr::mutate(
      corr = factor(unname(sparsity_corr[covcor]), unname(sparsity_corr)),
      structure = factor(sprintf("%s, ν = %g", corr, decay), panels))
}

#' l1 norm of the linear outcome coefs j^-nu scaled to unit var under the cell covariance
#' same cov as dgp3 .Sig in R/dgp.R
l1_norm <- function(s, nu, covcor) {
  j <- seq_len(s)
  a <- j^-nu
  sigma <- switch(covcor,
                  iid = diag(s),
                  ar1 = 0.5^abs(outer(j, j, "-")),
                  equi = 0.5 * diag(s) + 0.5)
  sum(a) / sqrt(drop(t(a) %*% sigma %*% a))
}

#' IPW-balnet linear bias as pct of the imbalance bound (path ends x l1 norm) (tbl-app-sp-bound)
#' correct score p/n >= 0.5 only. lam_end rows hold the logged value as error since tau = 0
sparsity_bound <- function(res) {
  res |>
    dplyr::filter(spec == "correct", outcome == "linear", pn >= 0.5,
                  estimator %in% c("balnet0", "lam_end1", "lam_end0")) |>
    dplyr::select(n, pn, s, decay, covcor, corr, sim, estimator, error) |>
    tidyr::pivot_wider(names_from = estimator, values_from = error) |>
    dplyr::summarise(bias = mean(balnet0),
                     mcse = sd(balnet0) / sqrt(dplyr::n()),
                     lam = mean(lam_end1) + mean(lam_end0),
                     .by = c(n, pn, s, decay, covcor, corr)) |>
    dplyr::mutate(norm = purrr::pmap_dbl(list(s, decay, covcor), l1_norm),
                  pct = -100 * bias / (lam * norm),
                  pct_mcse = 100 * mcse / (lam * norm))
}

#' share of reps where a balnet path truncates per n c and spec (fig-overlap-nc-*)
#' truncates = either arm ends above the 1e-4 target (> 1.01e-4 as ov-nc-data)
nc_truncation <- function(res) {
  res |>
    dplyr::filter(estimator %in% arms_lam, outcome == "linear") |>
    dplyr::summarise(trunc = any(error > 1.01e-4),
                     .by = c(spec, n, overlap, sim)) |>
    dplyr::summarise(trunc = mean(trunc), .by = c(spec, n, overlap))
}

#' share of reps short of exact balance per sparsity cell. same rule as nc_truncation()
sparsity_truncation <- function(res) {
  res |>
    dplyr::filter(estimator %in% arms_lam, outcome == "linear") |>
    dplyr::summarise(trunc = any(error > 1.01e-4),
                     .by = c(spec, n, pn, s, decay, covcor, corr, sim)) |>
    dplyr::summarise(reps = dplyr::n(), share = mean(trunc),
                     .by = c(spec, n, pn, s, decay, covcor, corr)) |>
    dplyr::mutate(share_mcse = sqrt(share * (1 - share) / reps))
}

#' mean balnet path ends per spec and cell w/ MCSEs. stored by both sparsity and p/n sections
#' their s = 4 nu = 1 iid cells share a ps so the ends must agree (v_across() in verify.R)
path_ends <- function(res) {
  by <- setdiff(names(res), c("sim", "estimator", "outcome", "error"))
  res |>
    dplyr::filter(estimator %in% arms_lam, outcome == "linear") |>
    tidyr::pivot_wider(names_from = estimator, values_from = error) |>
    dplyr::summarise(
      reps = dplyr::n(),
      lam1 = mean(lam_end1),
      lam1_mcse = stats::sd(lam_end1) / sqrt(reps),
      lam0 = mean(lam_end0),
      lam0_mcse = stats::sd(lam_end0) / sqrt(reps),
      lam_sum = mean(lam_end1 + lam_end0),
      lam_sum_mcse = stats::sd(lam_end1 + lam_end0) / sqrt(reps),
      .by = dplyr::all_of(by))
}

#' cost of misspec. per cell RMSE misspec / RMSE correct paired by rep then averaged like Rbar
#' ratio = exp(mean log R_c) = 1 - Rbar / 100. needs build.R and rbar.R so only runs in a build
misspec_cost <- function(res, arms, sets, miss_max = 0.01) {
  # one spec "cost" where each arm@score is its own arm so misspec vs correct is a normal pair
  est <- res |>
    dplyr::filter(estimator %in% arms) |>
    dplyr::mutate(estimator = paste0(estimator, "@", spec), spec = "cost")
  design <- est |>
    dplyr::select(-sim, -estimator, -outcome, -error) |>
    dplyr::distinct()
  purrr::map(arms, \(arm) {
    pd <- pair_draws(est, paste0(arm, "@misspec"), paste0(arm, "@correct"),
                     miss_max)
    cells <- pd$cells |>
      dplyr::left_join(pair_cell_measures(pd$draws),
                       by = c("spec", "level", "outcome")) |>
      dplyr::mutate(common = kept) |>
      dplyr::left_join(design, by = c("spec", "level"))
    average_sets(level_matrices(pd$draws), cells, list(cost = "cost"),
                 sets) |>
      dplyr::filter(theta == "pair") |>
      dplyr::mutate(arm = arm, ratio = exp(mean_log_R),
                    ratio_mcse = exp(mean_log_R) * mean_log_R_mcse,
                    .before = 1)
  }) |>
    purrr::list_rbind()
}

#' copies of a sparsity set one per cov structure. named <name>_<covcor>
within_covcor <- function(set, name) {
  keep <- if (is.null(set$keep)) \(d) rep(TRUE, nrow(d)) else set$keep
  by <- setdiff(set$by, "corr")
  if (!length(by)) by <- NULL
  copies <- lapply(names(sparsity_corr), \(v) {
    list(keep = \(d) keep(d) & d$covcor == v, by = by)
  })
  stats::setNames(copies, paste0(name, "_", names(sparsity_corr)))
}

# sparsity regimes. exact balance (p/n = 0.1) and truncated (p/n >= 0.5)
# exact_s only at n = 1000 since s = 64 is not run at n = 500 there
sparsity_regimes <- list(
  exact_s = list(keep = \(d) d$pn == 0.1 & d$n == 1000, by = "s"),
  exact = list(keep = \(d) d$pn == 0.1, by = c("decay", "corr")),
  truncated = list(keep = \(d) d$pn >= 0.5, by = c("s", "decay", "corr"))
)

# misspec cost (sparsity). rows average over n decay and outcome within regime structure and s
# exact uses n = 1000 only like exact_s so the s rows do not mix s with n
cost_arms <- c("balnet0", "glmnetcv_hajek", "abw_balnet0", "aipw_glmnetcv")
cost_regimes <- list(
  exact = list(keep = \(d) d$pn == 0.1 & d$n == 1000, by = c("s", "corr")),
  truncated = list(keep = \(d) d$pn >= 0.5, by = c("s", "corr"))
)
sparsity_cost_sets <- c(
  cost_regimes,
  list(exact_all = list(keep = cost_regimes$exact$keep),
       truncated_all = list(keep = cost_regimes$truncated$keep)),
  within_covcor(cost_regimes$exact, "exact"),
  within_covcor(cost_regimes$truncated, "truncated")
)

# sections. check_arm rows must match rep by rep across a spec's files (add-on pairing)
# spec_arms must match across specs (shared draws which the pooled MCSE assumes)

sections <- list(

  # correct spec (tbl-q1-adv fig-q1-adv tbl-q1 tbl-q1-ols)
  correct_n = list(
    title = "Correct specification, n = 500 to 400k",
    files = list(correct = list(
      "aipw_correct_n400k_dgp2.csv.gz" = arms_main,
      "aipw_add_correct_n400k_dgp2.csv.gz" = arms_add
    )),
    check_arm = "oracle_hajek",
    pairs = pairs_main, reps = 1e4, x = "n",
    sets = list(
      # drops n = 400k where glmnet jumps (appendix B)
      to_200k = list(keep = \(d) d$n <= 2e5)
    )
  ),

  # misspec (tbl-n-abw-m fig-n-abw-m)
  misspec_n = list(
    title = "Misspecified propensity score, n = 500 to 50k",
    files = list(misspec = list(
      "aipw_misspec_n_dgp2.csv.gz" = c(arms_main, arms_true),
      "aipw_add_misspec_n_dgp2.csv.gz" = arms_add
    )),
    check_arm = "oracle_hajek",
    pairs = dplyr::distinct(dplyr::bind_rows(pairs_main, pairs_abw)),
    reps = 1e4, x = "n"
  ),

  # both specs on the shared grid n = 500 to 50k (tbl-n-adv). same draws and 1e4 reps in both
  # true-OR pair left out as it was only run misspec so it cannot go in the pooled block
  n_both = list(
    title = "Both specifications, n = 500 to 50k",
    files = list(
      correct = list(
        "aipw_correct_n400k_dgp2.csv.gz" = arms_main,
        "aipw_add_correct_n400k_dgp2.csv.gz" = arms_add
      ),
      misspec = list(
        "aipw_misspec_n_dgp2.csv.gz" = arms_main,
        "aipw_add_misspec_n_dgp2.csv.gz" = arms_add
      )
    ),
    keep = \(d) d$n <= 5e4,
    check_arm = "oracle_hajek",
    spec_arms = spec_invariant,
    # tbl-n-adv also needs ABW-balnet vs balnet and glmnet (Hájek) and balnet vs AIPW-glmnet
    pairs = dplyr::distinct(dplyr::bind_rows(
      dplyr::filter(pairs_main, A != "abw_balnet0_ortrue"),
      dplyr::filter(pairs_abw, B != "glmnetcv_ht"),
      tribble(~A, ~B, "balnet0", "aipw_glmnetcv")
    )),
    reps = 1e4, x = "n"
  ),

  # along n (fig-n-misspec). ABW-balnet only for the caption. arms_true for the ps-only misspec row
  # n = 1M = 10 batches of 1000 reps read as 1 cell. fine grid file is 546 MB so build on the PC
  misspec_npath = list(
    title = "Misspecified propensity score, n = 250 to 200k and n = 1M",
    files = list(misspec = list(
      "aipw_full_misspec_nfine_dgp2.csv.gz" =
        c(arms_npath, "abw_balnet0", arms_true),
      "aipw_misspec_n_dgp2n1mill.csv.gz" =
        c(arms_npath, "abw_balnet0", arms_true)
    )),
    pairs = pairs_s1[0, ], reps = 1e4, x = "n"
  ),

  # overlap (tbl-ov-adv -lo -hi fig-ov-adv). n = 1000 from overlap7 and 10 000 from overlap8
  overlap = list(
    title = "Overlap c = 0.05 to 2, n = 1000 and 10 000",
    files = list(
      correct = list("aipw_correct_overlap7_dgp2.csv.gz" = arms_main,
                     "aipw_correct_overlap8_dgp2.csv.gz" = arms_main),
      misspec = list("aipw_misspec_overlap7_dgp2.csv.gz" = arms_main,
                     "aipw_misspec_overlap8_dgp2.csv.gz" = arms_main)
    ),
    spec_arms = spec_invariant,
    pairs = pairs_s1, reps = 1e4, x = "overlap", by = "n",
    sets = list(
      poor = list(keep = \(d) d$overlap <= 0.25),   # tbl-ov-adv-lo
      good = list(keep = \(d) d$overlap >= 0.5),    # tbl-ov-adv-hi
      by_c = list(by = "overlap")                   # fig-ov-adv
    )
  ),

  # overlap n x c grid (fig-overlap-nc-200 -840 -5000)
  nc = list(
    title = "n x c grid, n = 100 to 5000, c = 0.01 to 1, 1000 replications",
    files = list(
      correct = list("aipw_correct_n_by_overlap_dgp2.csv.gz" =
                       c(arms_nc, spec_invariant, arms_lam)),
      misspec = list("aipw_misspec_n_by_overlap_dgp2.csv.gz" =
                       c(arms_nc, spec_invariant, arms_lam))
    ),
    spec_arms = spec_invariant,
    pairs = pairs_s1[0, ], reps = 1e3, x = "overlap", by = "n",
    extra = list(trunc = nc_truncation)
  ),

  # high dim p/n grid (tbl-pn-adv fig-p-adv)
  p_axis = list(
    title = "p/n = 0.05 to 4, n = 500 and 1000, s = 4",
    # path ends too. sparsity grid s = 4 nu = 1 iid cells share this grid's ps (path_ends())
    files = list(
      correct = list("aipw_correct_p_dgp2.csv.gz" = c(arms_main, arms_lam)),
      misspec = list("aipw_misspec_p_dgp2.csv.gz" =
                       c(arms_main, arms_true, arms_lam))
    ),
    spec_arms = spec_invariant,
    derive = \(d) dplyr::mutate(d, pn = p / n),
    pairs = pairs_s1, reps = 1e4, x = "pn", by = "n",
    sets = list(
      exact = list(keep = \(d) d$pn < 0.5),         # exact balance holds
      truncated = list(keep = \(d) d$pn >= 0.5)     # balnet path truncates
    ),
    extra = list(path_ends = path_ends)
  ),

  # noise covariates and denser profile (tbl-sp-regimes fig-sp-adv tbl-sp-adv tbl-app-sp-bound)
  # also tbl-app-avg-sp-s -exact -trunc
  sparsity = list(
    title = "Sparsity and correlation, correct score and Tan's transformation",
    # equi has its own 34 cell grid w/ dgp3 (dgp2 + equi cov). misspec = Tan's transform of obs X
    # (dgp3 tan_obs) on the correct files' draws. runs/aipw/sparsity_axis.R
    files = list(
      correct = list(
        "aipw_correct_sparsity_hd_dgp2.csv.gz" = c(arms_main, arms_lam),
        "aipw_correct_sparsity_hd_equi_dgp3.csv.gz" = c(arms_main, arms_lam)
      ),
      misspec = list(
        "aipw_misspec_tanobs_sparsity_hd_dgp3.csv.gz" =
          c(arms_main, arms_true, arms_lam),
        "aipw_misspec_tanobs_sparsity_hd_equi_dgp3.csv.gz" =
          c(arms_main, arms_true, arms_lam)
      )
    ),
    spec_arms = spec_invariant,
    derive = sparsity_labels,
    pairs = pairs_s1, reps = 1e3, x = "s", by = "structure",
    sets = c(
      sparsity_regimes,
      # each regime over all its cells
      list(exact_all = list(keep = sparsity_regimes$exact$keep),
           truncated_all = list(keep = sparsity_regimes$truncated$keep)),
      # within each cov structure since trends differ (balnet lead narrows w/ s under equi)
      within_covcor(sparsity_regimes$exact_s, "exact_s"),
      # all cells one avg per cov structure. for the text only and no table shows it
      list(by_corr = list(by = "corr")),
      within_covcor(sparsity_regimes$exact, "exact"),
      within_covcor(sparsity_regimes$truncated, "truncated")
    ),
    extra = list(
      bound = sparsity_bound,
      trunc = sparsity_truncation,
      path_ends = path_ends,
      cost = \(res) misspec_cost(res, cost_arms, sparsity_cost_sets)
    ),
    cost = list(arms = cost_arms, sets = sparsity_cost_sets)  # verify.R recomputes cost from these
  ),

  # outcome signal (tbl-sig-adv -weak -strong fig-sig-adv)
  signal = list(
    title = "Outcome signal strength 0.02 to 1, n = 1000",
    files = list(
      correct = list("aipw_correct_strength_dgp2.csv.gz" = arms_main),
      misspec = list("aipw_misspec_strength_dgp2.csv.gz" =
                       c(arms_main, arms_true))
    ),
    spec_arms = spec_invariant,
    pairs = pairs_s1, reps = 1e4, x = "strength",
    # split at kappa 0.1 vs 0.2 where balnet rank vs the modelling ests flips (correct spec)
    sets = list(
      weak = list(keep = \(d) d$strength <= 0.1),    # tbl-sig-adv-weak
      strong = list(keep = \(d) d$strength >= 0.2),  # tbl-sig-adv-strong
      by_kappa = list(by = "strength")              # tbl-app-avg-bykappa-*
    )
  ),

  # appendix H treated share (tbl-ts-adv)
  treat_share = list(
    title = "Treated share, n = 1000 and 10 000",
    files = list(
      correct = list("aipw_correct_treat_prop_dgp2.csv.gz" = arms_main),
      misspec = list("aipw_misspec_treat_prop_dgp2.csv.gz" =
                       c(arms_main, arms_true))
    ),
    spec_arms = spec_invariant,
    pairs = pairs_s1, reps = 1e4, x = "treat_prop", by = "n"
  ),

  # appendix 2 IPW-balnet at endpoints 1e-2 to 1e-7 vs OLS (tbl-app-tol). 1e-4 is correct_n balnet0
  # tol run (runs/aipw/aipw_tol_dgp2.R) stops at 50k. the 3 files share draws by n and sim to there
  tol_n = list(
    title = "Path endpoints 1e-2 to 1e-7, correct specification, n = 500 to 50k",
    files = list(correct = list(
      "aipw_correct_n400k_dgp2.csv.gz" = "balnet0",
      "aipw_add_correct_n400k_dgp2.csv.gz" = "ols",
      "aipw_tol_correct_n_dgp2.csv.gz" = arms_tol
    )),
    keep = \(d) d$n <= 5e4,
    check_arm = "oracle_hajek",
    pairs = pairs_tol, reps = 1e4, x = "n"
  )
)

# appendix tables. same fields as sections plus dir full diag approx simple display
# full$B = MAE bootstrap resamples (0 = no MAE MCSE). seed restarts per spec like run_perf() did

# diag rows logged by estimate_all() in R/estimators_ipw.R. weight sets balnet glmnet true
arms_diag <- c(
  as.vector(outer(as.vector(outer(c("ess", "wmax", "smd"),
                                  c("", "_glm", "_true"), paste0)),
                  c("1", "0"), paste0)),
  "nnz_ps", "nnz_bal1", "nnz_bal0", "nnz_or1", "nnz_or0",
  "lam_end1", "lam_end0", "time_bal", "time_ps", "time_or"
)

spec_names <- c(correct = "Correct specification",
                misspec = "Misspecified propensity score")

appendix_sections <- list(

  # appendix misspec (tbl-app-n-perf-m-* tbl-app-n-diag-m-*
  # tbl-apx-bias-n-m tbl-apx-bias-n-m-avg)
  app_misspec_n = list(
    title = "Appendix: misspecified propensity score, n = 500 to 50k",
    files = sections$misspec_n$files,
    reps = 1e4,
    full = list(drop = c("oracle_ht", "ols", "unadjusted"), B = 1000,
                seed = 2026),
    diag = TRUE,
    approx = list(A = "abw_balnet0", B = "ols"),
    display = list(x = "n", x_fmt = "thousands", x_title = "*n* (thousands)",
                   x_name = "n", spec = spec_names["misspec"])
  ),

  # appendix p/n grid (app-pn)
  app_p_axis = list(
    title = "Appendix: p/n = 0.05 to 4, n = 500 and 1000",
    files = sections$p_axis$files,
    derive = sections$p_axis$derive,
    reps = 1e4,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "pn", by = "n", x_fmt = "plain", x_title = "*p*/*n*",
                   x_name = "p/n", by_fmt = "n = %s", spec = spec_names)
  ),

  # appendix sparsity grid (app-sp)
  app_sparsity = list(
    title = "Appendix: sparsity and correlation, correct specification",
    files = list(correct = list(
      "aipw_correct_sparsity_hd_dgp2.csv.gz" = arms_main
    )),
    reps = 1e3,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "s", x_fmt = "plain", x_title = "s", x_name = "s",
                   design_labs = c(pn = "p/n", decay = "Decay",
                                   covcor = "Covariates"),
                   spec = spec_names["correct"])
  ),

  # appendix sparsity equi cells correct score. own entry since adding them to app_sparsity
  # would move every MAE MCSE there (bootstrap goes in grid_key order and these keys interleave)
  app_sparsity_equi = list(
    title = paste("Appendix: sparsity, equicorrelated covariates, correct",
                  "specification"),
    files = list(correct = list(
      "aipw_correct_sparsity_hd_equi_dgp3.csv.gz" = arms_main
    )),
    reps = 1e3,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "s", x_fmt = "plain", x_title = "s", x_name = "s",
                   design_labs = c(pn = "p/n", decay = "Decay",
                                   covcor = "Covariates"),
                   spec = spec_names["correct"])
  ),

  # appendix sparsity under Tan's transformation. every cell iid AR(1) and equi
  app_sparsity_tan = list(
    title = "Appendix: sparsity, Tan's transformation",
    files = list(misspec = list(
      "aipw_misspec_tanobs_sparsity_hd_dgp3.csv.gz" = c(arms_main, arms_true),
      "aipw_misspec_tanobs_sparsity_hd_equi_dgp3.csv.gz" =
        c(arms_main, arms_true)
    )),
    reps = 1e3,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "s", x_fmt = "plain", x_title = "s", x_name = "s",
                   design_labs = c(pn = "p/n", decay = "Decay",
                                   covcor = "Covariates"),
                   spec = spec_names["misspec"])
  ),

  # appendix n x c grid (app-nc-perf). no MAE so no bootstrap
  app_nc = list(
    title = "Appendix: n x c grid, both specifications",
    files = list(
      correct = list("aipw_correct_n_by_overlap_dgp2.csv.gz" = arms_main),
      misspec = list("aipw_misspec_n_by_overlap_dgp2.csv.gz" =
                       c(arms_main, arms_true))
    ),
    reps = 1e3,
    full = list(drop = "oracle_ht", B = 0, seed = 2026),
    display = list(x = "overlap", x_fmt = "plain", x_title = "Overlap *c*",
                   x_name = "c", spec = spec_names)
  ),

  # appendix 1 IPW-balnet to 800k (tbl-balnet-n800k). old run runs/ipw/correctspec_dgp2_n800k.R
  # seeds 202 (n <= 200k) 203 (400k) 204 (800k). linear and quad1 only. not paired w/ correct_n
  app_ipw_800k = list(
    title = "Appendix: IPW-balnet to n = 800 000 (IPW-set run)",
    dir = "ipw",  # folder in results/ (default aipw)
    files = list(correct = list(
      "correctspecDGP2.csv.gz" = c("balnet0", "balnetcv"),
      "correctspecDGP400k.csv.gz" = c("balnet0", "balnetcv"),
      "correctspecDGP800k.csv.gz" = c("balnet0", "balnetcv")
    )),
    reps = 1e3,
    simple = TRUE,  # bias EmpSE RMSE only and every rep counted
    full = list(B = 1000, seed = 2026),             # both arms none dropped
    display = list(x = "n", x_fmt = "thousands", x_title = "*n* (thousands)",
                   x_name = "n", spec = spec_names["correct"])
  ),

  # appendix A correct spec (tbl-app-q1-full-* tbl-app-q1-diag-* tbl-q1-approx tbl-q1-approx-avg)
  # OLS and the other add-on arms from the add-on file on correct_n draws (paired there)
  app_correct_n = list(
    title = "Appendix: correct specification, n = 500 to 400k",
    files = sections$correct_n$files,
    reps = 1e4,
    full = list(B = 1000, seed = 2026),             # every arm none dropped
    diag = TRUE,
    approx = list(A = "abw_balnet0", B = "ols"),
    display = list(x = "n", x_fmt = "thousands", x_title = "*n* (thousands)",
                   x_name = "n", spec = spec_names["correct"])
  ),

  # appendix approx bias ABW-balnet vs OLS n = 1000 (tbl-apx-bias-ov-*). OLS from add-on runs
  # (runs/aipw/aipw_addon_dgp2.R). different files so check_arm must match or nothing is stored
  app_overlap = list(
    title = "Appendix: overlap c = 0.05 to 2, n = 1000",
    files = list(
      correct = list("aipw_correct_overlap7_dgp2.csv.gz" = "abw_balnet0",
                     "aipw_diag_correct_overlap7_dgp2.csv.gz" = "ols"),
      misspec = list("aipw_misspec_overlap7_dgp2.csv.gz" = "abw_balnet0",
                     "aipw_full_misspec_overlap7_dgp2.csv.gz" = "ols")
    ),
    check_arm = "oracle_hajek",
    reps = 1e4,
    approx = list(A = "abw_balnet0", B = "ols"),
    display = list(x = "overlap", x_fmt = "plain", x_title = "Overlap *c*",
                   x_name = "c", spec = spec_names)
  ),

  # appendix C overlap full perf and diag minus Oracle (HT) (tbl-app-ov-full-* tbl-app-ov-diag-*)
  # overlap section files log the diagnostics too
  app_overlap_full = list(
    title = "Appendix: overlap c = 0.05 to 2, n = 1000 and 10 000",
    files = sections$overlap$files,
    reps = 1e4,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "overlap", by = "n", x_fmt = "plain",
                   x_title = "Overlap *c*", x_name = "c", by_fmt = "n = %s",
                   spec = spec_names)
  ),

  # appendix G signal full perf and diag minus Oracle (HT) (tbl-app-sig-full-* tbl-app-sig-diag-*)
  # signal section files log the diagnostics too
  app_signal = list(
    title = "Appendix: outcome signal strength 0.02 to 1, n = 1000",
    files = sections$signal$files,
    reps = 1e4,
    full = list(drop = "oracle_ht", B = 1000, seed = 2026),
    diag = TRUE,
    display = list(x = "strength", x_fmt = "plain",
                   x_title = "Signal strength *κ*", x_name = "κ",
                   spec = spec_names)
  ),

  # appendix approx bias on the p/n grid (tbl-apx-bias-pn-*). OLS from diag twins on same draws
  # OLS is NA once p + 1 > the smaller arm (p/n >= 0.5) so those cells are left out
  app_p_approx = list(
    title = "Appendix: approximation bias, p/n = 0.05 to 4, n = 500 and 1000",
    files = list(
      correct = list("aipw_correct_p_dgp2.csv.gz" = "abw_balnet0",
                     "aipw_diag_correct_p_dgp2.csv.gz" = "ols"),
      misspec = list("aipw_misspec_p_dgp2.csv.gz" = "abw_balnet0",
                     "aipw_diag_misspec_p_dgp2.csv.gz" = "ols")
    ),
    derive = sections$p_axis$derive,
    check_arm = "oracle_hajek",
    reps = 1e4,
    approx = list(A = "abw_balnet0", B = "ols"),
    display = list(x = "pn", by = "n", x_fmt = "plain", x_title = "*p*/*n*",
                   x_name = "p/n", by_fmt = "n = %s", spec = spec_names)
  ),

  # appendix approx bias fine n grid misspec 250 to 200k 29 levels (tbl-apx-bias-nfine-m*)
  # both arms in one 546 MB file so build this entry on the PC
  app_misspec_nfine = list(
    title = "Appendix: misspecified propensity score, fine n grid to 200k",
    files = list(misspec = list(
      "aipw_full_misspec_nfine_dgp2.csv.gz" = c("abw_balnet0", "ols")
    )),
    reps = 1e4,
    approx = list(A = "abw_balnet0", B = "ols"),
    display = list(x = "n", x_fmt = "thousands", x_title = "*n* (thousands)",
                   x_name = "n", spec = spec_names["misspec"])
  )
)
