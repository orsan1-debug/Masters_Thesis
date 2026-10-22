# study 2 choosing lambda on the balnet path. results/cv/our dgp runs to one stored summary each
# tau = 0 in every grid so each est is its own error. measures only from rbar.R

library(dplyr)
library(tidyr)
library(purrr)

# runs

# est rows <est>_<arm>. est ipw or abw. arm = a selector or fixed (path end) or path_01 to 30
# lam1_ lam0_ = lambda each arm chose (treated and control fit). smd* ess* not summarised here
cv_selectors <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
cv_outcomes <- c("linear", "quad1", "exp")
cv_path_k <- 30   # points on each cell's lambda grid. path_30 is the end
# interior optimum = path min below this share of the path end rmse
# stored `interior` col keeps the 8 Oct rule (gain_floor > 0.05)
cv_interior_at <- 0.95

# per run files (in results/cv/our dgp) ps spec and reps. merge = fold runs that only add cv rows
# on the same draws. pc_only = files over 400 MB so built on the pc
cv_sections <- list(
  cv_overlap = list(
    title = "CV, overlap axis",
    files = c(main = "r2500/cv_overlap_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500)),
  cv_noise = list(
    title = "CV, outcome strength axis",
    files = c(main = "r2500/cv_noise_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500)),
  cv_noise_x_overlap = list(
    title = "CV, outcome strength by overlap",
    files = c(main = "r2500/cv_noise_x_overlap_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500)),
  cv_misspec = list(
    title = "CV, misspecified propensity score",
    files = c(main = "r2500/cv_misspec_dgp2.csv.gz"),
    spec = "misspec", reps = c(main = 2500)),
  cv_folds = list(
    title = "CV, number of folds",
    files = c(k10 = "r2500/cv_folds_k10_dgp2.csv.gz",
              k3 = "r2500/cv_folds_k3_dgp2.csv.gz",
              k5 = "r2500/cv_folds_k5_dgp2.csv.gz"),
    merge = c(k3 = "_k3", k5 = "_k5"),
    spec = "correct", reps = c(k10 = 2500, k3 = 2500, k5 = 2500)),
  cv_sparsity = list(
    title = "CV, dense confounding by covariance",
    files = c(main = "r2500/cv_sparsity_dgp3.csv.gz"),
    spec = "correct", reps = c(main = 2500)),
  cv_sparsity_misspec = list(
    title = "CV, dense confounding by covariance, misspecified score",
    files = c(main = "r2500/cv_sparsity_misspec_dgp3.csv.gz"),
    spec = "misspec", reps = c(main = 2500)),
  cv_nbig = list(
    title = "CV, large n by overlap",
    files = c(n5000 = "r2500/cv_n5000_dgp2.csv.gz",
              n10000 = "r1000/cv_n10000_dgp2.csv.gz"),
    spec = "correct", reps = c(n5000 = 2500, n10000 = 1000)),
  cv_alpha_x_overlap = list(
    title = "CV, penalty alpha by overlap",
    files = c(main = "r2500/cv_alpha_x_overlap_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500), pc_only = TRUE),
  cv_n = list(
    title = "CV, n by overlap by alpha",
    files = c(main = "r2500/cv_n_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500), pc_only = TRUE),
  cv_p = list(
    title = "CV, p by overlap by alpha",
    files = c(main = "r2500/cv_p_dgp2.csv.gz"),
    spec = "correct", reps = c(main = 2500), pc_only = TRUE)
)

#' selector arms of a run. the five plus the cv.balnet ones of each merged file w/ suffix
cv_section_selectors <- function(section) {
  cv_only <- grep("^cv\\.", cv_selectors, value = TRUE)
  c(cv_selectors,
    unlist(map(unname(section$merge), \(suffix) paste0(cv_only, suffix))))
}

# read

#' one run file as stored. err read as text (empty unless the rep failed)
cv_read <- function(file) {
  data.table::fread(here::here("results", "cv", "our dgp", file),
                    colClasses = list(character = "err"))
}

#' lambda grid a run file was scored on
cv_read_lamgrid <- function(file) {
  path <- here::here("results", "cv", "our dgp",
                     sub("\\.csv\\.gz$", "_lamgrid.csv", file))
  tibble::as_tibble(data.table::fread(path))
}

#' error cols. one per outcome or one per outcome and strength (linear_s1 etc)
cv_error_cols <- function(d) {
  if (all(cv_outcomes %in% names(d))) return(cv_outcomes)
  grep("^(linear|quad1|exp)_s", names(d), value = TRUE)
}

#' design cols. everything that is not an id or est row or error
cv_design_cols <- function(d) {
  setdiff(names(d), c("cell", "sim", "estimator", "err", cv_error_cols(d)))
}

#' ipw and abw rows of one cell in long form. one row per draw est row outcome and strength
#' level pastes the design cols and strength so rbar.R treats each setting as one cell
cv_long <- function(d, spec) {
  design <- cv_design_cols(d)
  cols <- cv_error_cols(d)
  d <- d[grepl("^(ipw|abw)_", d$estimator), c("cell", design, "sim",
                                               "estimator", cols),
         with = FALSE]
  long <- data.table::melt(d, measure.vars = cols, variable.name = "column",
                           value.name = "error", variable.factor = FALSE)
  long <- tibble::as_tibble(long)
  if (identical(cols, cv_outcomes)) {
    long <- rename(long, outcome = column)
  } else {
    long <- long |>
      separate_wider_delim(column, "_s", names = c("outcome", "strength")) |>
      mutate(strength = as.numeric(strength))
  }
  long |>
    mutate(
      spec = spec,
      level = do.call(paste, c(map(pick(all_of(c(design, "strength"))),
                                   key_text), sep = "_")),
      est = sub("_.*", "", estimator),
      arm = sub("^(ipw|abw)_", "", estimator),
      k = suppressWarnings(as.integer(sub("^path_", "", arm))),
      arm = if_else(is.na(k), arm, "path"),
      # failed rep leaves the estimate empty
      error = suppressWarnings(as.numeric(error)),
      .before = 1) |>
    relocate(spec, level, cell, all_of(c(design, "strength")), sim, est, arm,
             k, estimator, outcome, error)
}

# checks
# each check is one row (check_row() in build.R). nothing stored for a run that fails one

#' checks on one run file before anything is computed from it
cv_check_file <- function(d, lamgrid, reps, name) {
  per_cell <- d[, .(n = .N, sims = data.table::uniqueN(sim)),
                by = .(cell, estimator)]
  cols <- cv_error_cols(d)
  lam <- d[grepl("^lam[01]_", d$estimator), cols, with = FALSE]
  # failed rep leaves every col empty which is fine. empty in only some cols is not
  row_spread <- \(v) {
    ok <- is.finite(v)
    if (!any(ok)) 0 else if (!all(ok)) Inf else diff(range(v))
  }
  lam_spread <- if (nrow(lam)) max(apply(as.matrix(lam), 1, row_spread)) else 0
  # path_30 is the grid end so it must equal fixed draw for draw
  # the cells table measures the path end by `fixed`
  ends <- map_dbl(c("ipw", "abw"), \(e) {
    a <- d[d$estimator == paste0(e, "_fixed"), c("cell", "sim", cols),
           with = FALSE]
    b <- d[d$estimator == sprintf("%s_path_%02d", e, cv_path_k),
           c("cell", "sim", cols), with = FALSE]
    data.table::setkey(a, cell, sim)
    data.table::setkey(b, cell, sim)
    if (nrow(a) != nrow(b) || !all(a$sim == b$sim)) return(Inf)
    max(abs(as.matrix(a[, cols, with = FALSE]) -
              as.matrix(b[, cols, with = FALSE])), na.rm = TRUE)
  })
  failed <- sum(nzchar(d$err))
  grid_k <- lamgrid |> count(cell)
  bind_rows(
    check_row(sprintf("%s: %d replications of every row in every cell",
                      name, reps),
              all(per_cell$n == reps & per_cell$sims == reps),
              sprintf("%d cells, %d estimator rows, n per cell %s",
                      data.table::uniqueN(d$cell),
                      data.table::uniqueN(d$estimator),
                      paste(range(per_cell$n), collapse = " to "))),
    check_row(sprintf("%s: replications that failed (err not empty)", name),
              TRUE, sprintf("%d rows", failed)),
    check_row(sprintf("%s: path_%02d equals the fixed arm in every draw",
                      name, cv_path_k),
              all(ends < 1e-10),
              sprintf("largest difference %.1e (IPW-balnet), %.1e (ABW)",
                      ends[1], ends[2])),
    check_row(sprintf("%s: chosen lambdas the same in every error column",
                      name),
              lam_spread == 0, sprintf("largest spread %.1e", lam_spread)),
    check_row(sprintf("%s: a lambda grid of %d points for every cell", name,
                      cv_path_k),
              setequal(grid_k$cell, unique(d$cell)) &&
                all(grid_k$n == cv_path_k),
              sprintf("%d cells in the grid file", nrow(grid_k))))
}

#' cv.balnet rows of a merged file (estimates and lambdas) renamed w/ its suffix
cv_merge_rows <- function(d, suffix) {
  d <- d[grepl("^(ipw|abw|lam[01])_cv\\.", d$estimator), ]
  d$estimator <- paste0(d$estimator, suffix)
  d
}

#' check a merged file has the same draws. path and path end must equal the first file's
cv_check_merge <- function(d, dm, name) {
  shared <- grep("^(ipw|abw)_(fixed|path_\\d+)$", unique(d$estimator),
                 value = TRUE)
  cols <- cv_error_cols(d)
  key <- c("cell", "sim", "estimator")
  a <- d[d$estimator %in% shared, c(key, cols), with = FALSE]
  b <- dm[dm$estimator %in% shared, c(key, cols), with = FALSE]
  data.table::setkeyv(a, key)
  data.table::setkeyv(b, key)
  same_rows <- nrow(a) == nrow(b) && isTRUE(all.equal(a[, ..key], b[, ..key],
                                                      check.attributes = FALSE))
  gap <- if (same_rows) {
    max(abs(as.matrix(a[, ..cols]) - as.matrix(b[, ..cols])), na.rm = TRUE)
  } else {
    Inf
  }
  check_row(sprintf(
    "%s: path and path end equal the first file's in every draw", name),
            same_rows && gap < 1e-10,
            sprintf("%d rows, largest difference %.1e", nrow(b), gap))
}

# measures

#' chosen arm vs path end on the same draws. R_c and mcse from rbar.R as log2 ratios
cv_against_end <- function(long, choice, end) {
  rows <- bind_rows(
    long |>
      semi_join(choice, by = c("spec", "level", "outcome", "estimator")) |>
      mutate(estimator = "chosen"),
    long |>
      filter(estimator == end) |>
      mutate(estimator = "end"))
  pd <- pair_draws(rows, "chosen", "end")
  pair_cell_measures(pd$draws) |>
    transmute(spec, level, outcome,
              log2_R = log2(R_c), log2_R_mcse = log_R_c_mcse / log(2))
}

#' cells table of one est. rules of the 8 Oct summaries
#' path min is over the whole grid incl the end so gain_floor is never negative
cv_cells_est <- function(perf, long, e, selectors = cv_selectors) {
  p <- filter(perf, est == e)
  l <- filter(long, est == e)
  keys <- c("spec", "level", "outcome")
  end <- p |>
    filter(arm == "fixed") |>
    select(all_of(keys), fixed = rmse, se_f = mcse_rmse)
  path_min <- p |>
    filter(arm == "path") |>
    arrange(k) |>
    slice_min(rmse, n = 1, with_ties = FALSE, by = all_of(keys)) |>
    select(all_of(keys), estimator, k_opt = k, path_min = rmse,
           se_m = mcse_rmse)
  best <- p |>
    filter(arm %in% selectors) |>
    mutate(arm = factor(arm, selectors)) |>
    arrange(arm) |>
    slice_min(rmse, n = 1, with_ties = FALSE, by = all_of(keys)) |>
    transmute(across(all_of(keys)), estimator, best = as.character(arm),
              best_rmse = rmse, se_b = mcse_rmse)
  end_row <- paste0(e, "_fixed")
  ratio_best <- cv_against_end(l, best, end_row) |>
    rename(log2_R_best = log2_R, log2_R_best_mcse = log2_R_mcse)
  ratio_min <- cv_against_end(l, path_min, end_row) |>
    rename(log2_R_min = log2_R, log2_R_min_mcse = log2_R_mcse)
  end |>
    left_join(select(path_min, -estimator), by = keys) |>
    left_join(select(best, -estimator), by = keys) |>
    left_join(ratio_best, by = keys) |>
    left_join(ratio_min, by = keys) |>
    mutate(est = e, .after = outcome) |>
    mutate(gain_floor = fixed / path_min - 1,
           interior = gain_floor > 0.05,
           beats_floor = best_rmse < fixed - 2 * sqrt(se_b^2 + se_f^2),
           sel_vs_fixed = best_rmse / fixed - 1,
           sel_vs_min = best_rmse / path_min - 1,
           reaches = beats_floor & sel_vs_min < 0.05)
}

#' pairs within one est. each selector and grid point vs the path end and cv vs boot
#' cv.smd vs boot.smd and cv.inf vs boot.inf. other fold counts vs the same boot arm
cv_pairs <- function(selectors) {
  cv_crit <- grep("^cv\\.(smd|inf)", selectors, value = TRUE)
  bind_rows(
    tibble(kind = "selector", A = selectors, B = "fixed"),
    tibble(kind = "path", A = sprintf("path_%02d", seq_len(cv_path_k - 1)),
           B = "fixed"),
    tibble(kind = "cv_boot", A = cv_crit,
           B = sub("^cv\\.(smd|inf).*$", "boot.\\1", cv_crit)))
}

#' sets of one design cell's strengths and outcomes to average over for one est
#' all or interior or outside (cv_interior_at) or one strength (kappa)
cv_sets <- function(perf, e, strengths) {
  cells <- perf |>
    filter(est == e) |>
    summarise(interior = min(rmse[arm == "path"]) <
                cv_interior_at * rmse[arm == "fixed"],
              .by = c(spec, level, outcome)) |>
    inner_join(strengths, by = c("spec", "level", "outcome"))
  bind_rows(
    mutate(cells, set = "all", value = NA_real_),
    mutate(filter(cells, interior), set = "interior", value = NA_real_),
    mutate(filter(cells, !interior), set = "outside", value = NA_real_),
    mutate(cells, set = "kappa", value = strength)) |>
    select(spec, level, outcome, set, value)
}

#' one arm's errors in pair_draws() form so its log rmse sums over cells like log R_c
#' for ratios between design cells that don't share draws (penalty forms)
cv_arm_draws <- function(rows, miss_max) {
  rows |>
    select(spec, level, sim, outcome, error) |>
    filter(mean(!is.finite(error)) <= miss_max,
           .by = c(spec, level, outcome)) |>
    filter(is.finite(error)) |>
    mutate(influence_log_R = influence_log_mse(error) / 2,
           influence_s_A = 0, influence_s_B = 0,
           .by = c(spec, level, outcome))
}

#' sums of log R_c (pair) or log rmse (arm) over the sets of one design cell w/ mcse
#' strengths and outcomes share every draw so for the mcse they count as cells of one level
cv_set_sums <- function(draws, values, sets) {
  one_level <- \(x) mutate(x, outcome = paste(level, outcome, sep = "/"),
                           level = "design")
  by_level <- level_matrices(one_level(draws))
  sets |>
    inner_join(values, by = c("spec", "level", "outcome")) |>
    summarise(
      n_cells = n(),
      sum_log = sum(log),
      sum_log_mcse = n() * mean_mcse_over_cells(
        by_level, one_level(tibble(spec, level, outcome))
      )[["influence_log_R"]],
      .by = c(set, value))
}

#' pairs and set sums of one est in one design cell (plus the arms fixed and path_min)
#' design cells have own seed streams so sums over cells or runs add w/ their variances
cv_compare_est <- function(perf, long, e, selectors, strengths, miss_max) {
  rows <- split(filter(long, est == e), ~estimator)
  sets <- cv_sets(perf, e, strengths)
  keys <- c("spec", "level", "outcome")
  compared <- pmap(cv_pairs(selectors), \(kind, A, B) {
    a <- paste0(e, "_", A)
    b <- paste0(e, "_", B)
    pd <- pair_draws(bind_rows(rows[[a]], rows[[b]]), a, b, miss_max)
    cells <- pair_cell_measures(pd$draws) |>
      select(all_of(keys), reps, rmse_A, rmse_B, R_c, log_R_c_mcse, won,
             lost) |>
      mutate(kind = kind, A = A, B = B, .before = 1)
    sums <- if (kind != "path") {
      cv_set_sums(pd$draws, transmute(cells, across(all_of(keys)),
                                      log = log(R_c)), sets) |>
        mutate(kind = kind, A = A, B = B, .before = 1)
    }
    list(pairs = cells, sums = sums)
  })
  opt <- perf |>
    filter(est == e, arm == "path") |>
    arrange(k) |>
    slice_min(rmse, n = 1, with_ties = FALSE, by = all_of(keys)) |>
    select(all_of(keys), estimator)
  arms <- list(fixed = rows[[paste0(e, "_fixed")]],
               path_min = semi_join(filter(long, est == e, arm == "path"),
                                    opt, by = c(keys, "estimator")))
  arm_sums <- imap(arms, \(r, A) {
    draws <- cv_arm_draws(r, miss_max)
    values <- summarise(draws, log = log(mean(error^2)) / 2,
                        .by = all_of(keys))
    cv_set_sums(draws, values, sets) |>
      mutate(kind = "arm", A = A, B = NA_character_, .before = 1)
  })
  list(pairs = list_rbind(map(compared, "pairs")) |> mutate(est = e),
       sums = bind_rows(list_rbind(map(compared, "sums")),
                        list_rbind(arm_sums)) |>
         mutate(est = e))
}

#' perf cells picks pairs and sums of one cell of a run file
cv_one_cell <- function(d, spec, miss_max, selectors = cv_selectors) {
  long <- cv_long(d, spec)
  arms <- distinct(long, estimator, est, arm, k)
  perf <- cell_performance(long, miss_max) |>
    left_join(arms, by = "estimator")
  cells <- map(c("ipw", "abw"), \(e)
    cv_cells_est(perf, long, e, selectors)) |>
    list_rbind()
  design <- cv_design_cols(d)
  first_col <- cv_error_cols(d)[1]
  picks <- tibble::as_tibble(d[grepl("^lam[01]_", d$estimator),
                               c("cell", design, "estimator", first_col),
                               with = FALSE]) |>
    rename(lambda = all_of(first_col)) |>
    summarise(median_lambda = median(lambda, na.rm = TRUE),
              .by = c(cell, all_of(design), estimator)) |>
    mutate(side = if_else(startsWith(estimator, "lam1"), "treated",
                          "control"),
           arm = sub("^lam[01]_", "", estimator), .before = estimator) |>
    select(-estimator)
  level_design <- distinct(long, spec, level, cell,
                           across(all_of(c(design, "strength"))))
  strengths <- distinct(long, spec, level, outcome, strength)
  compared <- map(c("ipw", "abw"), \(e)
    cv_compare_est(perf, long, e, selectors, strengths, miss_max))
  # a sum spans the strengths so it only gets the design cols
  design_cell <- distinct(long, cell, across(all_of(design)))
  list(perf = left_join(perf, level_design, by = c("spec", "level")),
       cells = left_join(cells, level_design, by = c("spec", "level")),
       picks = picks,
       pairs = list_rbind(map(compared, "pairs")) |>
         left_join(level_design, by = c("spec", "level")),
       sums = list_rbind(map(compared, "sums")) |>
         cross_join(design_cell))
}

# build

#' build one run's summary. read and check files then do one cell at a time
#' (a whole r2500 file of the n or p grid is ~70 million rows in long form)
build_cv_section <- function(section, id, miss_max = 0.01) {
  merged <- names(section$merge)
  selectors <- cv_section_selectors(section)
  parts <- setdiff(names(section$files), merged)
  out <- map(set_names(parts), \(name) {
    file <- section$files[[name]]
    d <- cv_read(file)
    lamgrid <- cv_read_lamgrid(file)
    checks <- cv_check_file(d, lamgrid, section$reps[[name]], name)
    for (m in merged) {
      dm <- cv_read(section$files[[m]])
      checks <- bind_rows(
        checks,
        cv_check_file(dm, cv_read_lamgrid(section$files[[m]]),
                      section$reps[[m]], m),
        cv_check_merge(d, dm, m))
      d <- rbind(d, cv_merge_rows(dm, section$merge[[m]]))
      rm(dm)
    }
    cells <- map(sort(unique(d$cell)), \(cl)
      cv_one_cell(d[d$cell == cl, ], section$spec, miss_max, selectors))
    tag <- \(x) mutate(x, part = name, .before = 1)
    list(checks = checks,
         perf = tag(list_rbind(map(cells, "perf"))),
         cells = tag(list_rbind(map(cells, "cells"))),
         picks = tag(list_rbind(map(cells, "picks"))),
         pairs = tag(list_rbind(map(cells, "pairs"))),
         sums = tag(list_rbind(map(cells, "sums"))),
         lamgrid = tag(lamgrid))
  })
  lamgrid <- list_rbind(map(out, "lamgrid"))
  # file relative to info$dir as read_summary() in thesis.R checks it
  files <- tibble(part = names(section$files), file = unname(section$files),
                  reps = unname(section$reps[names(section$files)])) |>
    mutate(path = here::here("results", "cv", "our dgp", file),
           size = file.size(path), mtime = file.mtime(path),
           md5 = unname(tools::md5sum(path))) |>
    select(-path)
  list(
    info = list(id = id, title = section$title, spec = section$spec,
                dir = file.path("cv", "our dgp"), files = files,
                selectors = selectors, path_points = cv_path_k,
                pairs = cv_pairs(selectors), interior_at = cv_interior_at,
                miss_max = miss_max, tau = 0,
                built = Sys.time(), r_version = R.version.string,
                packages = map_chr(c("dplyr", "tidyr", "data.table"),
                                   \(p) as.character(packageVersion(p)))),
    checks = list_rbind(map(out, "checks")),
    perf = list_rbind(map(out, "perf")) |>
      left_join(lamgrid, by = c("part", "cell", "k")),
    cells = list_rbind(map(out, "cells")) |>
      left_join(rename(lamgrid, k_opt = k, lambda_opt = lambda),
                by = c("part", "cell", "k_opt")),
    picks = list_rbind(map(out, "picks")),
    lamgrid = lamgrid,
    # path pairs get their point's lambda like the path rows of perf
    pairs = list_rbind(map(out, "pairs")) |>
      mutate(k = suppressWarnings(as.integer(sub("^path_", "", A)))) |>
      left_join(lamgrid, by = c("part", "cell", "k")),
    sums = list_rbind(map(out, "sums"))
  )
}

# wong-chan vs dgp2 (wc_compare). Wong and Chan 2018 design as Wang and Zubizarreta 2020 ran it
# extra study 2 check. balnet only vs dgp2 misspec cv grid (read from its stored summary)

# matched on true ps logit sd (wc 1.005 x 1 2 4. dgp2 1 / c w/ c = 1 .5 .25) and R2_e
# R2_e = share of y(0) var explained by the true ps. each design keeps its own endpoint rule

wc_design <- tidyr::expand_grid(sigma = c(1, 10, 30, 100), mult = c(1, 2, 4))

wc_sections <- list(
  wc_compare = list(
    title = "Wong-Chan against DGP 2 at matched overlap and confounding",
    files = set_names(sprintf("wc_noise_s%g_c%g.rds", wc_design$sigma,
                              wc_design$mult),
                      sprintf("s%g_c%g", wc_design$sigma, wc_design$mult)),
    reps = 1000,
    reached = 0.95,          # grid point counts once this share of reps reaches it
    dgp2 = "cv_misspec",     # stored summary of the twin grid
    overlaps = c(1, 0.5, 0.25))
)

# the two study 2 rules as s2_cell_flags() in thesis.R applies them
wc_interior_gain <- 0.05
wc_win_z <- 3

#' R2_e of a linear y(0) under a linear ps index w/ iid std normal covariates
#' ps is monotone in the index so Var E[Y(0) | e(X)] = cov^2 of signal and index / var(index)
wc_r2e_linear <- function(a, b, noise_var, scale = 1) {
  explained <- scale^2 * sum(a * b)^2 / sum(b^2)
  explained / (scale^2 * sum(a^2) + noise_var)
}

#' R2_e of wc model B = Z1 Z2^3 Z3^2 Z4 + Z4 |Z1|^(1/2). 1st term has mean 0 given Z1 and Z4
#' integrates over u = Z1 + 0.1 Z4 and v = 0.1 Z1 - Z4 (indep). Var B = 45 + sqrt(2 / pi)
wc_r2e_b <- function(sigma) {
  s <- sqrt(1.01)
  cond_mean <- function(u) {
    vapply(u, \(ui) {
      integrate(\(v) (0.1 * ui - v) / 1.01 * abs((ui + 0.1 * v) / 1.01)^0.5 *
                  dnorm(v, sd = s), -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1))
  }
  explained <- integrate(\(u) cond_mean(u)^2 * dnorm(u, sd = s), -Inf, Inf,
                         rel.tol = 1e-10)$value
  explained / (45 + sqrt(2 / pi) + sigma^2)
}

#' R2_e of each design point. wc model A per sigma and model B and dgp2 linear per kappa
wc_r2e <- function(series, level) {
  a_wc <- -0.5 * c(27.4, 13.7, 13.7, 13.7)
  b_wc <- c(1, 0, 0, 0.1)
  a_d2 <- c(1, 0.5, 0.5, 0.5)
  b_d2 <- 1 / (1:4)
  out <- numeric(length(series))
  for (i in seq_along(series)) {
    out[i] <- switch(series[i],
      wc_a = wc_r2e_linear(a_wc, b_wc, level[i]^2),
      wc_b = wc_r2e_b(level[i]),
      dgp2 = wc_r2e_linear(a_d2, b_d2, 1, scale = level[i]))
  }
  out
}

#' one wc file in cv_long() form. truth is the sample ATT for model A and 0 for B
#' early stopping paths repeat their last est so points < `reached` of reps reach are dropped
wc_read <- function(file, reached = 0.95) {
  x <- readRDS(here::here("results", "cv", "wz_replication", file))
  res <- Filter(\(r) is.null(r$err), x$res)
  m <- length(res)
  lam <- x$lam
  lam_end <- vapply(res, \(r) as.numeric(r$lam_end), numeric(1))
  share <- colMeans(outer(lam_end, lam * (1 + 1e-6), "<="))
  keep <- which(share >= reached)
  k_end <- max(keep)
  satt <- vapply(res, \(r) as.numeric(r$satt), numeric(1))
  sigma <- as.numeric(x$cell$sigma)
  mult <- as.numeric(x$cell$c)
  long <- map(c("A", "B"), \(o) {
    truth <- if (o == "A") satt else 0
    path <- t(vapply(res, \(r) r$est_path[o, ], numeric(length(lam)))) - truth
    sel <- t(vapply(res, \(r) r$est_sel[o, cv_selectors],
                    numeric(length(cv_selectors)))) - truth
    colnames(path) <- sprintf("path_%03d", seq_along(lam))
    colnames(sel) <- cv_selectors
    wide <- cbind(path[, keep, drop = FALSE], fixed = path[, k_end], sel)
    tibble(sim = rep(seq_len(m), times = ncol(wide)),
           arm = rep(colnames(wide), each = m),
           outcome = o,
           error = as.vector(wide))
  }) |>
    list_rbind() |>
    mutate(spec = "wc", level = sprintf("s%g_c%g", sigma, mult),
           sigma = sigma, mult = mult, est = "wc",
           estimator = paste0("wc_", arm),
           k = suppressWarnings(as.integer(sub("^path_", "", arm))),
           arm = if_else(is.na(k), arm, "path"), .before = 1) |>
    relocate(spec, level, sigma, mult, sim, est, arm, k, estimator, outcome,
             error)
  list(long = long, k_end = k_end, lam = lam, n_rep = as.numeric(x$n_rep),
       failed = length(x$res) - m)
}

#' checks on one wc file before anything is computed from it
wc_check_file <- function(r, reps, name) {
  per_arm <- r$long |>
    summarise(n = n(), sims = n_distinct(sim), .by = c(outcome, estimator))
  bind_rows(
    check_row(sprintf("%s: %d replications of every arm, none failed", name,
                      reps),
              r$n_rep == reps && r$failed == 0 &&
                all(per_arm$n == reps & per_arm$sims == reps),
              sprintf("%d run, %d failed, %d arms per outcome", r$n_rep,
                      r$failed, n_distinct(per_arm$estimator))),
    check_row(sprintf("%s: path endpoint reached by 95%% of replications",
                      name),
              r$k_end >= 1 && r$k_end <= length(r$lam),
              sprintf("endpoint at grid point %d of %d", r$k_end,
                      length(r$lam))),
    check_row(sprintf("%s: the five shared selectors present", name),
              all(paste0("wc_", cv_selectors) %in% per_arm$estimator),
              paste(cv_selectors, collapse = ", ")))
}

#' each kept grid point and selector vs the path end on the same draws (R_c and mcse)
wc_pairs <- function(long, miss_max = 0.01) {
  arms <- long |>
    filter(estimator != "wc_fixed") |>
    distinct(estimator, k)
  map2(arms$estimator, arms$k, \(a, k) {
    pd <- pair_draws(filter(long, estimator %in% c(a, "wc_fixed")), a,
                     "wc_fixed", miss_max)
    pair_cell_measures(pd$draws) |>
      transmute(level, outcome, A = sub("^wc_", "", a), k = k, reps, rmse_A,
                rmse_B, R_c, log_R_c_mcse)
  }) |>
    list_rbind()
}

#' one point per level and outcome w/ the same rules for both designs (log ratios to the end)
#' best fixed lambda = lowest rmse grid point or the end if nothing beats it (1 + 1e-12)
wc_points <- function(perf, pairs, ends) {
  keys <- c("level", "outcome")
  end <- perf |>
    filter(arm == "fixed") |>
    select(all_of(keys), m, rmse_end = rmse, mcse_end = mcse_rmse)
  opt <- perf |>
    filter(arm == "path") |>
    arrange(k) |>
    slice_min(rmse, n = 1, with_ties = FALSE, by = all_of(keys)) |>
    select(all_of(keys), k_opt = k, path_min = rmse)
  gain <- pairs |>
    filter(!is.na(k)) |>
    inner_join(select(opt, all_of(keys), k = k_opt), by = c(keys, "k")) |>
    transmute(across(all_of(keys)), log_gain = log(R_c),
              mcse_gain = log_R_c_mcse)
  best <- perf |>
    filter(arm %in% cv_selectors) |>
    mutate(arm = factor(arm, cv_selectors)) |>
    arrange(arm) |>
    inner_join(end, by = keys) |>
    mutate(wins = rmse < rmse_end - wc_win_z * sqrt(mcse_rmse^2 + mcse_end^2)) |>
    summarise(best = as.character(arm[which.min(rmse)]), win = any(wins),
              .by = all_of(keys))
  sel <- pairs |>
    filter(is.na(k)) |>
    inner_join(select(best, all_of(keys), A = best), by = c(keys, "A")) |>
    transmute(across(all_of(keys)), log_sel = log(R_c),
              mcse_sel = log_R_c_mcse)
  # left join. no gain pair if the end is the best fixed lambda and the last grid point (dgp2)
  out <- end |>
    inner_join(opt, by = keys) |>
    left_join(gain, by = keys) |>
    inner_join(best, by = keys) |>
    inner_join(sel, by = keys) |>
    mutate(k_end = unname(ends[level]),
           at_end = rmse_end <= path_min * (1 + 1e-12),
           where = case_when(at_end ~ "endpoint", k_opt == 1 ~ "lambda_max",
                             .default = "interior"),
           k_opt = if_else(at_end, as.integer(k_end), as.integer(k_opt)),
           log_gain = if_else(at_end, 0, log_gain),
           mcse_gain = if_else(at_end, 0, mcse_gain),
           interior = path_min < (1 - wc_interior_gain) * rmse_end)
  stopifnot(nrow(out) == nrow(end), all(is.finite(out$log_gain)),
            all(is.finite(out$log_sel)))
  select(out, all_of(keys), m, k_end, k_opt, where, log_gain, mcse_gain,
         best, log_sel, mcse_sel, interior, win, rmse_end, mcse_end)
}

#' dgp2 side from the stored cv_misspec summary. ipw-balnet linear outcome at matched overlaps
#' reshaped for wc_points(). path end is fixed so grid point 30 (which repeats it) is the end
wc_dgp2_points <- function(d2, overlaps) {
  level_of <- \(overlap, strength) sprintf("c%g_k%g", overlap, strength)
  perf <- d2$perf |>
    filter(est == "ipw", outcome == "linear", overlap %in% overlaps) |>
    transmute(level = level_of(overlap, strength), outcome, arm, k, m, rmse,
              mcse_rmse, overlap, strength)
  pairs <- d2$pairs |>
    filter(est == "ipw", outcome == "linear", overlap %in% overlaps,
           kind %in% c("path", "selector"), B == "fixed") |>
    transmute(level = level_of(overlap, strength), outcome, A, k, R_c,
              log_R_c_mcse)
  ends <- set_names(rep(d2$info$path_points, n_distinct(perf$level)),
                    unique(perf$level))
  points <- wc_points(perf, pairs, ends) |>
    inner_join(distinct(perf, level, overlap, strength), by = "level")
  list(points = points, perf = perf, pairs = pairs)
}

#' each wc point vs the dgp2 points bracketing its R2_e at the same logit sd. model B vs kappa 0
#' brackets since sigma 30 is ~halfway between kappa 1/4 and 1/2 (log). agree = matches either
wc_twins <- function(points) {
  wc <- filter(points, design != "DGP 2")
  d2 <- filter(points, design == "DGP 2")
  map(seq_len(nrow(wc)), \(i) {
    w <- wc[i, ]
    pool <- filter(d2, sd_logit == w$sd_logit)
    tw <- if (w$series == "wc_b") {
      filter(pool, level == 0)
    } else {
      pos <- filter(pool, r2e > 0)
      bind_rows(slice_max(filter(pos, r2e <= w$r2e), r2e, n = 1),
                slice_min(filter(pos, r2e >= w$r2e), r2e, n = 1)) |>
        distinct()
    }
    tw <- arrange(tw, r2e)
    tibble(sd_logit = w$sd_logit, series = w$series, sigma = w$level,
           r2e = w$r2e,
           kappa_lo = min(tw$level), kappa_hi = max(tw$level),
           r2e_lo = min(tw$r2e), r2e_hi = max(tw$r2e),
           wc_gain = w$log_gain, wc_gain_mcse = w$mcse_gain,
           d2_gain_lo = min(tw$log_gain), d2_gain_hi = max(tw$log_gain),
           d2_gain_mcse = max(tw$mcse_gain),
           wc_sel = w$log_sel, wc_sel_mcse = w$mcse_sel,
           d2_sel_lo = min(tw$log_sel), d2_sel_hi = max(tw$log_sel),
           d2_sel_mcse = max(tw$mcse_sel),
           wc_interior = w$interior, wc_win = w$win,
           d2_interior_lo = tw$interior[1], d2_interior_hi = tw$interior[nrow(tw)],
           d2_win_lo = tw$win[1], d2_win_hi = tw$win[nrow(tw)],
           interior_agree = w$interior %in% tw$interior,
           win_agree = w$win %in% tw$win)
  }) |>
    list_rbind() |>
    arrange(sd_logit, r2e)
}

#' build the wc entry. read and check the 12 files then both designs by one route and match
build_wc_section <- function(section, id, miss_max = 0.01) {
  reads <- imap(section$files, \(file, name) wc_read(file, section$reached))
  checks <- imap(reads, \(r, name) wc_check_file(r, section$reps, name)) |>
    list_rbind()
  long <- list_rbind(map(reads, "long"))
  ends <- set_names(map_dbl(reads, "k_end"),
                    map_chr(reads, \(r) r$long$level[1]))
  perf <- cell_performance(long, miss_max) |>
    left_join(distinct(long, level, estimator, arm, k, sigma, mult),
              by = c("level", "estimator"))
  pairs <- wc_pairs(long, miss_max)
  wc <- wc_points(perf, pairs, ends) |>
    inner_join(distinct(perf, level, sigma, mult), by = "level") |>
    mutate(design = "Wong–Chan",
           series = if_else(outcome == "A", "wc_a", "wc_b"),
           sd_logit = mult, level_value = sigma)
  d2_path <- here::here("output", "thesis_summaries",
                        paste0(section$dgp2, ".rds"))
  d2 <- readRDS(d2_path)
  twin <- wc_dgp2_points(d2, section$overlaps)
  dgp2 <- twin$points |>
    mutate(design = "DGP 2", series = "dgp2", sd_logit = 1 / overlap,
           level_value = strength)
  points <- bind_rows(wc, dgp2) |>
    mutate(level_key = level, level = level_value,
           r2e = wc_r2e(series, level)) |>
    select(design, series, outcome, level, sd_logit, r2e, m, k_end, k_opt,
           where, log_gain, mcse_gain, best, log_sel, mcse_sel, interior,
           win, rmse_end, mcse_end, level_key) |>
    arrange(sd_logit, design, r2e)
  # model B only at sigma = 1 (the paper's design). at larger sigma its R2_e < 1e-4
  compared <- filter(points, series != "wc_b" | level == 1)
  twins <- wc_twins(compared)
  checks <- bind_rows(
    checks,
    check_row("DGP 2 twin: the stored summary passed its own checks",
              isTRUE(d2$info$verified) && all(d2$checks$ok),
              sprintf("%s built %s", section$dgp2,
                      format(d2$info$built, "%Y-%m-%d %H:%M"))),
    check_row("DGP 2 twin: every matched overlap level present",
              setequal(unique(dgp2$overlap), section$overlaps),
              paste(sort(unique(dgp2$overlap)), collapse = ", ")),
    check_row("every Wong-Chan point matched to a DGP 2 bracket",
              nrow(twins) == nrow(filter(compared, design != "DGP 2")) &&
                all(is.finite(twins$r2e_lo)),
              sprintf("%d matched pairs", nrow(twins))))
  files <- tibble(part = names(section$files), file = unname(section$files),
                  reps = section$reps) |>
    mutate(path = here::here("results", "cv", "wz_replication", file),
           size = file.size(path), mtime = file.mtime(path),
           md5 = unname(tools::md5sum(path))) |>
    select(-path)
  list(
    info = list(id = id, title = section$title, dir = file.path("cv", "wz_replication"),
                files = files, reached = section$reached,
                selectors = cv_selectors, overlaps = section$overlaps,
                interior_gain = wc_interior_gain, win_z = wc_win_z,
                dgp2 = list(id = section$dgp2, built = d2$info$built,
                            md5 = unname(tools::md5sum(d2_path))),
                miss_max = miss_max,
                built = Sys.time(), r_version = R.version.string,
                packages = map_chr(c("dplyr", "tidyr"),
                                   \(p) as.character(packageVersion(p)))),
    checks = checks,
    points = points,
    twins = twins,
    perf = perf,
    pairs = pairs
  )
}
