# make_registry.R
# scans results/*.csv.gz, derives the factual columns from the data, joins
# registry.csv and prints anything that doesn't match. run from repo root.
# writes registry_derived.csv. facts come from the files; judgment calls
# (status, notes, dates) stay in registry.csv.

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

results_dir <- "results"
registry_path <- "registry.csv"

perf_est <- c("balnet0", "balnet05", "balnet10", "balnetcv", "balnetrate",
              "glmcv_normal", "glmnetcv", "oracle", "normalised_oracle")

derive_one <- function(path) {
  df <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  fac <- setdiff(names(df), c("sim", "estimator", "tau_hat"))
  ests <- sort(unique(df$estimator))
  reps <- df |>
    group_by(across(all_of(fac))) |>
    summarise(r = n_distinct(sim), .groups = "drop")
  lam <- df |> filter(estimator %in% c("lam_end0", "lam_end1"))
  errdf <- df |> filter(estimator == "err")
  perf_tau <- df$tau_hat[df$estimator %in% perf_est]

  # shared draws: lam_end0 equal across outcomes (only works for DGP2-era files)
  shared <- NA
  if ("outcome" %in% fac && nrow(lam) > 0 && !all(is.na(lam$tau_hat))) {
    key <- setdiff(fac, "outcome")
    w <- df |>
      filter(estimator == "lam_end0") |>
      pivot_wider(id_cols = all_of(c(key, "sim")),
                  names_from = outcome, values_from = tau_hat)
    oc <- setdiff(names(w), c(key, "sim"))
    if (length(oc) > 1) {
      shared <- all(map_lgl(oc[-1], ~ identical(w[[oc[1]]], w[[.x]])))
    }
  }

  lvl <- function(col) {
    if (col %in% fac) paste(sort(unique(df[[col]])), collapse = ";") else NA_character_
  }

  tibble(
    batch_file = basename(path),
    d_rows = nrow(df),
    d_cells = nrow(reps),
    d_reps_min = min(reps$r),
    d_reps_max = max(reps$r),
    d_estimator_rows = length(ests),
    d_n = lvl("n"), d_p = lvl("p"), d_s = lvl("s"),
    d_outcome = lvl("outcome"), d_misspec = lvl("misspec"),
    d_covcor = lvl("covcor"), d_overlap = lvl("overlap"),
    d_signs = lvl("signs"), d_decay = lvl("decay"),
    d_outcome_set = lvl("outcome_set"),
    d_na_frac_all = mean(is.na(df$tau_hat)),
    d_na_frac_perf = if (length(perf_tau)) mean(is.na(perf_tau)) else NA_real_,
    d_err_frac = if (nrow(errdf)) mean(errdf$tau_hat == 1, na.rm = TRUE) else NA_real_,
    d_lambda_logged = nrow(lam) > 0,
    d_lam_end_min = if (nrow(lam) && !all(is.na(lam$tau_hat))) {
      min(lam$tau_hat, na.rm = TRUE)
    } else NA_real_,
    d_shared_draws = shared
  )
}

files <- list.files(results_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
derived <- map_dfr(files, derive_one)

meta <- read_csv(registry_path, show_col_types = FALSE)
reg <- full_join(meta, derived, by = "batch_file")
write_csv(reg, "registry_derived.csv")

# mismatch checks
flags <- character()

no_meta <- setdiff(derived$batch_file, meta$batch_file)
if (length(no_meta)) {
  flags <- c(flags, paste("file not in registry:", no_meta))
}
planned <- meta$batch_file[str_detect(coalesce(meta$status, ""), "planned")]
no_file <- setdiff(setdiff(meta$batch_file, derived$batch_file), planned)
if (length(no_file)) {
  flags <- c(flags, paste("registry row with no file (not planned):", no_file))
}

bad <- reg |> filter(!is.na(d_err_frac) & d_err_frac > 0 | coalesce(d_na_frac_perf, 0) > 0)
if (nrow(bad)) {
  flags <- c(flags, sprintf("err/NA in %s (err_frac=%.3f perf_na=%.3f)",
                            bad$batch_file, coalesce(bad$d_err_frac, 0),
                            coalesce(bad$d_na_frac_perf, 0)))
}

reps_meta <- suppressWarnings(as.numeric(reg$reps_per_cell))
rm_idx <- which(!is.na(reps_meta) & !is.na(reg$d_reps_min) &
                (reps_meta != reg$d_reps_min | reg$d_reps_min != reg$d_reps_max))
if (length(rm_idx)) {
  flags <- c(flags, paste("reps mismatch:", reg$batch_file[rm_idx]))
}

floor_meta <- suppressWarnings(as.numeric(str_extract(reg$floor_max_imbalance,
                                                      "[0-9.]+e?-?[0-9]*")))
fm_idx <- which(!is.na(floor_meta) & !is.na(reg$d_lam_end_min) &
                abs(floor_meta - reg$d_lam_end_min) > 1e-12)
if (length(fm_idx)) {
  flags <- c(flags, sprintf("floor mismatch: %s (recorded %s, observed %g)",
                            reg$batch_file[fm_idx],
                            reg$floor_max_imbalance[fm_idx],
                            reg$d_lam_end_min[fm_idx]))
}

sign_set <- function(x) sort(strsplit(str_remove(x, " \\(.*\\)$"), ";")[[1]])
sm_idx <- which(map_lgl(seq_len(nrow(reg)), function(i) {
  a <- reg$signs[i]; b <- reg$d_signs[i]
  if (is.na(a) || is.na(b) || str_detect(a, "GAP|planned|inferred|FINGERPRINT")) return(FALSE)
  !setequal(sign_set(a), sign_set(b))
}))
if (length(sm_idx)) {
  flags <- c(flags, paste("signs mismatch:", reg$batch_file[sm_idx]))
}

cat(if (length(flags)) paste(flags, collapse = "\n") else "no mismatches", "\n")

# seed sharing check: oracle rows depend only on the draws, so exact equality
# between two files means shared seeds for those cells
check_seed_sharing <- function(file_a, file_b,
                               keys = c("n", "p", "outcome", "overlap", "misspec"),
                               est = c("oracle", "normalised_oracle")) {
  a <- read_csv(file.path(results_dir, file_a), show_col_types = FALSE) |>
    filter(estimator %in% est)
  b <- read_csv(file.path(results_dir, file_b), show_col_types = FALSE) |>
    filter(estimator %in% est)
  keys <- intersect(keys, intersect(names(a), names(b)))
  m <- inner_join(a, b, by = c(keys, "sim", "estimator"), suffix = c("_a", "_b"))
  if (!nrow(m)) return("no shared cells")
  m |>
    mutate(eq = tau_hat_a == tau_hat_b) |>
    group_by(across(all_of(keys))) |>
    summarise(frac_equal = mean(eq), n_rows = n(), .groups = "drop") |>
    arrange(desc(frac_equal))
}
# Example: check_seed_sharing("E1_results.csv.gz", "E1_resultsREDUX.csv.gz")
