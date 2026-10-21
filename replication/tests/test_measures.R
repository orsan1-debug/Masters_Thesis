# measures recomputed from the raw sample size results = stored misspec_n correct_n and coverage
# skipped when the results files or the stored summaries are absent

library(testthat)
source(here::here("replication", "R", "measures.R"))

stored_outcome <- c(linear = "linear", quadratic = "quad1", exponential = "exp")
results_file <- function(f) here::here("results", "aipw", f)

#' main file + add-on file. the add-on oracle_hajek rows repeat the main file so dropped
load_spec <- function(main, add, spec) {
  if (!file.exists(results_file(main)) || !file.exists(results_file(add)))
    skip(paste("results file not present:", main))
  bind_rows(read_results(results_file(main)),
            read_results(results_file(add)) |> filter(estimator != "oracle_hajek")) |>
    mutate(spec = spec)
}

check_against_stored <- function(results, summary_id, pairs) {
  path <- here::here("output", "thesis_summaries", paste0(summary_id, ".rds"))
  if (!file.exists(path)) skip(paste("stored summary not present:", summary_id))
  stored <- readRDS(path)

  perf <- performance(results) |> mutate(outcome = stored_outcome[outcome])
  joined <- inner_join(perf, as_tibble(stored$perf),
                       by = c("spec", "level", "estimator", "outcome"),
                       suffix = c("", ".stored"))
  expect_equal(nrow(joined), nrow(stored$perf))
  for (v in c("m", "missing", "bias", "empse", "rmse", "mcse_bias", "mcse_empse",
              "mcse_rmse", "share", "mcse_share", "dropped"))
    expect_equal(joined[[v]], joined[[paste0(v, ".stored")]], tolerance = 1e-10,
                 label = paste(summary_id, v))

  for (i in seq_len(nrow(pairs))) {
    A <- pairs$A[i]; B <- pairs$B[i]
    pd    <- pair_draws(results, A, B)
    cells <- pair_cells(pd$draws) |> mutate(outcome = stored_outcome[outcome])
    sc <- as_tibble(stored$cells) |> filter(A == !!A, B == !!B, kept)
    j <- inner_join(cells, sc, by = c("spec", "level", "outcome"), suffix = c("", ".stored"))
    expect_equal(nrow(j), nrow(sc), label = paste(summary_id, A, B, "cells"))
    for (v in c("reps", "rmse_A", "rmse_B", "R_c", "log_R_c_mcse", "s_A", "s_A_mcse",
                "s_B", "s_B_mcse", "won", "lost"))
      expect_equal(j[[v]], j[[paste0(v, ".stored")]], tolerance = 1e-10,
                   label = paste(summary_id, A, B, v))

    avg <- relative_rmse(results, A, B) |>
      mutate(outcome = ifelse(outcome == "pooled", "pooled", stored_outcome[outcome]))
    sa <- as_tibble(stored$avg) |>
      filter(A == !!A, B == !!B, theta == "pair", set == "all", is.na(factor))
    j <- inner_join(avg, sa, by = c("block", "outcome"), suffix = c("", ".stored"))
    expect_equal(nrow(j), 4L, label = paste(summary_id, A, B, "averages"))
    for (v in c("n_cells", "R_bar", "R_bar_mcse", "R_min", "R_max", "won", "lost",
                "s_bar_A", "s_bar_A_mcse", "s_bar_B", "s_bar_B_mcse"))
      expect_equal(j[[v]], j[[paste0(v, ".stored")]], tolerance = 1e-10,
                   label = paste(summary_id, A, B, v))
  }
  invisible(perf)
}

pairs_main <- tibble::tribble(
  ~A,              ~B,
  "balnet0",       "glmnetcv_hajek",
  "abw_balnet0",   "aipw_glmnetcv",
  "aipw_glmnetcv", "glmnetcv_hajek",
  "balnet0",       "oracle_hajek",
  "balnet0",       "ols")

test_that("misspecified sample-size grid: every stored number is reproduced", {
  results <- load_spec("aipw_misspec_n_dgp2.csv.gz", "aipw_add_misspec_n_dgp2.csv.gz", "misspec")
  check_against_stored(results, "misspec_n",
                       bind_rows(pairs_main,
                                 tibble(A = c("abw_balnet0_ortrue", "abw_balnet0"),
                                        B = c("aipw_glmnetcv_ortrue", "balnet0"))))
})

test_that("correctly specified sample-size grid: every stored number is reproduced", {
  results <- load_spec("aipw_correct_n400k_dgp2.csv.gz", "aipw_add_correct_n400k_dgp2.csv.gz",
                       "correct")
  check_against_stored(results, "correct_n", pairs_main)
})

test_that("coverage at the empirical SE matches the stored table", {
  path <- here::here("output", "thesis_summaries", "coverage_empse.csv")
  if (!file.exists(path)) skip("stored coverage table not present")
  stored <- read.csv(path) |>
    filter(grid == "n", spec == "misspec")
  results <- load_spec("aipw_misspec_n_dgp2.csv.gz", "aipw_add_misspec_n_dgp2.csv.gz", "misspec")
  perf <- performance(results) |>
    mutate(outcome = stored_outcome[outcome]) |>
    semi_join(stored, by = c("n", "overlap", "estimator", "outcome"))
  j <- inner_join(perf, stored, by = c("n", "overlap", "estimator", "outcome"),
                  suffix = c("", ".stored"))
  expect_gt(nrow(j), 0)
  expect_equal(j$coverage, j$coverage.stored, tolerance = 1e-10)
  expect_equal(j$mcse_coverage, j$mcse, tolerance = 1e-10)
  expect_equal(j$bias, j$bias.stored, tolerance = 1e-10)
  expect_equal(j$empse, j$empse.stored, tolerance = 1e-10)
})

test_that("the measures are right on a case with a known answer", {
  m <- cell_measures(c(1, -1, 3, -3))
  expect_equal(m$bias, 0)
  expect_equal(m$rmse, sqrt(5))
  expect_equal(m$empse, sqrt(20 / 3))
  expect_equal(m$coverage, 1)
  expect_equal(m$share, 0)
  m <- cell_measures(rep(2, 10))
  expect_equal(m$rmse, 2); expect_equal(m$share, 1); expect_equal(m$coverage, 0)
  r <- tibble(spec = "s", level = rep(c("a", "b"), each = 400), sim = rep(1:400, 2),
              outcome = "linear") |>
    mutate(balnet0 = rnorm(800, sd = 1), other = 2 * balnet0) |>
    pivot_longer(c(balnet0, other), names_to = "estimator", values_to = "error")
  expect_equal(relative_rmse(r, "balnet0", "other")$R_bar[1], 50, tolerance = 1e-10)
})
