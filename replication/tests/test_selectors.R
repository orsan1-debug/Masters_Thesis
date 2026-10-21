# estimate_study2() and pilot_lambda_grid() vs the record (R/estimators_cv_dgp2.R) on one draw
# lasso and enet. also last grid point = endpoint and selected lambdas lie on the path

library(testthat)
source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))
source(here::here("replication", "R", "selectors.R"))

record <- new.env()
sys.source(here::here("R", "estimators_ipw.R"), envir = record)
sys.source(here::here("R", "estimators_cv_dgp2.R"), envir = record)

#' the record names outcomes linear quad1 exp so its strength cols read linear_s1 quad1_s1 exp_s1
as_record_draw <- function(d) {
  to_record <- function(x) x |>
    sub("^quadratic", "quad1", x = _) |>
    sub("^exponential", "exp", x = _)
  colnames(d$Y) <- to_record(colnames(d$Y))
  colnames(d$m) <- to_record(colnames(d$m))
  d
}

expect_same_study2 <- function(alpha, seed) {
  set.seed(seed)
  d <- share_strengths(generate_sample(500, p = 30, c = 0.5))
  settings <- path_settings(alpha)
  min_ratio <- if (alpha == 1) 1e-4 else settings$lambda.min.ratio
  grid <- pilot_lambda_grid(d, alpha = alpha, min_ratio = min_ratio)
  grid_record <- record$make_lam_grid(d, alpha = alpha, min_ratio = min_ratio)
  expect_equal(grid, grid_record, tolerance = 1e-12)

  set.seed(seed + 1)
  new <- estimate_study2(d, grid, alpha = alpha, min_ratio = min_ratio, B = 100)
  set.seed(seed + 1)
  old <- record$estimate_cv(as_record_draw(d), lam_grid = grid, alpha = alpha,
                            min_ratio = min_ratio, B = 100)
  colnames(old) <- colnames(new)
  expect_identical(rownames(new), rownames(old))
  expect_equal(new, old, tolerance = 1e-12)
  list(draw = d, grid = grid, estimates = new)
}

test_that("lasso: same grid, selected lambdas, estimates and diagnostics as the record", {
  r <- expect_same_study2(alpha = 1, seed = 301)
  est <- r$estimates
  expect_equal(est["ipw_path_30", ], est["ipw_fixed", ], tolerance = 1e-10)
  expect_equal(est["abw_path_30", ], est["abw_fixed", ], tolerance = 1e-10)
  for (sel in selector_names) for (arm in c("lam1_", "lam0_")) {
    lam <- est[paste0(arm, sel), 1]
    expect_gte(lam, est[paste0(arm, "fixed"), 1] - 1e-12)
    expect_lte(lam, max(r$grid) + 1e-12)
  }
})

test_that("elastic net (alpha = 0.5): same as the record, path ending at a fraction of lambda_max", {
  r <- expect_same_study2(alpha = 0.5, seed = 302)
  expect_equal(r$grid[length(r$grid)], 0)
})

test_that("share_strengths() keeps the noise and scales the signal", {
  set.seed(303)
  d <- generate_sample(300, p = 10, outcomes = c("linear", "quadratic"))
  s <- share_strengths(d, strengths = c(1, 0.5, 0))
  expect_equal(colnames(s$Y), c("linear_s1", "quadratic_s1", "linear_s0.5", "quadratic_s0.5",
                                "linear_s0", "quadratic_s0"))
  eps <- d$Y - d$m
  expect_equal(unname(s$Y[, "linear_s0"]), unname(eps[, "linear"]), tolerance = 1e-12)
  expect_equal(unname(s$Y[, "linear_s0.5"]), unname(0.5 * d$m[, "linear"] + eps[, "linear"]),
               tolerance = 1e-12)
  expect_equal(s$Y[, "quadratic_s1"], d$Y[, "quadratic"], tolerance = 1e-12)
})
