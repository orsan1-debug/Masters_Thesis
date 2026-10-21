# generate_wong_chan() and estimate_wong_chan() vs dgp_wc_overlap() and wc_estimate() on one seed
# plus a stored rep regenerated. w/o results/ the stored cell is read from tests/fixtures/

library(testthat)
source(here::here("replication", "R", "wong_chan.R"))

record <- new.env()
sys.source(here::here("R", "dgp.R"), envir = record)
sys.source(here::here("R", "estimators_cv.R"), envir = record)
sys.source(here::here("runs", "cv", "wz_replication", "path_floor.R"), envir = record)
sys.source(here::here("runs", "cv", "wz_replication", "wc_helpers.R"), envir = record)

#' a stored file or if results/ is absent its rows in tests/fixtures/ (from make_fixtures.R)
stored_path <- function(...) {
  path <- here::here("results", ...)
  ifelse(file.exists(path), path, here::here("replication", "tests", "fixtures", ...))
}

test_that("the draw is the record's, at every overlap and noise level", {
  for (mult in c(1, 2, 4)) for (sigma in c(1, 30)) {
    set.seed(401); new <- generate_wong_chan(800, overlap = mult, sigma = sigma)
    set.seed(401); old <- record$dgp_wc_overlap(800, overlap = mult, sigma = sigma)
    names(old)[names(old) == "eta"] <- "logit"
    expect_identical(new[names(old)], old)
  }
})

test_that("one replication of every arm equals the record's", {
  set.seed(402); d <- generate_wong_chan(2000, overlap = 2, sigma = 10)
  basis <- wc_basis(d$X)
  grid_new <- wc_lambda_grid(basis, d$W)
  grid_old <- record$wc_lam_grid(basis, d$W)
  expect_equal(grid_new, grid_old, tolerance = 1e-12)
  tolerances <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  set.seed(403); new <- estimate_wong_chan(d, basis, grid_new, tolerances)
  set.seed(403); old <- record$wc_estimate(d, basis, grid_old, tolerances)
  expect_identical(names(new), names(old))
  expect_equal(new, old, tolerance = 1e-12)
})

test_that("design facts: logit sd about 1.005 times the multiplier, half treated, no effect in B", {
  set.seed(404); d <- generate_wong_chan(2e5, overlap = 2)
  expect_equal(sd(d$logit), 2 * sqrt(1 + 0.01), tolerance = 0.01)
  expect_equal(mean(d$e), 0.5, tolerance = 0.005)
  # outcome B has no effect so its noise has the same mean in both arms
  Z <- d$X_true
  eps_B <- d$Y[, "B"] - (Z[, 1] * Z[, 2]^3 * Z[, 3]^2 * Z[, 4] + Z[, 4] * abs(Z[, 1])^0.5)
  expect_lt(abs(mean(eps_B[d$W == 1]) - mean(eps_B[d$W == 0])), 0.02)
  g <- 27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]
  expect_equal(d$Y[, "A"] - eps_B, 210 + (1.5 * d$W - 0.5) * g, tolerance = 1e-10)
  # outcome A's sample ATT = mean of the individual effects on the treated
  expect_equal(mean(d$tau_i[d$W == 1]), mean(1.5 * (27.4 * d$X_true[d$W == 1, 1] +
               13.7 * rowSums(d$X_true[d$W == 1, 2:4]))), tolerance = 1e-10)
})

test_that("the first replication of a stored cell is regenerated: sigma = 10, c = 2", {
  f <- stored_path("cv", "wz_replication", "wc_noise_s10_c2.rds")
  if (!file.exists(f)) skip("stored cell not present")
  source(here::here("replication", "R", "simulate.R"))
  stored <- readRDS(f)
  streams <- seed_streams(wc_master_seed, 1)
  assign(".Random.seed", streams[[1]], envir = globalenv())
  d <- generate_wong_chan(5000, overlap = 2, sigma = 10)
  basis <- wc_basis(d$X)
  lambda_grid <- wc_lambda_grid(basis, d$W)
  expect_equal(lambda_grid, stored$lam, tolerance = 1e-10)
  assign(".Random.seed", streams[[1]], envir = globalenv())
  d <- generate_wong_chan(5000, overlap = 2, sigma = 10)
  new <- estimate_wong_chan(d, basis, lambda_grid, stored$grid)
  old <- stored$res[[1]]
  expect_equal(new$satt, old$satt, tolerance = 1e-10)
  expect_equal(new$est_path, old$est_path, tolerance = 1e-8)
  expect_equal(new$lam_end, old$lam_end, tolerance = 1e-10)
  expect_equal(new$est_sel[, c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")],
               old$est_sel[, c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")],
               tolerance = 1e-8)
  expect_equal(new$lam_sel[1:5], old$lam_sel[1:5], tolerance = 1e-10)
})
