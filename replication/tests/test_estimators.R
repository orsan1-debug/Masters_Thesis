# estimate_study1() vs estimate_all() (R/estimators_ipw.R) on the same draw row by row. no timings
# then what the thesis relies on (oracle unbiased. path end. Hajek = HT. ABW near OLS. MLE shrinks)

library(testthat)
source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))

record <- new.env()
sys.source(here::here("R", "estimators_ipw.R"), envir = record)

# part A. same numbers as the record

as_record_draw <- function(d) {
  record_names <- c(linear = "linear", quadratic = "quad1", exponential = "exp")
  colnames(d$Y) <- colnames(d$m) <- unname(record_names[colnames(d$Y)])
  d
}
timing_rows <- c("time_bal", "time_ps", "time_or")

#' seed after the dots so s cannot partially match it
expect_same_estimates <- function(..., seed) {
  set.seed(seed); d <- generate_sample(...)
  set.seed(seed + 1); new <- estimate_study1(d)
  set.seed(seed + 1); old <- record$estimate_all(as_record_draw(d))
  colnames(old) <- colnames(new)
  keep <- setdiff(rownames(old), timing_rows)
  expect_identical(rownames(new), rownames(old))
  expect_equal(new[keep, ], old[keep, ], tolerance = 1e-12)
}

test_that("correct specification, base case", {
  expect_same_estimates(n = 400, p = 20, seed = 101)
})

test_that("misspecified score: the refit-on-X_true rows included", {
  expect_same_estimates(n = 400, p = 20, misspecification = "kang_schafer", seed = 102)
})

test_that("weak overlap and a noisy outcome", {
  expect_same_estimates(n = 400, p = 20, c = 0.25, kappa = 0.05, seed = 103)
})

test_that("p past n", {
  expect_same_estimates(n = 150, p = 300, outcomes = c("linear", "quadratic"), seed = 104)
})

test_that("denser profile under Tan's transformation, AR(1) covariates", {
  expect_same_estimates(n = 300, p = 60, s = 16, nu = 0.5, covariance = "ar1",
                        loadings = "all_active", misspecification = "tan",
                        outcomes = c("linear", "quadratic"), seed = 105)
})

# part B. what the thesis relies on

test_that("the Oracle recovers tau = 0 over replications, both forms", {
  reps <- 200
  est <- replicate(reps, {
    d <- generate_sample(500, p = 10, outcomes = "linear")
    c(hajek = unname(ate_hajek(d$Y, d$W, d$e, d$e0)),
      ht    = unname(ate_horvitz_thompson(d$Y, d$W, d$e, d$e0)))
  })
  for (form in c("hajek", "ht"))
    expect_lt(abs(mean(est[form, ])), 3 * sd(est[form, ]) / sqrt(reps))
})

test_that("at p/n = 0.05 the balnet path reaches its target and balances every covariate", {
  set.seed(201)
  d   <- generate_sample(1000, p = 50)
  bal <- fit_balancing(d$X, d$W)
  expect_true(all(bal$lambda_end <= 1e-4 + 1e-12))
  xbar <- colMeans(d$X)
  sdx  <- apply(d$X, 2, sd) * sqrt(999 / 1000)
  for (arm in c("treated", "control"))
    expect_lt(weight_diagnostics(bal$weights[[arm]], d$X, xbar, sdx)[["smd"]], 1.5e-4)
})

test_that("from p/n = 0.5 the path stops short of its target", {
  set.seed(202)
  d   <- generate_sample(400, p = 200)
  bal <- fit_balancing(d$X, d$W)
  expect_true(all(bal$lambda_end > 1e-4))
})

test_that("Hajek and Horvitz-Thompson coincide for balancing weights", {
  set.seed(203)
  d   <- generate_sample(1000, p = 50)
  w   <- fit_balancing(d$X, d$W)$weights
  expect_equal(sum(w$treated), 1000, tolerance = 1e-8)
  expect_equal(sum(w$control), 1000, tolerance = 1e-8)
  expect_equal(colMeans(d$Y * drop(w$treated - w$control)),
               ate_weighted(d$Y, w$treated, w$control), tolerance = 1e-8)
})

test_that("at exact balance ABW is close to OLS, the collapse (Bruns-Smith et al., 2026)", {
  set.seed(204)
  d   <- generate_sample(2000, p = 20, outcomes = "linear")
  bal <- fit_balancing(d$X, d$W)
  om  <- fit_outcome_models(d$X, d$Y, d$W)
  abw <- ate_augmented(d$Y, bal$weights$treated, bal$weights$control, om$m1, om$m0)
  Xc  <- cbind(1, d$X)
  b   <- function(arm) lm.fit(Xc[d$W == arm, ], d$Y[d$W == arm, 1])$coefficients
  ols <- drop(colMeans(Xc) %*% (b(1) - b(0)))
  expect_lt(abs(abw - ols), 0.05)
  expect_lt(abs(abw - ols), abs(ate_weighted(d$Y, bal$weights$treated, bal$weights$control) - ols) + 0.05)
})

test_that("the MLE weights carry the lasso's shrinkage: e_hat is pulled toward one half", {
  set.seed(205)
  d   <- generate_sample(1000, p = 50, c = 1)
  mle <- fit_mle(d$X, d$W)
  expect_lt(sd(qlogis(mle$e_hat)), sd(d$logit))     # flatter fitted logit than the truth
  expect_gt(cor(mle$e_hat, d$e), 0.9)                # but the right ordering
})
