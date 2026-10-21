# generate_sample() vs R/dgp.R (dgp2 dgp3) draw for draw on every setting type the thesis runs
# then the design facts the thesis states (logit sd 1/c share 0.5 App Tables 6 and 20 tau = 0)

library(testthat)
source(here::here("replication", "R", "dgp.R"))

record <- new.env()
sys.source(here::here("R", "dgp.R"), envir = record)

# part A. same draws as the record

#' the record names the outcomes linear quad1 exp and the logit eta
as_thesis_names <- function(d) {
  thesis_names <- c(linear = "linear", quad1 = "quadratic", exp = "exponential")
  colnames(d$Y) <- colnames(d$m) <- unname(thesis_names[colnames(d$Y)])
  names(d)[names(d) == "eta"] <- "logit"
  d[c("Y", "W", "X", "X_true", "e", "e0", "logit", "m", "tau")]
}

expect_same_draws <- function(seed, new_args, old_fun, old_args) {
  set.seed(seed); new <- do.call(generate_sample, new_args)
  set.seed(seed); old <- as_thesis_names(do.call(old_fun, old_args))
  expect_identical(new[names(old)], old)
}

test_that("sample size, overlap and signal grids: dgp2, four confounders, both specifications", {
  for (cov in c("independent", "ar1")) for (mis in c(FALSE, TRUE)) for (c in c(0.15, 1))
    expect_same_draws(
      seed = 11,
      new_args = list(n = 400, p = 50, s = 4, nu = 1, covariance = cov, c = c,
                      kappa = 0.35,
                      misspecification = if (mis) "kang_schafer" else "none"),
      old_fun  = record$dgp2,
      old_args = list(n = 400, p = 50, s = 4, outcome = c("linear", "quad1", "exp"),
                      misspec = mis, covcor = if (cov == "ar1") "ar1" else "iid",
                      overlap = c, signs = "pos", decay_ps = 1, decay_out = 1,
                      outcome_set = "fixed4", strength = 0.35))
})

test_that("dimension grid: p past n", {
  expect_same_draws(
    seed = 12,
    new_args = list(n = 100, p = 400, misspecification = "kang_schafer",
                    outcomes = c("linear", "quadratic")),
    old_fun  = record$dgp2,
    old_args = list(n = 100, p = 400, outcome = c("linear", "quad1"), misspec = TRUE,
                    covcor = "iid", signs = "pos"))
})

test_that("denser grids, correct specification: dgp2 with loadings on all active covariates", {
  for (cov in c("independent", "ar1")) for (s in c(4, 16, 64)) for (nu in c(1, 0.5))
    expect_same_draws(
      seed = 13,
      new_args = list(n = 300, p = 150, s = s, nu = nu, covariance = cov,
                      outcomes = c("linear", "quadratic"), loadings = "all_active"),
      old_fun  = record$dgp2,
      old_args = list(n = 300, p = 150, s = s, outcome = c("linear", "quad1"),
                      covcor = if (cov == "ar1") "ar1" else "iid", signs = "pos",
                      decay_ps = nu, decay_out = nu, outcome_set = "track_s"))
})

test_that("denser grids, misspecified: dgp3 with Tan's transformation", {
  for (cov in c("independent", "ar1")) for (s in c(4, 64)) for (nu in c(1, 0.5))
    expect_same_draws(
      seed = 14,
      new_args = list(n = 300, p = 150, s = s, nu = nu, covariance = cov,
                      outcomes = c("linear", "quadratic"), loadings = "all_active",
                      misspecification = "tan"),
      old_fun  = record$dgp3,
      old_args = list(n = 300, p = 150, s = s, outcome = c("linear", "quad1"),
                      misspec = TRUE, misspec_form = "tan_obs",
                      covcor = if (cov == "ar1") "ar1" else "iid", signs = "pos",
                      decay_ps = nu, decay_out = nu, outcome_set = "track_s"))
})

test_that("denser grids, equicorrelated covariates: dgp3, both specifications", {
  for (mis in c(FALSE, TRUE)) for (s in c(4, 16)) for (nu in c(1, 0.5))
    expect_same_draws(
      seed = 15,
      new_args = list(n = 300, p = 150, s = s, nu = nu, covariance = "equicorrelated",
                      outcomes = c("linear", "quadratic"), loadings = "all_active",
                      misspecification = if (mis) "tan" else "none"),
      old_fun  = record$dgp3,
      old_args = list(n = 300, p = 150, s = s, outcome = c("linear", "quad1"),
                      misspec = mis, misspec_form = "tan_obs", covcor = "equi",
                      signs = "pos", decay_ps = nu, decay_out = nu,
                      outcome_set = "track_s"))
})

test_that("treated share other than one half: same intercept as the record", {
  expect_same_draws(
    seed = 16,
    new_args = list(n = 400, p = 20, c = 0.5, treated_share = 0.2),
    old_fun  = record$dgp2,
    old_args = list(n = 400, p = 20, outcome = c("linear", "quad1", "exp"),
                    covcor = "iid", overlap = 0.5, signs = "pos", treat_prop = 0.2))
})

# part B. the properties the thesis states

bare <- function(m) {
  m <- unname(m)
  attr(m, "scaled:center") <- NULL
  attr(m, "scaled:scale")  <- NULL
  m
}

big <- 2e5   # draws per check

test_that("the true logit has sd 1/c whatever s, nu and the covariance are", {
  for (cov in c("independent", "ar1", "equicorrelated")) for (s in c(4, 64))
    for (nu in c(1, 0.5)) for (c in c(0.25, 1, 2)) {
      set.seed(21)
      d <- generate_sample(big, p = 64, s = s, nu = nu, covariance = cov, c = c,
                           outcomes = "linear", loadings = "all_active")
      expect_equal(sd(d$logit), 1 / c, tolerance = 0.01)
      expect_equal(mean(d$logit), 0, tolerance = 0.01 / c)
    }
})

test_that("the treated share is 0.5, and other targets are hit", {
  for (c in c(0.05, 0.5, 2)) {
    set.seed(22)
    expect_equal(mean(generate_sample(big, p = 4, c = c, outcomes = "linear")$e),
                 0.5, tolerance = 0.003)
  }
  for (share in c(0.05, 0.2)) {
    set.seed(23)
    d <- generate_sample(big, p = 4, c = 0.5, treated_share = share, outcomes = "linear")
    expect_equal(mean(d$e), share, tolerance = 0.003)
    expect_equal(mean(d$W), share, tolerance = 0.005)
  }
})

test_that("Appendix Table 6: shares of units outside [0.1, 0.9], [0.05, 0.95] and [0.01, 0.99]", {
  table6 <- data.frame(c = c(2, 1, 0.75, 0.5, 0.25, 0.15, 0.05, 0.01),
                       out10 = c(NA, 2.8, 9.9, 27.2, 58.3, 74.2, 91.3, 98.2),
                       out05 = c(NA, 0.3, 2.7, 14.1, 46.2, 65.9, 88.3, 97.7),
                       out01 = c(NA, NA, 0.1, 2.2, 25.1, 49.1, 81.8, 96.3))
  for (i in seq_len(nrow(table6))) {
    c <- table6$c[i]
    set.seed(24)
    e <- generate_sample(big, p = 4, c = c, outcomes = "linear")$e
    emp10 <- 100 * mean(e < 0.10 | e > 0.90)
    emp05 <- 100 * mean(e < 0.05 | e > 0.95)
    expect_equal(100 * share_outside(c, 0.10), emp10, tolerance = 0.015)
    expect_equal(100 * share_outside(c, 0.05), emp05, tolerance = 0.015)
    if (is.na(table6$out10[i])) {
      expect_lt(emp10, 0.1); expect_lt(emp05, 0.1)   # printed as "< 0.1%"
    } else {
      expect_equal(round(100 * share_outside(c, 0.10), 1), table6$out10[i])
      expect_equal(round(100 * share_outside(c, 0.05), 1), table6$out05[i])
    }
    if (is.na(table6$out01[i])) {
      expect_lt(100 * share_outside(c, 0.01), 0.05)                  # printed as "< 0.1%"
    } else {
      expect_equal(round(100 * share_outside(c, 0.01), 1), table6$out01[i])
    }
  }
})

test_that("covariance structures are what the text says", {
  for (cov in c("independent", "ar1", "equicorrelated")) {
    set.seed(25)
    X <- generate_sample(big, p = 6, covariance = cov, outcomes = "linear")$X_true
    target <- switch(cov,
      independent    = diag(6),
      ar1            = 0.5^abs(outer(1:6, 1:6, "-")),
      equicorrelated = 0.5 * diag(6) + 0.5)
    expect_equal(cov(X), target, tolerance = 0.015)
  }
})

test_that("Appendix Table 20: the common factor's share of the index variance", {
  table20 <- data.frame(s = c(4, 16, 64), nu1 = c(0.75, 0.88, 0.93),
                        nu05 = c(0.79, 0.93, 0.98))
  for (i in 1:3) {
    expect_equal(round(common_factor_share(table20$s[i], 1), 2), table20$nu1[i])
    expect_equal(round(common_factor_share(table20$s[i], 0.5), 2), table20$nu05[i])
  }
  # and the draws. var(g) = .5 (sum w)^2 + .5 sum w^2 under equicorrelation
  set.seed(26)
  d <- generate_sample(big, p = 16, s = 16, nu = 0.5, covariance = "equicorrelated",
                       outcomes = "linear", loadings = "all_active")
  w <- (1:16)^(-0.5)
  g <- drop(d$X_true[, 1:16] %*% w)
  expect_equal(var(g), 0.5 * sum(w)^2 + 0.5 * sum(w^2), tolerance = 0.01)
})

test_that("misspecification changes the observed covariates and nothing else", {
  set.seed(27); correct <- generate_sample(2000, p = 30, c = 0.5); after_correct <- .Random.seed
  set.seed(27); ks      <- generate_sample(2000, p = 30, c = 0.5,
                                           misspecification = "kang_schafer")
  same <- c("Y", "W", "X_true", "e", "e0", "logit", "m")
  expect_identical(correct[same], ks[same])
  expect_identical(.Random.seed, after_correct)   # the folds drawn next are the same
  expect_identical(correct$X[, 5:30], ks$X[, 5:30])
  expect_false(any(correct$X[, 1:4] == ks$X[, 1:4]))

  Z <- correct$X_true
  eq3 <- cbind(exp(Z[, 1] / 2),
               Z[, 2] / (1 + exp(Z[, 1])) + 10,
               (0.04 * Z[, 1] * Z[, 3] + 0.6)^3,
               (Z[, 2] + Z[, 4] + 20)^2)
  expect_equal(bare(ks$X[, 1:4]), bare(scale(eq3)), tolerance = 1e-12)

  set.seed(28); correct <- generate_sample(2000, p = 60, s = 16, covariance = "ar1",
                                           loadings = "all_active"); after_correct <- .Random.seed
  set.seed(28); tan     <- generate_sample(2000, p = 60, s = 16, covariance = "ar1",
                                           loadings = "all_active",
                                           misspecification = "tan")
  expect_identical(correct[same], tan[same])
  expect_identical(.Random.seed, after_correct)
  expect_identical(correct$X[, 17:60], tan$X[, 17:60])
  Z <- correct$X_true[, 1:16]
  expect_equal(bare(tan$X[, 1:16]), bare(scale(Z + pmax(Z + 1, 0)^2)),
               tolerance = 1e-12)
})

test_that("the outcome is eq. (2) times kappa plus unit noise, with no treatment effect", {
  set.seed(29); d1 <- generate_sample(big, p = 4, kappa = 1)
  set.seed(29); d2 <- generate_sample(big, p = 4, kappa = 2)
  Z <- d1$X_true
  h <- cbind(linear      = Z[, 1] + 0.5 * Z[, 2] + 0.5 * Z[, 3] + 0.5 * Z[, 4],
             quadratic   = rowSums(pmax(Z[, 1:4], 0)^2),
             exponential = rowSums(exp(Z[, 1:4] / 2)))
  expect_equal(unname(d1$m), unname(h), tolerance = 1e-12)
  expect_equal(unname(d2$m), unname(2 * h), tolerance = 1e-12)
  eps <- d1$Y - d1$m
  expect_equal(eps, d2$Y - d2$m, tolerance = 1e-12)       # same noise scaled signal
  expect_equal(unname(colMeans(eps)), c(0, 0, 0), tolerance = 0.01)
  expect_equal(unname(apply(eps, 2, sd)), c(1, 1, 1), tolerance = 0.01)
  expect_equal(mean(eps[d1$W == 1, "linear"]) - mean(eps[d1$W == 0, "linear"]), 0,
               tolerance = 0.015)
})

test_that("the Oracle Hajek estimator with the true score recovers tau = 0 (pilot check)", {
  hajek <- function(d, y) with(d, sum(W * y / e) / sum(W / e) -
                                   sum((1 - W) * y / e0) / sum((1 - W) / e0))
  for (c in c(0.25, 0.5, 1, 2)) {
    set.seed(30)
    d <- generate_sample(big, p = 4, c = c)
    for (type in colnames(d$Y)) {
        y <- d$Y[, type]
      arm <- function(w, ya) sum(w^2 * (ya - sum(w * ya) / sum(w))^2) / sum(w)^2
      se <- sqrt(arm(1 / d$e[d$W == 1], y[d$W == 1]) +
                 arm(1 / d$e0[d$W == 0], y[d$W == 0]))
      expect_lt(abs(hajek(d, y)), 3 * se)
    }
  }
})
