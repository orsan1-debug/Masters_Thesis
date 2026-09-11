# Some tuning experiments using the balnet experimental branch (has "cv.boot.balnet")
# Estimating E[Y(1)]
# devtools::install_github("erikcs/balnet", ref = "cv-experimental", subdir = "r-package/balnet")
rm(list = ls())
set.seed(42)

library(balnet)


#' High-dimensional DGP for estimating a counterfactual treated mean.
#'
#' Features are AR(1)-correlated with parameter rho. The first `s_prop`
#' features drive the propensity; the first `s_y` features drive the
#' outcome. They overlap on the first min(s_prop, s_y) coordinates (the
#' "confounders"); any remaining propensity-only coordinates are
#' instruments and outcome-only coordinates are pure predictors.
#'
#' Target: mu1 = E[Y(1)]. By construction E[X] = 0, so mu1 = 0.
#' Only Y for treated units is used by mu1 estimators; Y for controls
#' is returned for convenience but plays no role in the target.
#'
#' Overlap is controlled via the propensity coefficient magnitude:
#'   Var(eta) grows quadratically in the scale; large Var(eta) drives
#'   propensities toward 0 and 1.
gen_data <- function(n = 1000, p = 100, rho = 0.5,
                     s_prop = 5, s_y = 5,
                     overlap = c("bad", "moderate", "good"),
                     sigma_y = 1,
                     seed = NULL) {
  overlap <- match.arg(overlap)
  if (!is.null(seed)) set.seed(seed)

  # AR(1) covariance across features
  Sigma <- rho ^ abs(outer(1:p, 1:p, `-`))
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)

  # Propensity coefficient scale controls overlap
  c_prop <- switch(overlap,
                   good     = 0.7,   # SD(eta) ~ 1,  most props in [0.2, 0.8]
                   moderate = 1.5,   # SD(eta) ~ 2,  some near 0/1
                   bad      = 2.5)   # SD(eta) ~ 3.7, many near 0/1
                   # bad      = 4)   # SD(eta) ~ 3.7, many near 0/1

  beta_prop <- rep(0, p)
  beta_prop[seq_len(s_prop)] <- c_prop / sqrt(s_prop)

  beta_y <- rep(0, p)
  beta_y[seq_len(s_y)] <- 1 / sqrt(s_y)

  eta  <- as.numeric(X %*% beta_prop)
  prop <- plogis(eta)
  W    <- rbinom(n, 1, prop)

  # Y is Y(1): potential outcome under treatment. mu1 = E[Y(1)] = 0.
  Y <- as.numeric(X %*% beta_y) + rnorm(n, sd = sigma_y)

  list(X = X, W = W, Y = Y, prop = prop, true = 0,
       active_prop = which(beta_prop != 0),
       active_y    = which(beta_y != 0))
}



n = 1000
true = 0 # always zero

## sparse lasso with lowest imbalance = best
alpha = 1
sy = 5
p = 100
sigma = 1

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- parallel::mclapply(seq_len(10), function(rep.i) {
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W

  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)

  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)

  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "imbalance.mean",
                       maxit = 1e4, tol = 1e-5)

  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "imbalance.inf",
                       maxit = 1e4, tol = 1e-5)

  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)

  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)

  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,

    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,

    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,

    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,

    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,

    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
# }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)



## dense lasso with interior min
alpha = 1
sy = 100
p = 100
sigma = 1

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- parallel::mclapply(seq_len(25), function(rep.i) {
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W

  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)

  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)

  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.mean",
                     maxit = 1e4, tol = 1e-5)

  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.inf",
                     maxit = 1e4, tol = 1e-5)

  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)

  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)

  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,

    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,

    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,

    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,

    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,

    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
  # }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)



## RMSE path "reversal" depending on outcome noise
## dense lasso with interior min
alpha = 1
sy = 100
p = 100
sigma = 5

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- parallel::mclapply(seq_len(25), function(rep.i) {
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W

  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)

  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)

  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.mean",
                     maxit = 1e4, tol = 1e-5)

  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.inf",
                     maxit = 1e4, tol = 1e-5)

  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)

  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)

  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,

    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,

    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,

    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,

    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,

    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
  # }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)






