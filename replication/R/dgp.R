# the thesis design (Methods 3.2). one call draws one sample. no treatment effect
# draws in the order of dgp2() and dgp3() in the root R/dgp.R so stored draws match (test_dgp.R)

#' one sample. Y W X (observed) X_true e e0 logit m = E[Y(0) | X] and tau = 0
#' logit = standardised index w/ weights j^(-nu) on the s active covs. e = plogis(alpha - logit)
generate_sample <- function(n, p, s = 4, nu = 1,
                            covariance = c("independent", "ar1", "equicorrelated"),
                            c = 1, kappa = 1,
                            outcomes = c("linear", "quadratic", "exponential"),
                            misspecification = c("none", "kang_schafer", "tan"),
                            loadings = c("first_four", "all_active"),
                            treated_share = 0.5) {
  covariance       <- match.arg(covariance)
  misspecification <- match.arg(misspecification)
  loadings         <- match.arg(loadings)
  outcomes         <- match.arg(outcomes, several.ok = TRUE)
  stopifnot(s >= 4, s <= p, c > 0, misspecification != "kang_schafer" || s == 4)

  # covariates. AR(1) w/ rho = 0.5 or one common factor for equicorrelated
  X <- matrix(rnorm(n * p), n, p)
  if (covariance == "ar1" && p >= 2)
    for (j in 2:p) X[, j] <- 0.5 * X[, j - 1] + sqrt(0.75) * X[, j]
  if (covariance == "equicorrelated")
    X <- sqrt(0.5) * rnorm(n) + sqrt(0.5) * X
  X_true <- X
  Sigma <- function(idx) switch(covariance,
    independent    = diag(length(idx)),
    ar1            = 0.5^abs(outer(idx, idx, "-")),
    equicorrelated = 0.5 * diag(length(idx)) + 0.5)
  active <- seq_len(s)

  # true score eq. (1). index is standardised so only c moves overlap
  w     <- 1 / active^nu
  logit <- drop(X[, active, drop = FALSE] %*% w) / (sqrt(drop(crossprod(w, Sigma(active) %*% w))) * c)
  alpha <- intercept_for_share(treated_share, c)
  e  <- plogis(alpha - logit)
  e0 <- plogis(logit - alpha)          # 1 - e without rounding to zero
  W  <- rbinom(n, 1, e)

  # outcome eq. (2). kappa h(X) + unit normal noise. h = 3 functions of the confounders (Tan 2020b)
  if (loadings == "first_four") {
    Xo <- X[, 1:4, drop = FALSE]
    a  <- c(1, 0.5, 0.5, 0.5)
    h  <- function(type) switch(type,
      linear      = drop(Xo %*% a),
      quadratic   = rowSums(pmax(Xo, 0)^2),
      exponential = rowSums(exp(Xo / 2)))
  } else {                             # all s active covs. unit var index. outcomes centred
    Xo <- X[, active, drop = FALSE]
    a  <- 1 / active^nu
    a  <- a / sqrt(drop(crossprod(a, Sigma(active) %*% a)))
    h  <- function(type) switch(type,
      linear      = drop(Xo %*% a),
      quadratic   = drop(pmax(Xo, 0)^2 %*% a) - 0.5 * sum(a),
      exponential = drop(exp(Xo / 2) %*% a) - exp(1 / 8) * sum(a))
  }
  m <- vapply(outcomes, function(type) kappa * h(type), numeric(n))
  Y <- m + rnorm(n)

  # what the analyst sees. standardised in sample (Tan 2020b). all draws are done by here
  # so both specs share the draw and leave the RNG at the same point
  if (misspecification == "kang_schafer") {   # Kang and Schafer 2007. needs s = 4
    Z <- X_true
    X[, 1] <- exp(Z[, 1] / 2)
    X[, 2] <- Z[, 2] / (1 + exp(Z[, 1])) + 10
    X[, 3] <- (0.04 * Z[, 1] * Z[, 3] + 0.6)^3
    X[, 4] <- (Z[, 2] + Z[, 4] + 20)^2
    X[, 1:4] <- scale(X[, 1:4])
  }
  if (misspecification == "tan")   # Tan 2020a. every active cov
    X[, active] <- scale(tan_transform(X_true[, active, drop = FALSE]))

  list(Y = Y, W = W, X = X, X_true = X_true, e = e, e0 = e0, logit = logit, m = m, tau = 0)
}

#' Tan's transform. monotone in x (Tan 2020a section 4.2)
tan_transform <- function(x) x + pmax(x + 1, 0)^2

#' intercept giving E[e(X)] = treated_share. zero at one half by symmetry
intercept_for_share <- function(treated_share, c) {
  if (treated_share == 0.5) return(0)
  mean_score <- function(alpha) integrate(function(z) plogis(alpha - z / c) * dnorm(z), -Inf, Inf)$value
  uniroot(function(alpha) mean_score(alpha) - treated_share, c(-50, 50))$root
}

#' share of units w/ a true score outside lower to 1 - lower (Appendix Table 6)
share_outside <- function(c, lower) 2 * pnorm(-c * qlogis(1 - lower))

#' share of the index variance from the common factor under equicorrelation (App Table 20)
common_factor_share <- function(s, nu) {
  w <- seq_len(s)^(-nu)
  sum(w)^2 / (sum(w)^2 + sum(w^2))
}
