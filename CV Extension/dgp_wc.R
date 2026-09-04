# dgp_wc.R -----------------------------------------------------------------
# Side study: Wong & Chan (2018) design as run by Wang & Zubizarreta (2020,
# Biometrika 107:93-105, Appendix D.4). Standalone; not part of the main
# analysis sequence.
# Ten latent N(0,1) covariates Z; propensity and outcomes are functions of
# Z. The analyst observes K&S-style transforms of Z (misspec = TRUE, the
# paper's design) or Z itself. Returns one draw with columns of Y sharing
# (Z, W, eps). tau = 0 is the PATE for both models; tau_i are the model A
# individual effects (model B has none), so SATE = mean(tau_i) and
# SATT = mean(tau_i[W == 1]); the true ATT is about -17.7.

dgp_wc <- function(n, outcome = c("A", "B"), misspec = TRUE) {
  outcome <- match.arg(outcome, several.ok = TRUE)
  
  Z   <- matrix(rnorm(n * 10), n, 10)
  eta <- -Z[, 1] - 0.1 * Z[, 4]
  e   <- 1 / (1 + exp(-eta))
  W   <- rbinom(n, 1, e)
  
  g   <- 27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]
  f   <- function(o) switch(o,
                            A = 210 + (1.5 * W - 0.5) * g,
                            B = Z[, 1] * Z[, 2]^3 * Z[, 3]^2 * Z[, 4] + Z[, 4] * abs(Z[, 1])^0.5
  )
  eps <- rnorm(n)
  Y   <- vapply(outcome, \(o) f(o) + eps, numeric(n))
  
  X <- Z
  if (misspec) {
    X[, 1] <- exp(Z[, 1] / 2)
    X[, 2] <- Z[, 2] / (1 + exp(Z[, 1]))
    X[, 3] <- (Z[, 1] * Z[, 3] / 25 + 0.6)^3
    X[, 4] <- (Z[, 2] + Z[, 4] + 20)^2
  }
  
  tau_i <- 1.5 * g                                   # model A; model B is 0
  list(Y = Y, W = W, X = X, X_true = Z, tau = 0, tau_i = tau_i,
       e = e, eta = eta, e0 = 1 - e)
}