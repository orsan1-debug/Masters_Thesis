
### NOTES ####

#DGP2 IS CURRENT, DGP 1 LEGACY.

# --- DGP 1 --- # legacy 

# Parameters:
#   X:  observables
#   e:  true  propensity score
#   Y:  outcome 
#   W:  treatment
#   tau: ATE       - normalised hajek estimators included for oracle and glmnet variables

#features :
#   n:                sample size
#   p:                number of covariates
#   overlap:                overlap scaler
#   mispeccified:     [True/False] http://arxiv.org/abs/1710.08074v1
#   outcome type:     linear/quadratic/exponential http://arxiv.org/abs/1710.08074v1
#   correlation:      [iid/ar1],  https://arxiv.org/pdf/1706.03461.pdf

dgp1 <- function(n, p, outcome = "linear", misspec = FALSE,
                 covcor = "iid", overlap = 1) {
  
  if (covcor == "iid") {
    X <- matrix(rnorm(n * p), n, p)
  } else if (covcor == "ar1") {
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = stats::toeplitz(0.5^seq(0, p - 1)))
  }
  
  e = 1 / (1 + exp((X[, 1] - 0.5*X[, 2] + 0.25*X[,3] + 0.1*X[,4]) / overlap))
  W <- rbinom(n, 1, e)
  tau <- 0
  
  eps <- rnorm(n)
  f <- list(linear = X[,1] + 0.5*X[,2] + 0.5*X[,3] + 0.5*X[,4],
            quad1  = rowSums(pmax(X[, 1:4], 0)^2),
            exp    = rowSums(exp(X[, 1:4] / 2)))
  Y <- vapply(f[outcome], function(m) m + tau * W + eps, numeric(n))
  
  if (misspec) {
    Xreal <- X
    X[,1] <- exp(0.5 * Xreal[,1])
    X[,2] <- 10 + Xreal[,2] / (1 + exp(Xreal[,1]))
    X[,3] <- (0.04 * Xreal[,1] * Xreal[,3] + 0.6)^3
    X[,4] <- (Xreal[,2] + Xreal[,4] + 20)^2
    X[,1:4] <- scale(X[,1:4])
  }
  
  list(Y = Y, W = W, X = X, tau = tau, e = e)
}


#{ --- DGP 2 --- 
#additional features 
#   s:                 =>4 . sparsity parameter, number of active covariates in propensity score (and outcome, when outcome_set = "track_s") [Zhao 2019]
#   signs:             sign pattern for propensity coefficients [pos/neg/mixed/ks]
#   decay_ps:          coeff decay exponent for propensity score; 1 = default. Higher = concentrated signal, lower = denser signal
#   decay_out:         coeff decay exponent for outcome; 1 = default, benchmark. 
#   treat_prop:       target treated proportion of population P(W=1)
#   outcome_set:      [fixed4/track_s] fixed4 = Tan outcomes on X_1..X_4 (default, reproduces DGP1).
#                     track_s = outcome loadings j^(-decay) on X_1..X_s, normalised + centered, all s covariates confounders
#                     strength = 1, signal to noise ration, increasing increases SNR (signal to noise ration on outcome). 

#   idx_ps/idx_out:   which covariates affect treatment / outcome; allows instruments and
#                     propensity-only covariates (Shortreed & Ertefaie 2017)  #### !!!!! - TEST - !!!!! IMPLEMENTATION !!!!!

# p/n ratio tested in grid.
# weak instruments added by increasing s>4, and setting outcome to fixed4.

#   eta sd-normalized to Var(eta) = 1/overlap^2.Isolates overlap changes to overlap knob.

dgp2 <- function(n, p, s = 4,
                          outcome     = "linear",
                          misspec     = FALSE,
                          covcor      = c("iid", "ar1"),
                          overlap     = 1,
                          signs       = c("pos", "neg", "mixed", "ks"),  #ks similar but not the same at decay = 1, (1,-.5, 0.333, 0.25)
                          decay_ps = 1, 
                 
                          decay_out = 1,
                          outcome_set = c("fixed4", "track_s"),
                          idx_ps = NULL,
                          idx_out = NULL,
                          treat_prop  = 0.5,
                          tau = 0,
                          strength = 1){
  outcome     <- match.arg(outcome, c("linear", "quad1", "exp"), several.ok = TRUE)   #shared error draws: no effect on estimates, reduces compute time
  covcor      <- match.arg(covcor)
  signs       <- match.arg(signs)
  outcome_set <- match.arg(outcome_set)                                        #catch bad inputs
  
  .Sig <- function(idx,covcor)
    if (covcor == "iid") diag(length(idx)) else 0.5^abs(outer(idx, idx, "-"))     #overlap normalization
  
  if (is.null(idx_ps))  idx_ps  <- seq_len(s) else s <- length(idx_ps)
  if (is.null(idx_out)) idx_out <- idx_ps[if (outcome_set == "fixed4") 1:4 else seq_len(s)]
  stopifnot(!anyDuplicated(idx_ps), !anyDuplicated(idx_out),
            c(idx_ps, idx_out) >= 1, c(idx_ps, idx_out) <= p,
            outcome_set != "fixed4" || length(idx_out) == 4, !misspec || s == 4)       # propensity-active columns
  
  X <- matrix(rnorm(n * p), n, p)
  if (covcor == "ar1" && p >= 2)
    for (j in 2:p) X[, j] <- 0.5 * X[, j - 1] + sqrt(0.75) * X[, j]
  
  X_true <- X
  
  # --- propensity: j^-decay, normalized ---
  beta <- switch(signs,
                 pos   = rep(1, s),
                 neg   = rep(-1, s),
                 mixed = rep(c(1, -1), length.out = s),
                 ks    = rep_len(c(1, -1, 1, 1), s)
  ) / seq_len(s)^decay_ps
  
  eta <- drop( (X[, idx_ps, drop = FALSE] %*% beta) /
                 (sqrt(drop(crossprod(beta, .Sig(idx_ps, covcor) %*% beta))) * overlap) )   #normalisation 
  
  alpha <- if (treat_prop == 0.5) 0 else       # treatment allocator shifter
    uniroot(function(a)
      integrate(function(z) plogis(a - z / overlap) * dnorm(z),
                -Inf, Inf)$value - treat_prop,
      c(-50, 50))$root       # c(-x,x) = increased x allows for more extreme overlap/treatment allocation conditions
  
  e  <- plogis(alpha - eta)     # P(W=1)
  e0 <- plogis(eta - alpha)     # P(W=0),   for extreme low overlap, stops oracle failing.    
  
  W <- rbinom(n, 1, e)
  
  # --- outcome signal ---
  Xo <- X[, idx_out, drop = FALSE]
  if (outcome_set == "fixed4") {
    a <- c(1, 0.5, 0.5, 0.5)
    f <- function(o) switch(o, linear = drop(Xo %*% a),
                            quad1 = rowSums(pmax(Xo, 0)^2), exp = rowSums(exp(Xo / 2)))
  } else {  # track_s
    a <- 1 / seq_along(idx_out)^decay_out
    a <- a / sqrt(drop(crossprod(a, .Sig(idx_out, covcor) %*% a)))
    f <- function(o) switch(o, linear = drop(Xo %*% a),
                            quad1 = drop(pmax(Xo, 0)^2 %*% a) - 0.5 * sum(a),
                            exp   = drop(exp(Xo / 2) %*% a) - exp(1/8) * sum(a))
  }
  
  eps <- rnorm(n)
  Y <- vapply(outcome, function(o) tau * W + strength * f(o) + eps, numeric(n))
  
  if (misspec) {
    j <- idx_ps[1:4]
    X[, j[1]] <- exp(0.5 * X_true[, j[1]])
    X[, j[2]] <- 10 + X_true[, j[2]] / (1 + exp(X_true[, j[1]]))
    X[, j[3]] <- (0.04 * X_true[, j[1]] * X_true[, j[3]] + 0.6)^3
    X[, j[4]] <- (X_true[, j[2]] + X_true[, j[4]] + 20)^2
    X[, j]    <- scale(X[, j])
  }
  
  list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta, e0 = e0)
}



## --- merged from R/dgp_wc.R (Phase 3, verbatim) ----------------------------

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

# Overlap / noise variant (our extension, not in either paper) ------------
# Identical to dgp_wc() except: the treatment logit is multiplied by
# overlap (1 reproduces dgp_wc(); larger pushes propensities toward 0 and
# 1, the device Wang & Zubizarreta use in their RHC study), and the outcome
# noise SD is sigma (1 reproduces dgp_wc()). Z, W and eps / sigma are
# identical to dgp_wc() for the same seed when overlap = 1.
dgp_wc_overlap <- function(n, overlap = 1, sigma = 1,
                           outcome = c("A", "B"), misspec = TRUE) {
  outcome <- match.arg(outcome, several.ok = TRUE)
  
  Z   <- matrix(rnorm(n * 10), n, 10)
  eta <- overlap * (-Z[, 1] - 0.1 * Z[, 4])
  e   <- 1 / (1 + exp(-eta))
  W   <- rbinom(n, 1, e)
  
  g   <- 27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]
  f   <- function(o) switch(o,
                            A = 210 + (1.5 * W - 0.5) * g,
                            B = Z[, 1] * Z[, 2]^3 * Z[, 3]^2 * Z[, 4] + Z[, 4] * abs(Z[, 1])^0.5
  )
  eps <- sigma * rnorm(n)
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

## --- merged from R/dgp_cv.R (Phase 3, verbatim) ----------------------------

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


# *** Parallel helper ***
#' added in for parallel runs on windows
run_par <- function(X, FUN, ...) {
  parallel::clusterSetRNGStream(cl)     
  parallel::clusterExport(cl, setdiff(ls(globalenv()), "cl"))
  parallel::parLapply(cl, X, FUN) # X reps per regieme
}
