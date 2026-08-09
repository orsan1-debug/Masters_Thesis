
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

#   idx_ps/idx_out:   which covariates affect treatment / outcome; allows instruments and
#                     propensity-only covariates (Shortreed & Ertefaie 2017)  #### !!!!! - TEST - !!!!! IMPLEMENTATION !!!!!

# p/n ratio tested in grid.
# weak instruments added by increasing s>4, and setting outcome to fixed4.

#   eta sd-normalized to Var(eta) = 1/overlap^2.Isolates overlap changes to overlap knob.
#note: removed strength (SNR knob, multiplier on outcome function), maybe relevent for AIPW.

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
                          tau = 0) {
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
  Y <- vapply(outcome, function(o) tau * W + f(o) + eps, numeric(n))
  
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

