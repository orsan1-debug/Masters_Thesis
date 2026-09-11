
# --- DGP 2 ---

# Parameters:
#   X:  observables
#   X_true: untransformed covariates (== X unless misspec = TRUE)
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
#   correlation:      [iid/ar1] from DGPs Erik sent,  https://arxiv.org/pdf/1706.03461.pdf
#   s:                sparsity parameter, number of active covariates in propensity score (and outcome, when outcome_set = "track_s") [Zhao 2019]
#   signs:            sign pattern for propensity coefficients [pos/neg/mixed/ks]


#   decay_ps:          coeff decay exponent for propensity score; 1 = default. Higher = concentrated signal, lower = denser signal
#   decay_out:         coeff decay exponent for outcome; 1 = default, benchmark. 

#   treat_prop:       target treated proportion of population P(W=1)

#   outcome_set:      [fixed4/track_s] fixed4 = Tan outcomes on X_1..X_4 (default, reproduces DGP1).
#                     track_s = outcome loadings j^(-decay) on X_1..X_s, normalised + centered, all s covariates confounders

#   idx_ps/idx_out:   which covariates affect treatment / outcome; allows instruments and
#                     prognostic-only covariates (Shortreed & Ertefaie 2017)  #### !!!!! - TEST - !!!!! IMPLEMENTATION !!!!!

mor

# p/n ratio tested in grid.
# weak instruments added by increasing s>4, and setting outcome to fixed4.

#   eta sd-normalized to Var(eta) = 1/overlap^2.Isolates overlap changes to overlap knob.


# -no removed strength (SNR knob, multiplier on outcome function), maybe relevent for AIPW.




# --- DGP2 (WIP) ---
.Sig <- function(idx, covcor)
  if (covcor == "iid") diag(length(idx)) else 0.5^abs(outer(idx, idx, "-"))            #overlap normalization

generate_data <- function(n, p, s = 4,
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
  outcome     <- match.arg(outcome, c("linear", "quad1", "exp"), several.ok = TRUE)
  covcor      <- match.arg(covcor)
  signs       <- match.arg(signs)
  outcome_set <- match.arg(outcome_set)                                        #catch bad inputs
  
  
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
  
  e <- 1 / (1 + exp(eta - alpha))
  
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
  Y <- drop(tau * W + sapply(outcome, f) + rnorm(n))
  
  if (misspec) {
    j <- idx_ps[1:4]
    X[, j[1]] <- exp(0.5 * X_true[, j[1]])
    X[, j[2]] <- 10 + X_true[, j[2]] / (1 + exp(X_true[, j[1]]))
    X[, j[3]] <- (0.04 * X_true[, j[1]] * X_true[, j[3]] + 0.6)^3
    X[, j[4]] <- (X_true[, j[2]] + X_true[, j[4]] + 20)^2
    X[, j]    <- scale(X[, j])
  }
  
  list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta)
}




