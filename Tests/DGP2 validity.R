
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
#                     propensity-only covariates (Shortreed & Ertefaie 2017)  #### !!!!! - TEST - !!!!! IMPLEMENTATION !!!!!

# p/n ratio tested in grid.
# weak instruments added by increasing s>4, and setting outcome to fixed4.

#   eta sd-normalized to Var(eta) = 1/overlap^2.Isolates overlap changes to overlap knob.


# -no removed strength (SNR knob, multiplier on outcome function), maybe relevent for AIPW.




# --- DGP2  ---
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

### - Validity Tests - #####


### ---- TREAT_PROP CALIBRATION ---- ###
# alpha is solved so E[e] = treat_prop, whatever the overlap.

grid <- expand.grid(treat_prop = c(0.2, 0.35, 0.5), overlap = c(0.25, 0.5, 1, 2))

for (i in 1:nrow(grid)) {
  set.seed(i)
  d <- generate_data(500000, p = 60, s = 4,
                     overlap    = grid$overlap[i],
                     treat_prop = grid$treat_prop[i])
  grid$mean_e[i] <- mean(d$e)
  grid$mean_W[i] <- mean(d$W)
  grid$se[i] <- sd(d$e) / sqrt(500000)
}
grid$z <- (grid$mean_e - grid$treat_prop) / grid$se
round(grid, 4)
max(abs(grid$mean_e - grid$treat_prop))   
max(abs(grid$mean_W - grid$treat_prop)) 
max(abs(grid$z))   # PASS

### ---- TESTING OVERLAP NORMALISER ---- ###

n <- 100000

grid <- expand.grid(overlap  = c(0.25, 0.5, 1, 2),
                    s        = c(4, 10, 25),
                    covcor   = c("iid", "ar1"),
                    decay_ps = c(0.5, 1, 2),
                    signs    = c("pos", "mixed", "ks"),
                    stringsAsFactors = FALSE)

grid$sd_eta <- NA
grid$target <- 1 / grid$overlap

for (i in 1:nrow(grid)) {
  set.seed(i)
  d <- generate_data(n, p = 60,
                     s        = grid$s[i],
                     overlap  = grid$overlap[i],
                     covcor   = grid$covcor[i],
                     decay_ps = grid$decay_ps[i],
                     signs    = grid$signs[i])
  grid$sd_eta[i] <- sd(d$eta)
}

grid$rel <- (grid$sd_eta - grid$target) / grid$target
grid$z   <- grid$rel * sqrt(2 * n)


mean(grid$rel)              # on target.
max(abs(grid$rel))          # worst cell 0.666%  from target
max(abs(grid$z))            # worst cell 2.98 MCSEs from target; within noise
subset(grid, abs(z) > 4)    # none above. 

### ----  p is invariant from propensity  ---- ###
set.seed(9); a <- generate_data(20000, p = 60,  s = 4)$eta
set.seed(9); b <- generate_data(20000, p = 300, s = 4)$eta
identical(a, b)             # TRUE


### ---- MISSPEC TOUCHES ONLY X ---- ###
set.seed(1); a <- generate_data(5000, p = 20, s = 4)
set.seed(1); b <- generate_data(5000, p = 20, s = 4, misspec = TRUE)

identical(a[c("Y","W","e","eta")], b[c("Y","W","e","eta")])   
all.equal(a$X[, 5:20], b$X[, 5:20])                            # PASS


### ---- OUTCOMES SHARE (X, W) ---- ###
set.seed(1); m <- generate_data(5000, p = 20, s = 4, outcome = c("linear","quad1","exp"))
set.seed(1); l <- generate_data(5000, p = 20, s = 4, outcome = "linear")

identical(m$W, l$W) && identical(m$X, l$X)     
all.equal(m$Y[, "linear"], l$Y[, "linear"])             # PASS

### ----SIGNS WORKS, POS/NEG SYMMETRY, AND POS MAKES cov(X,W(1)) negative ---- ###

set.seed(1); pos   <- generate_data(200000, p = 60, s = 4, signs = "pos")
set.seed(1); neg   <- generate_data(200000, p = 60, s = 4, signs = "neg")
set.seed(1); mixed <- generate_data(200000, p = 60, s = 4, signs = "mixed")

round(data.frame(
  pos   = c(cor(pos$X[,1],   pos$e),   cor(pos$X[,2],   pos$e)),
  neg   = c(cor(neg$X[,1],   neg$e),   cor(neg$X[,2],   neg$e)),
  mixed = c(cor(mixed$X[,1], mixed$e), cor(mixed$X[,2], mixed$e)),
  row.names = c("X1", "X2")
), 3)

max(abs(neg$eta + pos$eta))   # PASS


### ---- eta USES idx_ps COLUMNS, MISSPEC HITS ONLY THOSE ---- ###

idx <- c(5, 12, 30, 41)                     # lags 7, 18, 11
set.seed(1)
d <- generate_data(5000, p = 60, idx_ps = idx, covcor = "ar1",
                   overlap = 0.5, misspec = TRUE)

b <- (1 / seq_along(idx))
b <- b / (sqrt(drop(crossprod(b, 0.5^abs(outer(idx, idx, "-")) %*% b))) * 0.5)

all.equal(d$eta, drop(d$X_true[, idx] %*% b))          # TRUE
identical(which(colSums(d$X != d$X_true) > 0), as.integer(idx))   # PASS


### ---- Y LOADS ONLY ON idx_out, INSTRUMENTS EXCLUDED ---- ###
set.seed(1)
d <- generate_data(200000, p = 60, tau = 1,
                   idx_ps  = c(3, 9, 17, 25, 40, 55),
                   idx_out = c(3, 9, 46, 52))

true <- c(1, numeric(60))                              # W first
true[1 + c(3, 9, 46, 52)] <- c(1, 0.5, 0.5, 0.5)

m <- lm(d$Y ~ d$W + d$X_true)
max(abs((coef(m)[-1] - true) / sqrt(diag(vcov(m)))[-1]))   # PASS
max(abs(cor(d$eta, d$X_true[, c(17, 25, 40, 55)])))        # corr=.28 PASS



### ---- ORACLE RECOVERS tau, NAIVE DOES NOT ---- ###
set.seed(1)
d <- generate_data(200000, p = 60, s = 4, tau = 1)

hajek <- function(w, e, y)
  sum(w*y/e)/sum(w/e) - sum((1-w)*y/(1-e))/sum((1-w)/(1-e))

hajek(d$W, d$e, d$Y)                        # PASS
mean(d$Y[d$W == 1]) - mean(d$Y[d$W == 0])   # PASS (is biased)

# overlap = 0.25 checked separately: 200 seeds at n = 20000,
# mean 0.965 (MCSE 0.027), oracle Hajek SD approx 0.38.


set.seed(1)
d <- generate_data(2e5, p = 60, s = 6, outcome = c("linear","quad1","exp"),
                   outcome_set = "track_s")
a <- (1/1:6) / sqrt(sum((1/1:6)^2))
max(abs(coef(lm(d$Y[,"linear"] ~ d$X_true[,1:6]))[-1] - a))  # loadings recover a
sd(d$Y[,"linear"])                                            # ~ sqrt(2): unit signal variance
colMeans(d$Y)[c("quad1","exp")]                               # ~ 0; |.| < 3*sd(Y)/sqrt(n) ~ 0.01
# repeat one cell with covcor = "ar1", decay_out = 2


set.seed(2)
d <- generate_data(2e5, p = 60, s = 6, outcome = c("linear","quad1","exp"),
                   outcome_set = "track_s", covcor = "ar1", decay_out = 2)
S <- 0.5^abs(outer(1:6, 1:6, "-"))
a <- 1/(1:6)^2; a <- a / sqrt(drop(crossprod(a, S %*% a)))
max(abs(coef(lm(d$Y[,"linear"] ~ d$X_true[,1:6]))[-1] - a))  # < ~0.005
sd(d$Y[,"linear"])                                            # √2 (a'Σa = 1 by construction)
colMeans(d$Y)[c("quad1","exp")]                               # ~0; centering exact since ar1 marginals stay N(0,1)
