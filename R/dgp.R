# data generators. dgp2 is current. dgp3 = dgp2 + Tan 2020 misspec and equi cov
# dgp_wc_overlap = Wong and Chan design side study (not part of the main analysis)

#' dgp2. sparse logistic ps on first s covs (or idx_ps) w/ coefs sign * j^-decay_ps
#' Var(eta) = 1 / overlap^2 (only knob moving overlap). intercept gives E[e] = treat_prop
dgp2 <- function(n, p, s = 4,  # s >= 4 (Zhao 2019). s > 4 w/ fixed4 adds weak instruments
                 outcome     = "linear",  # Tan 2020 outcome models
                 misspec     = FALSE,
                 covcor      = c("iid", "ar1"),
                 overlap     = 1,
                 signs       = c("pos", "neg", "mixed", "ks"),
                 decay_ps = 1, 
                 
                 decay_out = 1,
                 outcome_set = c("fixed4", "track_s"),
                 idx_ps = NULL,  # instruments + outcome-only covs (Shortreed and Ertefaie 2017)
                 idx_out = NULL,
                 treat_prop  = 0.5,
                 tau = 0,
                 strength = 1,  # signal to noise knob
                 mu0 = 0){  # 0 keeps existing draws. 210 = Kang and Schafer level (HT vs Hajek)
  outcome     <- match.arg(outcome, c("linear", "quad1", "exp"), several.ok = TRUE)
  covcor      <- match.arg(covcor)
  signs       <- match.arg(signs)
  outcome_set <- match.arg(outcome_set)
  
  .Sig <- function(idx,covcor)  # for overlap normalisation
    if (covcor == "iid") diag(length(idx)) else 0.5^abs(outer(idx, idx, "-"))
  
  if (is.null(idx_ps))  idx_ps  <- seq_len(s) else s <- length(idx_ps)
  if (is.null(idx_out)) idx_out <- idx_ps[if (outcome_set == "fixed4") 1:4 else seq_len(s)]
  stopifnot(!anyDuplicated(idx_ps), !anyDuplicated(idx_out),
            c(idx_ps, idx_out) >= 1, c(idx_ps, idx_out) <= p,
            outcome_set != "fixed4" || length(idx_out) == 4, !misspec || s == 4)
  
  X <- matrix(rnorm(n * p), n, p)
  if (covcor == "ar1" && p >= 2)
    for (j in 2:p) X[, j] <- 0.5 * X[, j - 1] + sqrt(0.75) * X[, j]
  
  X_true <- X
  
  # propensity j^-decay normalised
  # signs ks similar but not the same at decay = 1 (1 -.5 .333 .25)
  beta <- switch(signs,
                 pos   = rep(1, s),
                 neg   = rep(-1, s),
                 mixed = rep(c(1, -1), length.out = s),
                 ks    = rep_len(c(1, -1, 1, 1), s)
  ) / seq_len(s)^decay_ps
  
  eta <- drop( (X[, idx_ps, drop = FALSE] %*% beta) /
                 (sqrt(drop(crossprod(beta, .Sig(idx_ps, covcor) %*% beta))) * overlap) )
  
  alpha <- if (treat_prop == 0.5) 0 else       # treatment allocator shifter
    uniroot(function(a)
      integrate(function(z) plogis(a - z / overlap) * dnorm(z),
                -Inf, Inf)$value - treat_prop,
      c(-50, 50))$root       # wider interval allows more extreme overlap / treat_prop
  
  e  <- plogis(alpha - eta)     # P(W=1)
  e0 <- plogis(eta - alpha)     # P(W=0) kept so the oracle doesn't fail at extreme low overlap
  
  W <- rbinom(n, 1, e)
  
  # outcome signal
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
  
  eps <- rnorm(n)  # one eps shared by all outcomes. no effect on estimates and saves compute
  Y <- vapply(outcome, function(o) mu0 + tau * W + strength * f(o) + eps, numeric(n))
  m <- vapply(outcome, function(o) mu0 + strength * f(o), numeric(n))   # true E[Y(0) | X]
  
  
  if (misspec) {  # Kang-Schafer transforms of the first 4 (s = 4). Y already built from X_true
    j <- idx_ps[1:4]
    X[, j[1]] <- exp(0.5 * X_true[, j[1]])
    X[, j[2]] <- 10 + X_true[, j[2]] / (1 + exp(X_true[, j[1]]))
    X[, j[3]] <- (0.04 * X_true[, j[1]] * X_true[, j[3]] + 0.6)^3
    X[, j[4]] <- (X_true[, j[2]] + X_true[, j[4]] + 20)^2
    X[, j]    <- scale(X[, j])
  }
  
  list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta, e0 = e0,
       m = m)
}

#' Wong and Chan 2018 design as run by Wang and Zubizarreta 2020 (App D.4). side study
#' overlap multiplies the treatment logit and sigma = noise sd (our extension not in either paper)
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
  if (misspec) {  # K&S style transforms. the paper's design
    X[, 1] <- exp(Z[, 1] / 2)
    X[, 2] <- Z[, 2] / (1 + exp(Z[, 1]))
    X[, 3] <- (Z[, 1] * Z[, 3] / 25 + 0.6)^3
    X[, 4] <- (Z[, 2] + Z[, 4] + 20)^2
  }
  
  tau_i <- 1.5 * g  # model A effects (B has none). PATE 0 for both. true ATT about -17.7
  list(Y = Y, W = W, X = X, X_true = Z, tau = 0, tau_i = tau_i,
       e = e, eta = eta, e0 = 1 - e)
}

# dgp3 and helpers

#' Tan 2020 AoS s4.2 transform x_dag = x + {(x + 1)_+}^2. monotone in x
.xdag <- function(x) x + pmax(x + 1, 0)^2
.xdag_mu <- 2 * pnorm(1) + dnorm(1)                         # E[x_dag] std normal X = 1.924660
.xdag_sd <- 3.390312                                        # sd[x_dag] from trunc normal moments

.xdag_corr_cache <- new.env()
#' corr of x_dag(Z1) and x_dag(Z2) for std normals w/ corr r. cached per r
#' replaces the Sigma entries in the tan_ps normalisation (ar1 r = .5^|j-k|. equi r = .5 -> .4675)
.xdag_corr <- function(r) {
  k <- format(r, digits = 17)
  if (!is.null(.xdag_corr_cache[[k]])) return(.xdag_corr_cache[[k]])
  v <- .xdag_corr_raw(r); assign(k, v, envir = .xdag_corr_cache); v
}

#' .xdag_corr() w/o the cache. exact inner expectation then 1-D outer integral
.xdag_corr_raw <- function(r) {
  if (r == 0) return(0)
  if (r == 1) return(1)
  s2 <- sqrt(1 - r^2)
  inner <- function(z) { m <- r * z + 1                     # E[x_dag(Z2) | Z1 = z]
    r * z + (m^2 + s2^2) * pnorm(m / s2) + m * s2 * dnorm(m / s2) }
  e12 <- integrate(function(z) .xdag(z) * inner(z) * dnorm(z), -Inf, Inf,
                   rel.tol = 1e-10)$value
  (e12 - .xdag_mu^2) / .xdag_sd^2
}

.alpha_cache <- new.env()
#' intercept so E[plogis(alpha - eta)] = treat_prop for skewed eta. solved once per cell + cached
#' ref draw has own seed and RNG kind so same alpha in every session. caller RNG restored
.alpha_skew <- function(key, make_eta, treat_prop, nref = 1e5) {
  if (!is.null(.alpha_cache[[key]])) return(.alpha_cache[[key]])
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(20261004, kind = "Mersenne-Twister", normal.kind = "Inversion")
  eta_ref <- make_eta(nref)
  a <- uniroot(function(a) mean(plogis(a - eta_ref)) - treat_prop, c(-50, 50))$root
  assign(key, a, envir = .alpha_cache)
  a
}

#' dgp3 = dgp2 + Tan 2020 x_dag misspec on all of idx_ps (any s) + equi covcor. same draws as dgp2
#' misspec_form ks = dgp2. tan_obs observes x_dag. tan_ps truth linear in x_dag and misspec sees X
dgp3 <- function(n, p, s = 4,
                 outcome     = "linear",
                 misspec     = FALSE,
                 misspec_form = c("ks", "tan_obs", "tan_ps"),
                 covcor      = c("iid", "ar1", "equi"),
                 overlap     = 1,
                 signs       = c("pos", "neg", "mixed", "ks"),
                 decay_ps = 1,
                 decay_out = 1,
                 outcome_set = c("fixed4", "track_s"),
                 idx_ps = NULL,
                 idx_out = NULL,
                 treat_prop  = 0.5,
                 tau = 0,
                 strength = 1,
                 mu0 = 0){
  outcome      <- match.arg(outcome, c("linear", "quad1", "exp"), several.ok = TRUE)
  misspec_form <- match.arg(misspec_form)
  covcor       <- match.arg(covcor)
  signs        <- match.arg(signs)
  outcome_set  <- match.arg(outcome_set)
  .Sig <- function(idx, covcor) switch(covcor,
    iid  = diag(length(idx)),
    ar1  = 0.5^abs(outer(idx, idx, "-")),
    equi = 0.5 * diag(length(idx)) + 0.5)                                    # Sigma_jk = .5 j != k
  if (is.null(idx_ps))  idx_ps  <- seq_len(s) else s <- length(idx_ps)
  if (is.null(idx_out)) idx_out <- idx_ps[if (outcome_set == "fixed4") 1:4 else seq_len(s)]
  stopifnot(!anyDuplicated(idx_ps), !anyDuplicated(idx_out),
            c(idx_ps, idx_out) >= 1, c(idx_ps, idx_out) <= p,
            outcome_set != "fixed4" || length(idx_out) == 4,
            !misspec || misspec_form != "ks" || s == 4)                      # s = 4 only for ks
  X <- matrix(rnorm(n * p), n, p)
  if (covcor == "ar1" && p >= 2)
    for (j in 2:p) X[, j] <- 0.5 * X[, j - 1] + sqrt(0.75) * X[, j]
  if (covcor == "equi")  # one common factor drawn after Z
    X <- sqrt(0.5) * rnorm(n) + sqrt(0.5) * X  # so iid and ar1 draws unchanged
  X_true <- X
  # propensity j^-decay normalised
  beta <- switch(signs,
                 pos   = rep(1, s),
                 neg   = rep(-1, s),
                 mixed = rep(c(1, -1), length.out = s),
                 ks    = rep_len(c(1, -1, 1, 1), s)
  ) / seq_len(s)^decay_ps
  if (misspec_form != "tan_ps") {                                            # dgp2 index
    eta <- drop( (X[, idx_ps, drop = FALSE] %*% beta) /
                   (sqrt(drop(crossprod(beta, .Sig(idx_ps, covcor) %*% beta))) * overlap) )
    alpha <- if (treat_prop == 0.5) 0 else
      uniroot(function(a)
        integrate(function(z) plogis(a - z / overlap) * dnorm(z),
                  -Inf, Inf)$value - treat_prop,
        c(-50, 50))$root
  } else {  # tan_ps. index linear in x_dag
    R <- .Sig(idx_ps, covcor)
    if (covcor != "iid") {  # use corr of x_dag for ar1 and equi
      u <- unique(as.vector(R)); R[] <- vapply(u, .xdag_corr, numeric(1))[match(R, u)]
    }
    sd_idx <- sqrt(drop(crossprod(beta, R %*% beta)))
    make_eta <- function(Xm) drop(((.xdag(Xm[, idx_ps, drop = FALSE]) - .xdag_mu) / .xdag_sd) %*% beta) /
      (sd_idx * overlap)
    eta <- make_eta(X)
    pr  <- max(idx_ps)  # ref draw only needs cols 1..max(idx_ps)
    key <- paste(covcor, overlap, signs, decay_ps, paste(idx_ps, collapse = "."), treat_prop)
    alpha <- .alpha_skew(key, function(m) {  # eta is skewed so E[e] = treat_prop needs a solve
      Xr <- matrix(rnorm(m * pr), m, pr)
      if (covcor == "ar1" && pr >= 2) for (k in 2:pr) Xr[, k] <- 0.5 * Xr[, k - 1] + sqrt(0.75) * Xr[, k]
      if (covcor == "equi") Xr <- sqrt(0.5) * rnorm(m) + sqrt(0.5) * Xr
      make_eta(Xr) }, treat_prop)
  }
  e  <- plogis(alpha - eta)     # P(W=1)
  e0 <- plogis(eta - alpha)     # P(W=0)
  W <- rbinom(n, 1, e)
  # outcome signal same as dgp2. always from X_true
  Xo <- X[, idx_out, drop = FALSE]
  if (outcome_set == "fixed4") {
    a <- c(1, 0.5, 0.5, 0.5)
    f <- function(o) switch(o, linear = drop(Xo %*% a),
                            quad1 = rowSums(pmax(Xo, 0)^2), exp = rowSums(exp(Xo / 2)))
  } else {  # track_s. normalised w/ .Sig so equi is covered
    a <- 1 / seq_along(idx_out)^decay_out
    a <- a / sqrt(drop(crossprod(a, .Sig(idx_out, covcor) %*% a)))
    f <- function(o) switch(o, linear = drop(Xo %*% a),
                            quad1 = drop(pmax(Xo, 0)^2 %*% a) - 0.5 * sum(a),
                            exp   = drop(exp(Xo / 2) %*% a) - exp(1/8) * sum(a))
  }
  eps <- rnorm(n)
  Y <- vapply(outcome, function(o) mu0 + tau * W + strength * f(o) + eps, numeric(n))
  m <- vapply(outcome, function(o) mu0 + strength * f(o), numeric(n))   # true E[Y(0) | X]
  # what the analyst observes. done after W e Y so correct and misspec draws match by rep
  j <- idx_ps
  if (misspec && misspec_form == "ks") {
    j <- idx_ps[1:4]
    X[, j[1]] <- exp(0.5 * X_true[, j[1]])
    X[, j[2]] <- 10 + X_true[, j[2]] / (1 + exp(X_true[, j[1]]))
    X[, j[3]] <- (0.04 * X_true[, j[1]] * X_true[, j[3]] + 0.6)^3
    X[, j[4]] <- (X_true[, j[2]] + X_true[, j[4]] + 20)^2
    X[, j]    <- scale(X[, j])
  }
  if ((misspec && misspec_form == "tan_obs") || (!misspec && misspec_form == "tan_ps"))
    X[, j] <- scale(.xdag(X_true[, j, drop = FALSE]))  # scaled per sample as K&S
  list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta, e0 = e0,
       m = m)
}
