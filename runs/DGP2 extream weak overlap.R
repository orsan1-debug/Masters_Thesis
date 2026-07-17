rm(list = ls())
library(balnet); library(glmnet); library(parallel)

RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
num.sim <- 1000
seeds <- vector("list", num.sim); seeds[[1]] <- .Random.seed
for (k in 2:num.sim) seeds[[k]] <- parallel::nextRNGStream(seeds[[k - 1]])



# --- Cluster setup ---

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
#   correlation:      [iid/ar1] from DGPs Erik sent,  https://arxiv.org/pdf/1706.03461.pdf
#   s =               sparsity parameter, number of active covariates 
#   s =               sparsity parameter, number of active covariates in propensity
#                     (and outcome, when outcome_set = "track_s") [Zhao 2019 §6.2]
#   signs =           sign pattern for propensity coefficients [pos/neg/mixed/ks]
#   strength =        outcome signal multiplier (Y = tau*W + strength * m(X) + ε)
#   decay =           coefficient decay exponent: beta_j ∝ sigma_j * j^(-decay).
#                     1 = harmonic (default, pinned outside E6), 0.5 = dense boundary;
#                     approx. sparse design per Belloni et al. (2014, ReStud)
#   outcome_set =     [fixed4/track_s] fixed4 = Tan outcomes on X_1..X_4 (default,
#                     reproduces DGP1). track_s = outcome loadings j^(-decay) on
#                     X_1..X_s, normalised + centered, all s covariates confounders (E6 only)
# Logged diagnostics (per rep):
#   lam_balcv1/0    CV-selected λ, treated/control arm (cv.balnet)
#   lam_end1/0      attained path endpoint λ, treated/control arm
#   lam_glmcv       CV-selected lambda.min (cv.glmnet)
#   trunc05         binary: path failed to reach λ = 0.05 (pos itivity diagnostic)
#   nnz_balcv1/0    nonzero coefficients at CV λ, treated/control arm
#   nnz_glm         nonzero coefficients at lambda.min (glmnet)
#   smd1/0_05       attained max |SMD| at λ = 0.05, treated/control arm (KKT check)
#   smd1/0_cv       attained max |SMD| at CV λ, treated/control arm
#   balnetrate      ATE from λ = √(log p / n), Wager (2024) §7.2 rate
#   err             convergence failure flag (1 = failed, 0 = ok)


# --- DGP 2 ---

# --- DGP2 (new, extended) ---
generate_data <- function(n, p, s = 4, outcome = "linear", misspec = FALSE,
                          covcor = "iid", overlap = 1, signs = "pos",
                          strength = 1, decay = 1, outcome_set = "fixed4") {
  
  if (covcor == "iid") {
    X <- matrix(rnorm(n * p), n, p)
  } else if (covcor == "ar1") {
    Z <- matrix(rnorm(n * p), n, p); X <- Z
    for (j in 2:p) X[, j] <- 0.5 * X[, j - 1] + sqrt(0.75) * Z[, j]
  }
  
  # --- propensity: j^-decay, normalized ---
  beta <- switch(signs,
                 pos   = rep(1, s),
                 neg   = rep(-1, s),
                 mixed = rep(c(1, -1), length.out = s),
                 ks    = rep_len(c(1, -1, 1, 1), s),
                 stop("unknown signs: ", signs)
  ) / seq_len(s)^decay
  
  Sigma_s <- if (covcor == "iid") diag(s) else stats::toeplitz(0.5^seq(0, s - 1))
  eta <- (X[, 1:s, drop = FALSE] %*% beta) /
    (sqrt(as.numeric(crossprod(beta, Sigma_s %*% beta))) * overlap)
  e <- 1 / (1 + exp(eta))
  
  W <- rbinom(n, 1, e)
  tau <- 0
  
  # --- outcome signal ---
  if (outcome_set == "fixed4") {
    m <- switch(outcome,
                linear = X[,1] + 0.5*X[,2] + 0.5*X[,3] + 0.5*X[,4],
                quad1  = rowSums(pmax(X[, 1:4], 0)^2),
                exp    = rowSums(exp(X[, 1:4] / 2)))
  } else {  # track_s
    a  <- 1 / seq_len(s)^decay
    a  <- a / sqrt(as.numeric(crossprod(a, Sigma_s %*% a)))
    Xs <- X[, 1:s, drop = FALSE]
    m  <- switch(outcome,
                 linear = drop(Xs %*% a),
                 quad1  = drop(pmax(Xs, 0)^2 %*% a) - 0.5 * sum(a),
                 exp    = drop(exp(Xs / 2) %*% a) - exp(1/8) * sum(a))
  }
  Y <- tau * W + strength * m + rnorm(n)
  
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
#--- grid
grid <- expand.grid(
  n           = c(1000, 5000),
  pn_ratio    = 0.05,
  s           = 4,
  outcome     = "linear",
  misspec     = FALSE,
  covcor      = "iid",
  overlap     = c(0.25, 0.20, 0.15),
  signs       = "pos",
  strength    = 1,
  decay       = 1,
  outcome_set = "fixed4",
  stringsAsFactors = FALSE
)
grid$p <- grid$n * grid$pn_ratio

# ---- cluster
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, { library(balnet); library(glmnet) })
clusterExport(cl, c("generate_data", "seeds"))

# --- Simulation Function ---

run_sim <- function(n, p, s = 4, outcome = "linear", misspec = FALSE, covcor = "iid",
                    overlap = 1, signs = "pos", strength = 1,
                    decay = 1, outcome_set = "fixed4", num.sim = 200) {
  clusterExport(cl, c("n","p","s","outcome","misspec","covcor","overlap",
                      "signs","strength","decay","outcome_set"),
                envir = environment())
  
  stat_names <- c("balnetcv","balnet0","balnet05","balnet10","balnetrate",
                  "glmnetcv","glmcv_normal","oracle","normalised_oracle",
                  "lam_balcv1","lam_balcv0","lam_end1","lam_end0","lam_glmcv",
                  "trunc05","nnz_balcv1","nnz_balcv0","nnz_glm",
                  "smd1_05","smd0_05","smd1_cv","smd0_cv",
                  "prev","emin","emax","nout05","nout01",
                  "cvloss_cv1","cvloss_cv0","cvloss_end1","cvloss_end0","err")
  
  parSapply(cl, 1:num.sim, function(i) {
    assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    tryCatch({
      data <- generate_data(n, p, s = s, outcome = outcome, misspec = misspec,
                            covcor = covcor, overlap = overlap,
                            signs = signs, strength = strength,
                            decay = decay, outcome_set = outcome_set)
      Y <- data$Y; W <- data$W; X <- data$X
      
      ## overlap / prevalence diagnostics (true e*)
      e_true <- data$e
      prev   <- mean(W)
      emin   <- min(e_true); emax <- max(e_true)
      nout05 <- sum(e_true < 0.05 | e_true > 0.95)
      nout01 <- sum(e_true < 0.01 | e_true > 0.99)
      
      ## balancing: one fit, extended floor
      fit_cv <- cv.balnet(X, W, nfolds = 5, max.imbalance = 1e-4)
      wts_cv <- balweights(fit_cv)
      wts_05 <- balweights(fit_cv, lambda = 0.05)
      wts_10 <- balweights(fit_cv, lambda = 0.10)
      wts_0  <- balweights(fit_cv, lambda = 0)
      lam_rate <- sqrt(log(p) / n)   # 
      
      # lam_end / trunc05 logs catch it
      wts_rate <- balweights(fit_cv, lambda = lam_rate)
      
      tau_balcv   <- mean(Y * (wts_cv$treated   - wts_cv$control))
      tau_bal0    <- mean(Y * (wts_0$treated    - wts_0$control))
      tau_bal05   <- mean(Y * (wts_05$treated   - wts_05$control))
      tau_bal10   <- mean(Y * (wts_10$treated   - wts_10$control))
      tau_balrate <- mean(Y * (wts_rate$treated - wts_rate$control))
      
      ## lambda logging
      lam_balcv1 <- fit_cv$lambda.min$treated
      lam_balcv0 <- fit_cv$lambda.min$control
      lam_end1   <- min(fit_cv$lambda$treated)
      lam_end0   <- min(fit_cv$lambda$control)
      trunc05    <- as.numeric(max(lam_end1, lam_end0) > 0.05)
      
      ## out-of-fold balance loss at lambda_cv vs endpoint (criterion flatness)
      cvget <- function(arm) {
        cv  <- tryCatch(fit_cv$cvm[[arm]],    error = function(e) NULL)
        lam <- tryCatch(fit_cv$lambda[[arm]], error = function(e) NULL)
        if (is.null(cv) || is.null(lam)) return(c(NA_real_, NA_real_))
        c(cv[which.min(abs(lam - fit_cv$lambda.min[[arm]]))], cv[which.min(lam)])
      }
      cv1 <- cvget("treated"); cv0 <- cvget("control")
      
      ## nonzero counts (accessor validated in DGP2)
      nnz_balcv1 <- sum(coef(fit_cv)$treated$betas != 0)
      nnz_balcv0 <- sum(coef(fit_cv)$control$betas != 0)
      
      ## attained max SMD (self-normalized; verify against raw KKT form once)
      xbar <- colMeans(X); sdx <- apply(X, 2, sd)
      wsmd <- function(w) max(abs(colSums(w * X) / sum(w) - xbar) / sdx)
      smd1_05 <- wsmd(wts_05$treated); smd0_05 <- wsmd(wts_05$control)
      smd1_cv <- wsmd(wts_cv$treated); smd0_cv <- wsmd(wts_cv$control)
      
      ## MLE
      fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = 5)
      e_hat <- as.numeric(predict(fit_glm, newx = X, s = "lambda.min", type = "response"))
      lam_glmcv <- fit_glm$lambda.min
      nnz_glm   <- sum(coef(fit_glm, s = "lambda.min")[-1] != 0)
      tau_ht    <- mean(W * Y / e_hat) - mean((1 - W) * Y / (1 - e_hat))
      tau_hajek <- sum(W * Y / e_hat) / sum(W / e_hat) -
        sum((1 - W) * Y / (1 - e_hat)) / sum((1 - W) / (1 - e_hat))
      
      ## oracles
      tau_oracle <- mean(W * Y / data$e) - mean((1 - W) * Y / (1 - data$e))
      tau_oracle_normal <- sum(W * Y / data$e) / sum(W / data$e) -
        sum((1 - W) * Y / (1 - data$e)) / sum((1 - W) / (1 - data$e))
      
      ## order must match stat_names
      setNames(c(tau_balcv, tau_bal0, tau_bal05, tau_bal10, tau_balrate,
                 tau_ht, tau_hajek, tau_oracle, tau_oracle_normal,
                 lam_balcv1, lam_balcv0, lam_end1, lam_end0, lam_glmcv,
                 trunc05, nnz_balcv1, nnz_balcv0, nnz_glm,
                 smd1_05, smd0_05, smd1_cv, smd0_cv,
                 prev, emin, emax, nout05, nout01,
                 cv1[1], cv0[1], cv1[2], cv0[2], 0), stat_names)
    }, error = function(e)
      setNames(c(rep(NA_real_, length(stat_names) - 1), 1), stat_names))
  })
}

#--- Run

out <- list()

for (i in seq_len(nrow(grid))) {
  cat(format(Sys.time()), "| grid", i, "of", nrow(grid), "\n")
  print(grid[i, ])
  
  results <- run_sim(n = grid$n[i], p = grid$p[i], s = grid$s[i],
                     outcome = grid$outcome[i], misspec = grid$misspec[i],
                     covcor = grid$covcor[i], overlap = grid$overlap[i],
                     signs = grid$signs[i], strength = grid$strength[i],
                     decay = grid$decay[i], outcome_set = grid$outcome_set[i],
                     num.sim = num.sim)
  
  out[[i]] <- data.frame(grid[i, ],
                         sim       = rep(seq_len(ncol(results)), each = nrow(results)),
                         estimator = rownames(results),
                         tau_hat   = c(results),
                         row.names = NULL, stringsAsFactors = FALSE)
  
  ## crash insurance: per-cell checkpoint
  saveRDS(out[[i]], sprintf("xlo_cell_%03d.rds", i))
}

stopCluster(cl)
out.df <- do.call(rbind, out)
rownames(out.df) <- NULL
write.csv(out.df, gzfile("extralowoverlap.csv.gz"), row.names = FALSE)
