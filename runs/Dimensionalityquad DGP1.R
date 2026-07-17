rm(list = ls())


library(balnet)
library(glmnet)
library(parallel)
library(MASS)



# --- Cluster setup ---

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, { library(balnet); library(glmnet); library(MASS) })
clusterSetRNGStream(cl, iseed = 2026)


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



# --- DGP 1 ---

generate_data = function(n, p, outcome = "linear", misspec = FALSE, covcor = "iid", overlap = 1) {
  
  if (covcor == "iid") {
    X <- matrix(rnorm(n * p), n, p)
  } else if (covcor == "ar1") {
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = stats::toeplitz(0.5^seq(0, p - 1)))
  }
  
  e = 1 / (1 + exp((X[, 1] - 0.5*X[, 2] + 0.25*X[,3] + 0.1*X[,4]) / overlap))
  W <- rbinom(n, 1, e)
  tau <- 0
  
  if (outcome == "linear") {
    Y <- tau * W + X[,1] + 0.5*X[,2] + 0.5*X[,3] + 0.5*X[,4] + rnorm(n)
  } else if (outcome == "quad1") {
    Y <- tau * W + rowSums(pmax(X[, 1:4], 0)^2) + rnorm(n)
  } else if (outcome == "exp") {
    Y <- tau * W + rowSums(exp(X[, 1:4] / 2)) + rnorm(n)
  }
  
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



grid <- expand.grid(
  n        = c(500, 1000),
  pn_ratio = c(0.1, 0.5, 1.0, 2.0),
  outcome  = c("quad1"),
  misspec  = FALSE,
  covcor   = "iid",
  overlap  = 1
)
grid$p <- grid$n * grid$pn_ratio

# --- Simulation Function ---

run_sim = function(n, p, outcome = "linear", misspec = FALSE, covcor = "iid", overlap = 1, num.sim = 1000) {
  clusterExport(cl, c("generate_data", "n", "p", "outcome", "misspec", "covcor", "overlap"), envir = environment())
  
  parSapply(cl, 1:num.sim, function(i) {
    data = generate_data(n, p, outcome = outcome, misspec = misspec, covcor = covcor, overlap = overlap)
    Y = data$Y; W = data$W; X = data$X
    #cv balnet
    fit_cv = cv.balnet(X, W, nfolds = 5)
    wts_cv = balweights(fit_cv)
    tau_balcv = mean(Y * (wts_cv$treated - wts_cv$control))
    # balnet lambda .05 and .01 - and 0.
    fit_path = balnet(X, W)
    wts_05 = balweights(fit_path, lambda = 0.05)
    tau_bal05 = mean(Y * (wts_05$treated - wts_05$control))
    wts_10 = balweights(fit_path, lambda = .1)
    tau_bal10 = mean(Y * (wts_10$treated - wts_10$control))
    wts_0 = balweights(fit_path, lambda = 0)
    tau_bal0 = mean(Y * (wts_0$treated - wts_0$control))
    #MLE with glmnet
    fit_glm = cv.glmnet(X, W, family = "binomial", nfolds = 5)
    e_hat = predict(fit_glm, newx = X, s = "lambda.min", type = "response")
    #standard ATE
    tau_ht = mean(W * Y / e_hat) - mean((1 - W) * Y / (1 - e_hat))
    #hajeck, normalised from Chattophady paper
    tau_hajek = sum(W * Y / e_hat) / sum(W / e_hat) - sum((1 - W) * Y / (1 - e_hat)) / sum((1 - W) / (1 - e_hat))
    #oracle ATE
    tau_oracle = mean(W * Y / data$e) - mean((1 - W) * Y / (1 - data$e))
    #oracle normalised from Chattophady paper
    tau_oracle_normal <- sum(W * Y / data$e) / sum(W / data$e) - sum((1 - W) * Y / (1 - data$e)) / sum((1 - W) / (1 - data$e))
    
    
    c(balnetcv = tau_balcv, balnet0 = tau_bal0, balnet05 = tau_bal05, balnet10 = tau_bal10, glmnetcv = tau_ht, 
      glmcv_normal = tau_hajek, oracle = tau_oracle, normalised_oracle = tau_oracle_normal )
  })
}



# --- Run ---

out <- list()

for (i in seq_len(nrow(grid))) {
  cat(format(Sys.time()), "| grid", i, "of", nrow(grid), "\n")
  print(grid[i, ])
  
  results <- run_sim(
    n       = grid$n[i],
    p       = grid$p[i],
    outcome = grid$outcome[i],
    misspec = grid$misspec[i],
    covcor  = grid$covcor[i],
    overlap = grid$overlap[i]
  )
  
  for (sim in seq_len(ncol(results))) {
    out[[length(out) + 1]] <- data.frame(
      n         = grid$n[i],
      p         = grid$p[i],
      outcome   = grid$outcome[i],
      misspec   = grid$misspec[i],
      covcor    = grid$covcor[i],
      overlap = grid$overlap[i],
      sim       = sim,
      estimator = rownames(results),
      tau_hat   = results[, sim],
      stringsAsFactors = FALSE
    )
  }
}

stopCluster(cl)
out.df <- do.call(rbind, out)
rownames(out.df) <- NULL
write.csv(out.df, gzfile("DimensionalityQuadDGP1.csv.gz"), row.names = FALSE)


str(out.df)
head(out.df)
table(out.df$estimator)
table(out.df$n, out.df$outcome)