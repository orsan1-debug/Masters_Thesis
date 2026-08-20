## R/estimators.R -----------------------------------------------------------
### LEGACY MAP
#  est_recode <- c(glmnetcv = "glmnetcv_ht", glmcv_normal = "glmnetcv_hajek", oracle = "oracle_ht", normalised_oracle = "oracle_hajek")

## --- ATE estimators -------------------------------------------------------

ate_ht <- function(Y, W, e1, e0)
  colMeans(Y * (W / e1 - (1 - W) / e0))

ate_hajek <- function(Y, W, e1, e0) {           # Hajek (Chattopadhyay et al. 2020)
  a1 <- W / e1
  a0 <- (1 - W) / e0
  colSums(Y * a1) / sum(a1) - colSums(Y * a0) / sum(a0)
}

ate_bal <- function(Y, w)                       # w = balweights() object
  colMeans(Y * drop(w$treated - w$control))

## --- diagnostics ----------------------------------------------------------

## max SMD of the weighted mean against the full-sample mean.
max_smd <- function(w, X, xbar, sdx)
  max(abs(colSums(drop(w) * X) / sum(w) - xbar) / sdx)

## CV loss at the selected lambda and at the end of the path.
cv_loss <- function(fit, arm) {
  cvi <- fit$`_cv.info`
  c(cvi$cv.mean[[arm]][cvi$idx.min[[arm]]],
    cvi$cv.mean[[arm]][which.min(fit$lambda[[arm]])])
}

## --- per-replication estimation -------------------------------------------

estimate_all <- function(data,
                         lambdas       = c("0" = 0, "05" = 0.05, "10" = 0.10),
                         nfolds        = 5,
                         max_imbalance = 1e-4,   # path floor.
                         cv_curve      = FALSE,
                         lambda_grid   = NULL) {   #check CV against optimal lambda.
  
  Y  <- as.matrix(data$Y)                        # outcomes share (X, W) draws
  W  <- data$W
  X  <- data$X
  n  <- nrow(X); p <- ncol(X); k <- ncol(Y)
  e1 <- data$e
  e0 <- if (is.null(data$e0)) 1 - e1 else data$e0  # only dgp2 supplies e0
  
  ## balancing: one CV fit. lambda floor is max_imbalance, never 0.
  fit_bal <- cv.balnet(X, W, nfolds = nfolds, max.imbalance = max_imbalance)
  w_cv    <- balweights(fit_bal)
  w_fix   <- lapply(lambdas, function(l) balweights(fit_bal, lambda = l))
  w_rate  <- balweights(fit_bal, lambda = sqrt(log(p) / n))  # Wager (2024) s.7.2
  
  ## check CV against optimal lambda: tau on a shared lambda grid, both arms.
  if (!is.null(lambda_grid)) {
    w_grid   <- lapply(lambda_grid, function(l) balweights(fit_bal, lambda = l))
    tau_path <- do.call(rbind, lapply(w_grid, ate_bal, Y = Y))   # grid x outcomes
  }
  
  ## MLE
  fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  e_hat   <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
  
  ## fixed-lambda ATE rows (balnet0, balnet05, balnet10)
  ate_fix <- do.call(rbind, lapply(w_fix, ate_bal, Y = Y))
  rownames(ate_fix) <- paste0("balnet", names(lambdas))
  
  
  
  ## scalar diagnostics (identical across outcome columns)
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((n - 1) / n)   # balnet's standardisation scale
  cv1  <- cv_loss(fit_bal, "treated")
  cv0  <- cv_loss(fit_bal, "control")
  
  smd_fix <- numeric(0)
  for (nm in names(lambdas)) {
    smd_fix[paste0("smd1_", nm)] <- max_smd(w_fix[[nm]]$treated, X, xbar, sdx)
    smd_fix[paste0("smd0_", nm)] <- max_smd(w_fix[[nm]]$control, X, xbar, sdx)
  }
  
  diags <- c(
    lam_balcv1 = fit_bal$lambda.min$treated,
    lam_balcv0 = fit_bal$lambda.min$control,
    lam_end1   = min(fit_bal$lambda$treated),
    lam_end0   = min(fit_bal$lambda$control),
    lam_glmcv  = fit_glm$lambda.min,
    trunc05    = as.numeric(min(fit_bal$lambda$treated) > 0.05 |
                              min(fit_bal$lambda$control) > 0.05),  # bal05 infeasible
    nnz_balcv1 = sum(coef(fit_bal)$treated[-1, ] != 0),
    nnz_balcv0 = sum(coef(fit_bal)$control[-1, ] != 0),
    nnz_glm    = sum(coef(fit_glm, s = "lambda.min")[-1] != 0),
    smd1_cv    = max_smd(w_cv$treated, X, xbar, sdx),
    smd0_cv    = max_smd(w_cv$control, X, xbar, sdx),
    smd_fix,
    prev       = mean(W),                          # true-e overlap diagnostics
    emin       = min(e1),
    emax       = max(e1),
    nout05     = sum(e1 < 0.05 | e1 > 0.95),
    nout01     = sum(e1 < 0.01 | e1 > 0.99),
    cvloss_cv1  = cv1[1], cvloss_cv0  = cv0[1],
    cvloss_end1 = cv1[2], cvloss_end0 = cv0[2]
  )
  ## estimator rows + diagnostics, one matrix per replication
  out <- rbind(
    balnetcv       = ate_bal(Y, w_cv),
    ate_fix,                                      # balnet0, balnet05, balnet10
    balnetrate     = ate_bal(Y, w_rate),
    glmnetcv_ht    = ate_ht(Y, W, e_hat, 1 - e_hat),     # screen only
    glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat),
    oracle_ht      = ate_ht(Y, W, e1, e0),
    oracle_hajek   = ate_hajek(Y, W, e1, e0),
    matrix(diags, length(diags), k, dimnames = list(names(diags), NULL))
  )
  colnames(out) <- colnames(Y)
  
  if (cv_curve)                                   # full CV curve, subset reps only
    attr(out, "cv_curve") <- list(lambda  = fit_bal$lambda,
                                  cv.mean = fit_bal$`_cv.info`$cv.mean)
  if (!is.null(lambda_grid))
    attr(out, "tau_path") <- list(lambda = lambda_grid, tau = tau_path)
  out
}
