## R/estimators.R -----------------------------------------------------------
### LEGACY MAP
#  est_recode <- c(glmnetcv = "glmnetcv_ht", glmcv_normal = "glmnetcv_hajek", oracle = "oracle_ht", normalised_oracle = "oracle_hajek")

## --- ATE estimators -------------------------------------------------------

#' Horvitz-Thompson ATE estimate
#'
#' Inverse-propensity weighted difference in means without normalisation, one
#' estimate per outcome column.
#'
#' @param Y Numeric matrix of outcomes, n x k.
#' @param W Binary treatment vector of length n.
#' @param e1 P(W = 1 | X), length n.
#' @param e0 P(W = 0 | X), length n.
#' @return Numeric vector of length k.
ate_ht <- function(Y, W, e1, e0)
  colMeans(Y * (W / e1 - (1 - W) / e0))

#' Hajek (normalised) ATE estimate
#'
#' Inverse-propensity weighted difference in means with the weights normalised
#' to sum to one within each arm, one estimate per outcome column.
#'
#' @param Y Numeric matrix of outcomes, n x k.
#' @param W Binary treatment vector of length n.
#' @param e1 P(W = 1 | X), length n.
#' @param e0 P(W = 0 | X), length n.
#' @return Numeric vector of length k.
ate_hajek <- function(Y, W, e1, e0) {           # Hajek (Chattopadhyay et al. 2020)
  a1 <- W / e1
  a0 <- (1 - W) / e0
  colSums(Y * a1) / sum(a1) - colSums(Y * a0) / sum(a0)
}

#' ATE from balancing weights
#'
#' Mean of Y times (treated weight minus control weight), one estimate per
#' outcome column.
#'
#' @param Y Numeric matrix of outcomes, n x k.
#' @param w A balweights() object with elements treated and control.
#' @return Numeric vector of length k.
ate_bal <- function(Y, w)                       # w = balweights() object
  colMeans(Y * drop(w$treated - w$control))

## --- diagnostics ----------------------------------------------------------

## max SMD of the weighted mean against the full-sample mean.
#' Maximum standardised mean difference of one weighted arm
#'
#' Largest absolute difference between the weighted column means and the
#' full-sample means, divided by the full-sample standard deviations.
#'
#' @param w Weight vector (or n x 1 matrix) for one arm.
#' @param X Covariate matrix, n x p.
#' @param xbar Full-sample column means of X.
#' @param sdx Full-sample column standard deviations of X.
#' @return A single number.
max_smd <- function(w, X, xbar, sdx)
  max(abs(colSums(drop(w) * X) / sum(w) - xbar) / sdx)

## CV loss at the selected lambda and at the end of the path.
#' Cross-validated balance loss at the selected and at the end-of-path lambda
#'
#' @param fit A cv.balnet fit; its `_cv.info` element is read.
#' @param arm "treated" or "control".
#' @return Numeric vector of length 2: the CV loss at lambda.min and at the
#'   smallest lambda on the path.
cv_loss <- function(fit, arm) {
  cvi <- fit$`_cv.info`
  c(cvi$cv.mean[[arm]][cvi$idx.min[[arm]]],
    cvi$cv.mean[[arm]][which.min(fit$lambda[[arm]])])
}

## --- per-replication estimation -------------------------------------------

#' One replication of every ipw estimator with its diagnostics
#'
#' Fits cv.balnet once and takes weights at the CV-selected lambda, at the
#' fixed lambdas and at sqrt(log(p) / n) (Wager 2024, s.7.2); fits cv.glmnet
#' for the propensity-score arms; adds the true-propensity oracle. The
#' estimates are stacked with scalar diagnostics (selected lambdas, path
#' endpoints, active-set sizes, attained max SMD, true-propensity overlap
#' summaries and the two CV losses), one column per outcome.
#'
#' @param data A draw from dgp1() or dgp2(): list with Y, W, X, e and,
#'   optionally, e0.
#' @param lambdas Named numeric vector of fixed lambdas; the names become the
#'   "balnet<name>" row labels.
#' @param nfolds Number of CV folds for cv.balnet and cv.glmnet.
#' @param max_imbalance Path floor passed to cv.balnet as max.imbalance.
#' @param cv_curve Logical; if TRUE the full CV loss curve is attached as the
#'   "cv_curve" attribute.
#' @param lambda_grid Optional lambda grid; if given, the ATE along the grid is
#'   attached as the "tau_path" attribute.
#' @return A numeric matrix (estimator and diagnostic rows x outcome columns)
#'   with optional attributes cv_curve and tau_path.
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
