## R/estimators_ipw.R -------------------------------------------------------
## References (Author Year, p. N) cite the PDF page of the project copy;
## [PKG] = balnet 0.0.4 package reference; [GK] = general knowledge.

## --- ATE estimators -------------------------------------------------------

#' Horvitz-Thompson ATE (Ding 2023 s.11.2.2, p. 186)
#'
#' @param Y Numeric matrix of outcomes, n x k.
#' @param W Binary treatment vector of length n.
#' @param e1 P(W = 1 | X), length n.
#' @param e0 P(W = 0 | X), length n.
#' @return Numeric vector of length k.
ate_ht <- function(Y, W, e1, e0)
  colMeans(Y * (W / e1 - (1 - W) / e0))

#' Hajek ATE, weights normalised within arm (Ding 2023 s.11.2.2, p. 186)
#'
#' @inheritParams ate_ht
#' @return Numeric vector of length k.
ate_hajek <- function(Y, W, e1, e0) {
  a1 <- W / e1
  a0 <- (1 - W) / e0
  colSums(Y * a1) / sum(a1) - colSums(Y * a0) / sum(a0)
}

# ATE from balnet weights, W/e on-arm and 0 off-arm; each arm sums to n, so
# HT = Hajek ([PKG] balweights; Tan 2020 RCAL eq. 10, p. 8).
#   Y: numeric matrix of outcomes, n x k
#   w: balweights() object with elements treated and control
ate_bal <- function(Y, w)
  colMeans(Y * drop(w$treated - w$control))

# Lasso outcome model per arm and outcome, predicted for all units (Ding 2023
# Def. 12.1, p. 197).
#   data:     a draw from dgp2()
#   nfolds:   CV folds for cv.glmnet
#   use_true: fit on X_true, the PS-wrong / OR-correct cell, estimator-side
#             as in Kang & Schafer 2007, p. 528
# Returns m1, m0 (n x k fitted means) and nnz1, nnz0 (non-zero coefficients
# per outcome, intercept excluded).
fit_or <- function(data, nfolds = 5, use_true = FALSE) {
  X <- if (use_true) data$X_true else data$X
  Y <- as.matrix(data$Y)
  fit_arm <- function(a) {
    idx <- data$W == a
    lapply(seq_len(ncol(Y)), function(j)
      cv.glmnet(X[idx, ], Y[idx, j], family = "gaussian", nfolds = nfolds))
  }
  pred <- function(fits)
    vapply(fits, function(f) predict(f, newx = X, s = "lambda.min")[, 1],
           numeric(nrow(X)))
  or1 <- fit_arm(1)
  or0 <- fit_arm(0)                    # arm 1 then 0: RNG order
  list(m1 = pred(or1), m0 = pred(or0),
       nnz1 = vapply(or1, nnz_lasso, numeric(1)),
       nnz0 = vapply(or0, nnz_lasso, numeric(1)))
}

#' AIPW ATE, plug-in mean of m1 - m0 plus Hajek-weighted residuals per arm
#' (Ding 2023 Problem 12.2, p. 204)
#'
#' @param Y Numeric matrix of outcomes, n x k.
#' @param w1,w0 Per-unit weights for each arm, zero off-arm.
#' @param m1,m0 Fitted outcome means, n x k.
#' @return Numeric vector of length k.
ate_aug <- function(Y, w1, w0, m1, m0) {
  w1 <- drop(w1)
  w0 <- drop(w0)
  colMeans(m1 - m0) +
    colSums((Y - m1) * w1) / sum(w1) - colSums((Y - m0) * w0) / sum(w0)
}

## --- diagnostics ----------------------------------------------------------

# One arm's weights, zero off-arm: max |SMD| against the full-sample mean
# (Sverdrup & Hastie 2026 eq. 10, p. 6), Kish ESS [GK], and the max weight
# with the arm scaled to sum n (Austin & Stuart 2015, p. 3663).
#   w:         weight vector or n x 1 matrix
#   X:         covariate matrix, n x p
#   xbar, sdx: full-sample column means and standard deviations of X
# Returns a named vector: smd, ess, wmax.
weight_diags <- function(w, X, xbar, sdx) {
  w <- drop(w)
  c(smd  = max(abs(colSums(w * X) / sum(w) - xbar) / sdx),
    ess  = sum(w)^2 / sum(w^2),
    wmax = max(w) * length(w) / sum(w))
}

#' Non-zero coefficients of a cv.glmnet fit at lambda.min, intercept
#' excluded [GK]
#'
#' @param fit A cv.glmnet object.
#' @return A single number.
nnz_lasso <- function(fit) sum(coef(fit, s = "lambda.min")[-1, ] != 0)

## --- per-replication estimation -------------------------------------------

#' One replication of every estimator and diagnostic on one draw; rows are
#' documented where they are built below
#'
#' @param data A draw from dgp2().
#' @param nfolds CV folds for cv.glmnet; 5 as Tan 2020 RCAL s.4, p. 20.
#' @param max_imbalance balnet path floor (max.imbalance); 1e-4 targets exact
#'   balance, the path truncating earlier where that is unattainable ([PKG]
#'   balnet(), role of lambda).
#' @return Numeric matrix, estimator and diagnostic rows x outcome columns.
estimate_all <- function(data, nfolds = 5, max_imbalance = 1e-4) {
  Y  <- as.matrix(data$Y)                        # outcomes share (X, W) draws
  W  <- data$W
  X  <- data$X
  e1 <- data$e
  e0 <- data$e0
  
  time_bal <- system.time(
    fit_bal <- balnet(X, W, max.imbalance = max_imbalance)
  )[["elapsed"]]
  w_bal <- balweights(fit_bal, lambda = 0)       # lambda = 0 selects the path end
  
  time_ps <- system.time(
    fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  )[["elapsed"]]
  e_hat <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
  
  time_or <- system.time(or <- fit_or(data, nfolds))[["elapsed"]]
  
  # per-arm weights, zero off-arm so length(w) = n; the true-weight SMD is the
  # benchmark (Ben-Michael et al. 2021, p. 18)
  w_arm <- list(
    bal1  = w_bal$treated, bal0  = w_bal$control,        # [PKG] balweights
    glm1  = W / e_hat,     glm0  = (1 - W) / (1 - e_hat),
    true1 = W / e1,        true0 = (1 - W) / e0
  )
  xbar <- colMeans(X)
  # population sd: balnet's SMD scale ([PKG] role of lambda)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))
  wd   <- vapply(w_arm, function(w) weight_diags(w, X, xbar, sdx), numeric(3))
  b_bal <- coef(fit_bal, lambda = 0)
  
  # Diagnostics and what they are for (Q1-Q8: research questions in
  # project_state)
  # smd/ess/wmax: balance and weight spread per weight set, the mechanism
  #   behind a win (Q3, Q6, Q8; Sverdrup & Hastie 2026 eq. 10, p. 6; Austin &
  #   Stuart 2015, p. 3663)
  # lam_end: where the path stopped; above max_imbalance balance was not
  #   reached and the OR has bias to remove (Q4, Q5; [PKG] role of lambda)
  # nnz_bal, nnz_ps: what each PS fit kept (Q5; [PKG] print.balnet; [GK])
  # time: two arm fits against one, cost only (Sverdrup & Hastie 2026
  #   Remark 1, p. 5)
  
  diags <- c(
    setNames(c(wd), sub("_bal", "", paste(
      rownames(wd), rep(colnames(wd), each = nrow(wd)), sep = "_"))),
    lam_end1 = min(fit_bal$lambda$treated),
    lam_end0 = min(fit_bal$lambda$control),
    nnz_bal1 = sum(b_bal$treated[-1, ] != 0),
    nnz_bal0 = sum(b_bal$control[-1, ] != 0),
    nnz_ps   = nnz_lasso(fit_glm),
    time_bal = time_bal,
    time_ps  = time_ps,
    time_or  = time_or
  )
  
  # Estimators and what they are for (Q1-Q8: research questions in
  # project_state)
  # balnet0: exact-balance weights, the method on trial (Q1, Q4, Q6; [PKG]
  #   balweights)
  # glmnetcv_hajek: the MLE it is judged against, ratio form (Q1-3, Q6;
  #   Tan 2020 RCAL p. 20)
  # glmnetcv_ht: the same weights unnormalised, so the gap to _hajek is the
  #   normalisation share (Q7; Imai & Ratkovic 2014 Table 1, p. 11)
  # aipw_glmnetcv: MLE weights + lasso OR, can augmentation rescue MLE
  #   (Q1, Q8; Ding 2023 Def. 12.1, p. 197)
  # abw_balnet0: balnet weights + the same OR, the arm Q5 expects to win in
  #   high dimension (Q1, Q5, Q8; Bruns-Smith 2025 eq. 7, p. 7)
  # oracle_hajek: true e, the floor every row is measured from [GK]
  # oracle_ht: true e unnormalised, normalisation cost alone (Q7; Ding 2023
  #   p. 186)
  # nnz_or: OR sparsity per arm and outcome (Q5; Bruns-Smith 2025 pp. 12-13)
  
  out <- rbind(
    balnet0        = ate_bal(Y, w_bal),
    glmnetcv_ht    = ate_ht(Y, W, e_hat, 1 - e_hat),
    glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat),
    aipw_glmnetcv  = ate_aug(Y, W / e_hat, (1 - W) / (1 - e_hat), or$m1, or$m0),
    abw_balnet0    = ate_aug(Y, w_bal$treated, w_bal$control, or$m1, or$m0),
    oracle_ht      = ate_ht(Y, W, e1, e0),
    oracle_hajek   = ate_hajek(Y, W, e1, e0),
    nnz_or1        = or$nnz1,
    nnz_or0        = or$nnz0,
    matrix(diags, length(diags), ncol(Y), dimnames = list(names(diags), NULL))
  )
  colnames(out) <- colnames(Y)
  out
}