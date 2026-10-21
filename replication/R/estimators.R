# study 1 estimators (Methods 3.3). balnet and MLE weights. ABW and AIPW add a lasso outcome model
# CV folds drawn in the order of R/estimators_ipw.R so stored results match (test_estimators.R)

library(balnet)
library(glmnet)

#' balancing weights at the path end (Sverdrup and Hastie 2026). each arm sums to n
#' the end is max imbalance 1e-4 or earlier where the path stops short
fit_balancing <- function(X, W, max_imbalance = 1e-4) {
  fit <- balnet(X, W, max.imbalance = max_imbalance)
  b   <- coef(fit, lambda = 0)                 # lambda = 0 selects the path endpoint
  list(fit        = fit,
       weights    = balweights(fit, lambda = 0),
       lambda_end = c(treated = min(fit$lambda$treated), control = min(fit$lambda$control)),
       nonzero    = c(treated = sum(b$treated[-1, ] != 0), control = sum(b$control[-1, ] != 0)))
}

#' MLE weights from a 5 fold CV lasso logistic reg at lambda.min (Friedman et al. 2010)
fit_mle <- function(X, W, nfolds = 5) {
  fit <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  list(fit     = fit,
       e_hat   = predict(fit, newx = X, s = "lambda.min", type = "response")[, 1],
       nonzero = nonzero_at_min(fit))
}

#' shared outcome model. one lasso fit per arm and outcome. treated arm first
fit_outcome_models <- function(X, Y, W, nfolds = 5) {
  Y <- as.matrix(Y)
  fit_arm <- function(arm) lapply(seq_len(ncol(Y)), function(j)
    cv.glmnet(X[W == arm, , drop = FALSE], Y[W == arm, j], nfolds = nfolds))
  predict_all <- function(fits)
    vapply(fits, function(f) predict(f, newx = X, s = "lambda.min")[, 1], numeric(nrow(X)))
  treated  <- fit_arm(1)
  controls <- fit_arm(0)
  list(m1 = predict_all(treated), m0 = predict_all(controls),
       nonzero1 = vapply(treated,  nonzero_at_min, numeric(1)),
       nonzero0 = vapply(controls, nonzero_at_min, numeric(1)))
}

nonzero_at_min <- function(fit) sum(coef(fit, s = "lambda.min")[-1, ] != 0)

# ATE estimators. one value per outcome column of Y
ate_horvitz_thompson <- function(Y, W, e1, e0) colMeans(Y * (W / e1 - (1 - W) / e0))

ate_hajek <- function(Y, W, e1, e0) ate_weighted(Y, W / e1, (1 - W) / e0)

ate_weighted <- function(Y, w1, w0) {
  w1 <- drop(w1); w0 <- drop(w0)
  colSums(Y * w1) / sum(w1) - colSums(Y * w0) / sum(w0)
}

#' augmented est. AIPW w/ MLE weights and ABW w/ balancing weights
ate_augmented <- function(Y, w1, w0, m1, m0) {
  w1 <- drop(w1); w0 <- drop(w0)
  colMeans(m1 - m0) + colSums((Y - m1) * w1) / sum(w1) - colSums((Y - m0) * w0) / sum(w0)
}

#' max |SMD| (population sd) ESS = (sum w)^2 / sum w^2 and max weight w/ the arm scaled to n
weight_diagnostics <- function(w, X, xbar, sdx) {
  w <- drop(w)
  c(smd  = max(abs(colSums(w * X) / sum(w) - xbar) / sdx),
    ess  = sum(w)^2 / sum(w^2),
    wmax = max(w) * length(w) / sum(w))
}

#' every study 1 estimator and diagnostic on one draw
estimate_study1 <- function(data, nfolds = 5, max_imbalance = 1e-4) {
  Y <- as.matrix(data$Y); W <- data$W; X <- data$X
  time_bal <- system.time(bal <- fit_balancing(X, W, max_imbalance))[["elapsed"]]
  time_ps  <- system.time(mle <- fit_mle(X, W, nfolds))[["elapsed"]]
  time_or  <- system.time(om  <- fit_outcome_models(X, Y, W, nfolds))[["elapsed"]]
  w_bal <- bal$weights
  e_hat <- mle$e_hat

  # misspec draws only. outcome model refit on X_true after every other fold draw (_ortrue rows)
  refit <- NULL
  if (!identical(X, data$X_true)) {
    om_true <- fit_outcome_models(data$X_true, Y, W, nfolds)
    refit <- rbind(
      aipw_glmnetcv_ortrue = ate_augmented(Y, W / e_hat, (1 - W) / (1 - e_hat), om_true$m1, om_true$m0),
      abw_balnet0_ortrue   = ate_augmented(Y, w_bal$treated, w_bal$control, om_true$m1, om_true$m0))
  }

  weight_sets <- list(bal1 = w_bal$treated, bal0 = w_bal$control,
                      glm1 = W / e_hat,     glm0 = (1 - W) / (1 - e_hat),
                      true1 = W / data$e,   true0 = (1 - W) / data$e0)
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))
  wd   <- vapply(weight_sets, function(w) weight_diagnostics(w, X, xbar, sdx), numeric(3))
  diagnostics <- c(
    setNames(c(wd), sub("_bal", "", paste(rownames(wd), rep(colnames(wd), each = 3), sep = "_"))),
    lam_end1 = bal$lambda_end[["treated"]], lam_end0 = bal$lambda_end[["control"]],
    nnz_bal1 = bal$nonzero[["treated"]],    nnz_bal0 = bal$nonzero[["control"]],
    nnz_ps   = mle$nonzero,
    time_bal = time_bal, time_ps = time_ps, time_or = time_or)

  out <- rbind(
    balnet0        = colMeans(Y * drop(w_bal$treated - w_bal$control)),   # arms sum to n so Hajek
    glmnetcv_ht    = ate_horvitz_thompson(Y, W, e_hat, 1 - e_hat),
    glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat),
    aipw_glmnetcv  = ate_augmented(Y, W / e_hat, (1 - W) / (1 - e_hat), om$m1, om$m0),
    abw_balnet0    = ate_augmented(Y, w_bal$treated, w_bal$control, om$m1, om$m0),
    oracle_ht      = ate_horvitz_thompson(Y, W, data$e, data$e0),
    oracle_hajek   = ate_hajek(Y, W, data$e, data$e0),
    refit,
    nnz_or1        = om$nonzero1,
    nnz_or0        = om$nonzero0,
    matrix(diagnostics, length(diagnostics), ncol(Y), dimnames = list(names(diagnostics), NULL)))
  colnames(out) <- colnames(Y)
  out
}

#' OLS unadjusted diff Oracle AIPW and Oracle Hajek (the add-on files). nothing random
estimate_reference <- function(data) {
  Y  <- as.matrix(data$Y); W <- data$W
  Xc <- cbind(1, data$X)
  ols_arm <- function(arm) lm.fit(Xc[W == arm, , drop = FALSE], Y[W == arm, , drop = FALSE])$coefficients
  ols <- if (min(table(W)) < ncol(Xc)) rep(NA_real_, ncol(Y)) else   # fewer units than coefs
    drop(colMeans(Xc) %*% (ols_arm(1) - ols_arm(0)))
  out <- rbind(
    ols          = ols,
    unadjusted   = colMeans(Y[W == 1, , drop = FALSE]) - colMeans(Y[W == 0, , drop = FALSE]),
    oracle_aipw  = data$tau + ate_hajek(Y - data$m - data$tau * W, W, data$e, data$e0),
    oracle_hajek = ate_hajek(Y, W, data$e, data$e0))
  colnames(out) <- colnames(Y)
  out
}

#' every study 1 row of one draw. OLS unadjusted and Oracle AIPW then estimate_study1()
estimate_study1_with_reference <- function(data, ...)
  rbind(estimate_reference(data)[c("ols", "unadjusted", "oracle_aipw"), , drop = FALSE],
        estimate_study1(data, ...))

# results row name -> name printed in the thesis (arm_labs in handoff/R/thesis.R)
estimator_labels <- c(
  balnet0              = "balnet",
  glmnetcv_hajek       = "glmnet (H\u00e1jek)",
  glmnetcv_ht          = "glmnet (HT)",
  abw_balnet0          = "ABW-balnet",
  aipw_glmnetcv        = "AIPW-glmnet",
  oracle_hajek         = "Oracle (H\u00e1jek)",
  oracle_ht            = "Oracle (HT)",
  abw_balnet0_ortrue   = "ABW-balnet (true OR)",
  aipw_glmnetcv_ortrue = "AIPW-glmnet (true OR)",
  ols                  = "OLS",
  unadjusted           = "Unadjusted",
  oracle_aipw          = "Oracle AIPW")
