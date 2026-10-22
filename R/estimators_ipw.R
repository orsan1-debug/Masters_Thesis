# ATE estimators and per draw estimate fns for runs/aipw
# p. refs are pdf pages. balnet docs = balnet 0.0.4 package reference

# ATE estimators

#' Horvitz-Thompson ATE (Ding 2023 s.11.2.2 p. 186)
ate_ht <- function(Y, W, e1, e0)
  colMeans(Y * (W / e1 - (1 - W) / e0))

#' Hajek ATE. weights normalised within arm (Ding 2023 s.11.2.2 p. 186)
ate_hajek <- function(Y, W, e1, e0) {
  a1 <- W / e1
  a0 <- (1 - W) / e0
  colSums(Y * a1) / sum(a1) - colSums(Y * a0) / sum(a0)
}

# ATE from balnet weights (W/e on arm and 0 off). each arm sums to n so HT = Hajek
# (balnet docs balweights and Tan 2020 RCAL eq. 10 p. 8)
ate_bal <- function(Y, w)
  colMeans(Y * drop(w$treated - w$control))

# lasso OR per arm and outcome predicted for all units (Ding 2023 Def. 12.1 p. 197)
# use_true fits on X_true for the PS wrong OR right cell (Kang and Schafer 2007 p. 528)
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
  or0 <- fit_arm(0)                    # arm 1 then 0 for RNG order
  list(m1 = pred(or1), m0 = pred(or0),
       nnz1 = vapply(or1, nnz_lasso, numeric(1)),
       nnz0 = vapply(or0, nnz_lasso, numeric(1)))
}

#' AIPW ATE. mean of m1 - m0 + Hajek weighted residuals per arm (Ding 2023 Problem 12.2 p. 204)
ate_aug <- function(Y, w1, w0, m1, m0) {
  w1 <- drop(w1)
  w0 <- drop(w0)
  colMeans(m1 - m0) +
    colSums((Y - m1) * w1) / sum(w1) - colSums((Y - m0) * w0) / sum(w0)
}

# diagnostics

# one arm's weights 0 off arm. max |SMD| vs full sample mean (Sverdrup and Hastie 2026 eq. 10 p. 6)
# Kish ESS and max weight w/ the arm scaled to sum n (Austin and Stuart 2015 p. 3663)
weight_diags <- function(w, X, xbar, sdx) {
  w <- drop(w)
  c(smd  = max(abs(colSums(w * X) / sum(w) - xbar) / sdx),
    ess  = sum(w)^2 / sum(w^2),
    wmax = max(w) * length(w) / sum(w))
}

#' non-zero coefs of a cv.glmnet fit at lambda.min. no intercept
nnz_lasso <- function(fit) sum(coef(fit, s = "lambda.min")[-1, ] != 0)

# per rep estimation

#' every estimator and diagnostic on one draw. nfolds 5 as Tan 2020 RCAL s.4 p. 20
#' max_imbalance 1e-4 aims at exact balance. path stops earlier if it cant get there (balnet docs)
estimate_all <- function(data, nfolds = 5, max_imbalance = 1e-4) {
  Y  <- as.matrix(data$Y)                        # all outcomes share the X and W draw
  W  <- data$W
  X  <- data$X
  e1 <- data$e
  e0 <- data$e0
  
  time_bal <- system.time(
    fit_bal <- balnet(X, W, max.imbalance = max_imbalance)
  )[["elapsed"]]
  w_bal <- balweights(fit_bal, lambda = 0)       # lambda = 0 is the path end
  
  time_ps <- system.time(
    fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  )[["elapsed"]]
  e_hat <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
  
  time_or <- system.time(or <- fit_or(data, nfolds))[["elapsed"]]
  
  # PS wrong OR right cell. OR refit on X_true on misspec draws only and not timed
  # done after all other RNG use so the rows above keep their draws (Kang and Schafer 2007 p. 528)
  rows_true <- if (identical(X, data$X_true)) NULL else {
    ot <- fit_or(data, nfolds, use_true = TRUE)
    rbind(
      aipw_glmnetcv_ortrue = ate_aug(Y, W / e_hat, (1 - W) / (1 - e_hat),
                                     ot$m1, ot$m0),
      abw_balnet0_ortrue   = ate_aug(Y, w_bal$treated, w_bal$control,
                                     ot$m1, ot$m0)
    )
  }
  
  # per arm weights 0 off arm so length(w) = n
  # true weight SMD is the benchmark (Ben-Michael et al. 2021 p. 18)
  w_arm <- list(
    bal1  = w_bal$treated, bal0  = w_bal$control,
    glm1  = W / e_hat,     glm0  = (1 - W) / (1 - e_hat),
    true1 = W / e1,        true0 = (1 - W) / e0
  )
  xbar <- colMeans(X)
  # population sd to match balnet's SMD scale
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))
  wd   <- vapply(w_arm, function(w) weight_diags(w, X, xbar, sdx), numeric(3))
  b_bal <- coef(fit_bal, lambda = 0)
  
  # smd ess wmax = balance and weight spread per set ie the mechanism behind a win
  # lam_end above max_imbalance = balance not reached so OR has bias to remove
  diags <- c(
    setNames(c(wd), sub("_bal", "", paste(
      rownames(wd), rep(colnames(wd), each = nrow(wd)), sep = "_"))),
    lam_end1 = min(fit_bal$lambda$treated),
    lam_end0 = min(fit_bal$lambda$control),
    nnz_bal1 = sum(b_bal$treated[-1, ] != 0),
    nnz_bal0 = sum(b_bal$control[-1, ] != 0),
    nnz_ps   = nnz_lasso(fit_glm),  # what each ps fit kept
    time_bal = time_bal,  # cost only. 2 arm fits vs 1 (Sverdrup and Hastie 2026 Remark 1 p. 5)
    time_ps  = time_ps,
    time_or  = time_or
  )
  
  out <- rbind(
    balnet0        = ate_bal(Y, w_bal),                 # exact balance wts. method on trial
    glmnetcv_ht    = ate_ht(Y, W, e_hat, 1 - e_hat),    # gap to _hajek = normalisation share
    glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat), # MLE to beat (Tan 2020 RCAL p. 20)
    # MLE wts + lasso OR. can augmentation rescue MLE (Ding 2023 Def. 12.1 p. 197)
    aipw_glmnetcv  = ate_aug(Y, W / e_hat, (1 - W) / (1 - e_hat), or$m1, or$m0),
    # balnet wts + same OR. expected to win in high dim (Bruns-Smith 2025 eq. 7 p. 7)
    abw_balnet0    = ate_aug(Y, w_bal$treated, w_bal$control, or$m1, or$m0),
    oracle_ht      = ate_ht(Y, W, e1, e0),              # true e unnormalised. normalisation cost
    oracle_hajek   = ate_hajek(Y, W, e1, e0),           # true e. the floor for every row
    rows_true,                                          # *_ortrue. 4th DR cell misspec only
    nnz_or1        = or$nnz1,                           # OR sparsity per arm and outcome
    nnz_or0        = or$nnz0,
    matrix(diags, length(diags), ncol(Y), dimnames = list(names(diags), NULL))
  )
  colnames(out) <- colnames(Y)
  out
}





# add-only arms

#' ols unadjusted and oracle AIPW on one draw. oracle_hajek only to check the pairing
#' no RNG so same grid num_sim and seed give the same draws. rows join stored results by n and sim
estimate_add <- function(data) {
  Y  <- as.matrix(data$Y)
  W  <- data$W
  Xc <- cbind(1, data$X)                 # observed cols as the fits see them
  
  # per arm OLS. mean of fitted values over all n is xbar %*% beta
  # (Kang and Schafer 2007 eq. 7 printed p. 530)
  b_ols <- function(a)
    lm.fit(Xc[W == a, ], Y[W == a, , drop = FALSE])$coefficients
  
  out <- rbind(
    ols          = drop(colMeans(Xc) %*% (b_ols(1) - b_ols(0))),
    unadjusted   = colMeans(Y[W == 1, , drop = FALSE]) -
      colMeans(Y[W == 0, , drop = FALSE]),
    oracle_aipw  = data$tau + ate_hajek(Y - data$m - data$tau * W, W,
                                        data$e, data$e0),   # true e and m. efficient var V*
    oracle_hajek = ate_hajek(Y, W, data$e, data$e0)   # same as estimate_all()
  )
  colnames(out) <- colnames(Y)
  out
}



# diagnostic arms

#' per covariate imbalance + arm means + ps fit error on one draw. no outcome model
#' same 2 ps fits in the same order as estimate_all() so the folds and lambda_CV match
estimate_diag <- function(data, nfolds = 5, max_imbalance = 1e-4) {
  Y  <- as.matrix(data$Y)
  W  <- data$W
  X  <- data$X
  e1 <- data$e
  e0 <- data$e0
  
  time_bal <- system.time(                                 # same as estimate_all()
    fit_bal <- balnet(X, W, max.imbalance = max_imbalance)
  )[["elapsed"]]
  w_bal   <- balweights(fit_bal, lambda = 0)
  p_bal   <- predict(fit_bal, X, lambda = 0)               # P(W = 1) per arm fit
  time_ps <- system.time(
    fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  )[["elapsed"]]
  e_hat   <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
  
  w_arm <- list(bal1 = w_bal$treated, bal0 = w_bal$control,
                glm1 = W / e_hat,     glm0 = (1 - W) / (1 - e_hat),
                true1 = W / e1,       true0 = (1 - W) / e0)
  
  # signed SMD of X1-X4 + max |SMD| over the rest per weight set in Xm
  # population sd (Sverdrup and Hastie 2026 eq. 10 p. 6)
  smd_rows <- function(Xm, prefix) {
    xbar <- colMeans(Xm)
    sdx  <- apply(Xm, 2, sd) * sqrt((nrow(Xm) - 1) / nrow(Xm))
    s <- vapply(w_arm, function(w) (colSums(drop(w) * Xm) / sum(w) - xbar) / sdx,
                numeric(ncol(Xm)))                                  # p x 6
    c(setNames(c(s[1:4, ]),
               paste0(prefix, "_", rep(colnames(s), each = 4), "_", 1:4)),
      setNames(apply(abs(s[-(1:4), , drop = FALSE]), 2, max),
               paste0(prefix, "rest_", colnames(s))))
  }
  
  # Tan's mse and msre per arm. arm sample reweighted by 1 / true e (Tan 2020 RCAL s.4 p. 24)
  fit_err <- function(p_hat, e, on) c(
    mse  = mean(on / e * (p_hat - e)^2),
    msre = mean(on / e * (e / p_hat - 1)^2))
  err <- cbind(bal1 = fit_err(p_bal$treated,     e1, W),
               bal0 = fit_err(1 - p_bal$control, e0, 1 - W),
               glm1 = fit_err(e_hat,             e1, W),
               glm0 = fit_err(1 - e_hat,         e0, 1 - W))   # 2 x 4
  
  # max |SMD| ESS wmax per weight set. same names as estimate_all()
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))
  wd   <- vapply(w_arm, function(w) weight_diags(w, X, xbar, sdx), numeric(3))
  b_bal <- coef(fit_bal, lambda = 0)
  
  diags <- c(
    setNames(c(wd), sub("_bal", "", paste(
      rownames(wd), rep(colnames(wd), each = nrow(wd)), sep = "_"))),
    smd_rows(X, "smd"),
    if (!identical(X, data$X_true)) smd_rows(data$X_true, "smdtrue"),
    setNames(c(err), paste0(rownames(err), "_", rep(colnames(err), each = 2))),
    lam_end1 = min(fit_bal$lambda$treated),
    lam_end0 = min(fit_bal$lambda$control),
    nnz_bal1 = sum(b_bal$treated[-1, ] != 0),
    nnz_bal0 = sum(b_bal$control[-1, ] != 0),
    nnz_ps   = nnz_lasso(fit_glm),
    time_bal = time_bal,
    time_ps  = time_ps
  )
  
  # Hajek arm means per weight set (k x 6) as rows
  mu <- t(vapply(w_arm, function(w) colSums(Y * drop(w)) / sum(w), numeric(ncol(Y))))
  
  out <- rbind(                                       # first 3 rows pair w/ the stored files
    balnet0        = ate_bal(Y, w_bal),
    glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat),
    oracle_hajek   = ate_hajek(Y, W, e1, e0),
    `rownames<-`(mu, paste0("mu_", rownames(mu))),
    mbar           = colMeans(as.matrix(data$m)),     # mean of true E[Y(0) | X]. ref for mu rows
    matrix(diags, length(diags), ncol(Y), dimnames = list(names(diags), NULL))
  )
  colnames(out) <- colnames(Y)
  out
}



# tolerance arms

#' balnet0 and abw_balnet0 at several path end tolerances on one draw
#' balnet has no RNG so RNG order = estimate_all() and pairing rows match the stored file
estimate_tol <- function(data, nfolds = 5,
                         tols = 10^-c(2, 2.5, 3, 3.5, 4.5, 5, 5.5, 6, 6.5, 7)) {
  # tols skip 1e-4 as thats the main file
  Y    <- as.matrix(data$Y)
  W    <- data$W
  X    <- data$X
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))
  
  fits <- lapply(tols, function(t) {
    time_bal <- system.time(
      fit <- balnet(X, W, max.imbalance = t)
    )[["elapsed"]]
    list(fit = fit, w = balweights(fit, lambda = 0), time = time_bal)
  })
  
  fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
  e_hat   <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
  or      <- fit_or(data, nfolds)
  ot      <- if (identical(X, data$X_true)) NULL else
    fit_or(data, nfolds, use_true = TRUE)          # misspec draws only
  
  tol_rows <- function(f, t) {
    b <- coef(f$fit, lambda = 0)
    diags <- c(setNames(weight_diags(f$w$treated, X, xbar, sdx),
                        c("smd1", "ess1", "wmax1")),
               setNames(weight_diags(f$w$control, X, xbar, sdx),
                        c("smd0", "ess0", "wmax0")),
               lam_end1 = min(f$fit$lambda$treated),
               lam_end0 = min(f$fit$lambda$control),
               nnz_bal1 = sum(b$treated[-1, ] != 0),
               nnz_bal0 = sum(b$control[-1, ] != 0),
               time_bal = f$time)
    out <- rbind(
      balnet0            = ate_bal(Y, f$w),
      abw_balnet0        = ate_aug(Y, f$w$treated, f$w$control, or$m1, or$m0),
      abw_balnet0_ortrue = if (!is.null(ot))
        ate_aug(Y, f$w$treated, f$w$control, ot$m1, ot$m0),
      matrix(diags, length(diags), ncol(Y),
             dimnames = list(names(diags), NULL)))
    rownames(out) <- paste0(rownames(out), "_tol", round(-log10(t), 1))
    out
  }
  
  out <- rbind(
    do.call(rbind, Map(tol_rows, fits, tols)),
    aipw_glmnetcv        = ate_aug(Y, W / e_hat, (1 - W) / (1 - e_hat),
                                   or$m1, or$m0),
    aipw_glmnetcv_ortrue = if (!is.null(ot))
      ate_aug(Y, W / e_hat, (1 - W) / (1 - e_hat), ot$m1, ot$m0),
    oracle_hajek         = ate_hajek(Y, W, data$e, data$e0))
  colnames(out) <- colnames(Y)
  out
}