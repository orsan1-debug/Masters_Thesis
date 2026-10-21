# CV twin of the basic_dgp selector grids on dgp2's overlap and strength axes. IPW-balnet and ABW
# needs the cv-experimental balnet build (cv.boot.balnet) + glmnet + R/estimators_ipw.R

#' lambda grid per cell from a pilot draw. log spaced to its path end + 1 pt below every path end
#' that last pt clamps each rep to its own path end so the last path row = the fixed arm
make_lam_grid <- function(data, n_lam = 30, alpha = 1, max_imbalance = 1e-4,
                          min_ratio = 1e-4, maxit = 1e4) {
  fit <- do.call(balnet, c(list(data$X, data$W, maxit = maxit),
                           path_args(alpha, max_imbalance, min_ratio)))
  top <- max(vapply(fit$lambda, max, numeric(1)))
  end <- min(vapply(fit$lambda, min, numeric(1)))   # pilot path end
  c(exp(seq(log(top), log(end), length.out = n_lam - 1)),
    if (alpha == 1) max_imbalance else 0)
}

#' penalty and path end args for balnet(). lasso ends at max.imbalance (max |SMD|)
#' balnet() refuses max.imbalance for alpha < 1 (lambda not an SMD) so lambda.min.ratio there
path_args <- function(alpha, max_imbalance = 1e-4, min_ratio = 1e-4)
  if (alpha == 1) list(alpha = 1, max.imbalance = max_imbalance) else
    list(alpha = alpha, lambda.min.ratio = min_ratio)

#' lambda picks of the 5 cv.balnet and cv.boot.balnet selectors from one full fit + the fold fits
#' balnet() has no RNG so separate calls refit the same paths. same picks w/ 11 fits not 35
select_lambdas <- function(X, W, foldid, B = 500, ...) {
  dots     <- list(...)
  full     <- balnet(X, W, ...)
  lam_full <- full[["_lambda"]]
  sw       <- full[["sample.weights"]]
  stopifnot(full[["target"]] == "ATE")

  # SMD ref = full sample means and sds as in both package fns (no sample weights so NULL)
  x_stats <- balnet:::col_stats(X, NULL, compute_sd = TRUE)
  x_stats$scale[x_stats$scale <= 0] <- 1

  # lambda w/ the min mean loss per arm as cv.balnet()
  pick <- function(losses)
    lapply(c(control = "control", treated = "treated"), function(a) {
      m <- matrix(unlist(lapply(losses, `[[`, a)), nrow = length(losses),
                  ncol = length(lam_full[[a]]), byrow = TRUE)
      lam_full[[a]][which.min(colMeans(m))]
    })

  # CV. one training fit per fold scored by all 3 criteria
  cv_loss <- list(cv.bloss = list(), cv.smd = list(), cv.inf = list())
  for (k in seq_len(max(foldid))) {
    test  <- foldid == k
    fit_k <- do.call(balnet, c(list(X = X[!test, , drop = FALSE],
                                    W = W[!test], standardize = ".inplace"),
                               dots, list(sample.weights = sw[!test])))
    held  <- list(fit_k, X[test, , drop = FALSE], W[test], sw[test], lam_full)
    cv_loss$cv.bloss[[k]] <- do.call(balnet:::get_balance_loss, held)
    cv_loss$cv.smd[[k]]   <- do.call(balnet:::get_imbalance,
      c(held, list(X.stats = x_stats, type.measure = "imbalance.mean")))
    cv_loss$cv.inf[[k]]   <- do.call(balnet:::get_imbalance,
      c(held, list(X.stats = x_stats, type.measure = "imbalance.inf")))
  }

  # bootstrap. full fit propensities scored on B random half samples as cv.boot.balnet()
  w_hat_full <- predict(full, X, lambda = lam_full, type = "response",
                        .simplify = FALSE)
  boot_loss <- function(measure) {
    tests <- replicate(B, sample.int(length(W), length(W) %/% 2),
                       simplify = FALSE)
    lapply(tests, function(test)
      balnet:::get_imbalance(full, X[test, , drop = FALSE], W[test], sw[test],
                             lam_full, X.stats = x_stats,
                             type.measure = measure,
                             W.hat = lapply(w_hat_full, function(m)
                               m[test, , drop = FALSE])))
  }
  boot_smd <- boot_loss("imbalance.mean")   # RNG order. mean SMD draws first
  boot_inf <- boot_loss("imbalance.inf")

  list(full = full,
       lambda_min = c(lapply(cv_loss, pick),
                      list(boot.smd = pick(boot_smd),
                           boot.inf = pick(boot_inf))))
}

#' one rep of the CV selector arms on one dgp2() draw. each weight set gives IPW-balnet and ABW
#' RNG order is folds then the 2 cv.boot.balnet subsample draws then fit_or()
estimate_cv <- function(data, lam_grid = data$lam_grid, nfolds = 10,  # 10 = cv.balnet default
                        alpha = 1, max_imbalance = 1e-4, min_ratio = 1e-4,
                        B = 500, or_folds = 5, maxit = 1e4, tol = 1e-7) {
  # maxit 1e4 not the 1e5 default. a truncated path burns extra iters failing at its last lambda
  # on 24 draws at overlap 0.05-0.5 paths and weights matched 1e5 exactly and ~10x faster
  Y    <- as.matrix(data$Y)
  W    <- data$W
  X    <- data$X
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))   # balnet SMD scale

  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X)))  # shared folds
  sel <- do.call(select_lambdas,
                 c(list(X, W, foldid, B = B, maxit = maxit, tol = tol),
                   path_args(alpha, max_imbalance, min_ratio)))
  or <- fit_or(data, nfolds = or_folds)

  # one balnet path serves every weight set. fixed = the path end
  full    <- sel$full
  lam_set <- c(sel$lambda_min,
               list(fixed = list(control = min(full$lambda$control),
                                 treated = min(full$lambda$treated))))
  w_set <- lapply(lam_set, function(l) balweights(full, lambda = l))

  ate_w <- function(w1, w0)                               # Hajek per arm (= HT at solved lambdas)
    colSums(Y * w1) / sum(w1) - colSums(Y * w0) / sum(w0)

  set_rows <- function(w, name) {
    w1  <- drop(w$treated)
    w0  <- drop(w$control)
    out <- rbind(ate_w(w1, w0), ate_aug(Y, w1, w0, or$m1, or$m0))
    rownames(out) <- paste0(c("ipw_", "abw_"), name)
    out
  }
  sets <- do.call(rbind, Map(set_rows, w_set, names(w_set)))

  # path rows. same lambda in both arms. balweights() interpolates and clamps outside the path
  wp  <- balweights(full, lambda = lam_grid)              # n x K per arm
  k   <- seq_along(lam_grid)
  ipw <- do.call(rbind, lapply(k, function(j)
    ate_w(wp$treated[, j], wp$control[, j])))
  abw <- do.call(rbind, lapply(k, function(j)
    ate_aug(Y, wp$treated[, j], wp$control[, j], or$m1, or$m0)))
  rownames(ipw) <- sprintf("ipw_path_%02d", k)
  rownames(abw) <- sprintf("abw_path_%02d", k)

  # selected lambda + max |SMD| + ESS per arm and weight set
  diag_set <- function(w, l, name) {
    d1 <- weight_diags(w$treated, X, xbar, sdx)
    d0 <- weight_diags(w$control, X, xbar, sdx)
    d  <- c(lam1 = l$treated, lam0 = l$control,
            smd1 = d1[["smd"]], smd0 = d0[["smd"]],
            ess1 = d1[["ess"]], ess0 = d0[["ess"]])
    setNames(d, paste0(names(d), "_", name))
  }
  diags <- unlist(unname(Map(diag_set, w_set, lam_set, names(w_set))))

  out <- rbind(sets, ipw, abw,
               matrix(diags, length(diags), ncol(Y),
                      dimnames = list(names(diags), NULL)))
  colnames(out) <- colnames(Y)
  out
}
