# study 2 (Methods 3.4 and 3.5.1). one full balnet fit per rep gives endpoint selectors and path
# draws in the order of R/estimators_cv_dgp2.R so stored results match (test_selectors.R)

library(balnet)
library(glmnet)

selector_names <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")

#' balnet() args for penalty and path end. lasso to imbalance 1e-4. else a fraction of lambda_max
path_settings <- function(alpha = 1, max_imbalance = 1e-4,
                          min_ratio = if (alpha == 0) 1e-4 else 1e-2) {
  if (alpha == 1) list(alpha = 1, max.imbalance = max_imbalance)
  else list(alpha = alpha, lambda.min.ratio = min_ratio)
}

#' a cell's fixed lambda grid from a pilot draw. 29 points log spaced lambda_max to the path end
#' then one point at or below every path end
pilot_lambda_grid <- function(data, n_lambda = 30, alpha = 1, max_imbalance = 1e-4,
                              min_ratio = if (alpha == 0) 1e-4 else 1e-2, maxit = 1e4) {
  fit <- do.call(balnet, c(list(data$X, data$W, maxit = maxit),
                           path_settings(alpha, max_imbalance, min_ratio)))
  top <- max(vapply(fit$lambda, max, numeric(1)))
  end <- min(vapply(fit$lambda, min, numeric(1)))
  c(exp(seq(log(top), log(end), length.out = n_lambda - 1)), if (alpha == 1) max_imbalance else 0)
}

#' the five selectors' lambdas per arm from one full fit and one set of fold fits
#' scores w/ balnet's own internals (as cv.balnet() and cv.boot.balnet()) so each path is fit once
select_lambdas <- function(X, W, foldid, B = 500, ...) {
  settings <- list(...)
  full     <- balnet(X, W, ...)
  lambdas  <- full[["_lambda"]]
  sw       <- full[["sample.weights"]]
  stopifnot(full[["target"]] == "ATE")
  x_stats <- balnet:::col_stats(X, NULL, compute_sd = TRUE)
  x_stats$scale[x_stats$scale <= 0] <- 1

  # lambda w/ the lowest mean score per arm
  lowest <- function(scores)
    lapply(c(control = "control", treated = "treated"), function(arm) {
      m <- matrix(unlist(lapply(scores, `[[`, arm)), nrow = length(scores),
                  ncol = length(lambdas[[arm]]), byrow = TRUE)
      lambdas[[arm]][which.min(colMeans(m))]
    })

  # cv selectors refit w/o each fold and score every lambda on the held out fold
  cv <- list(cv.bloss = list(), cv.smd = list(), cv.inf = list())
  for (k in seq_len(max(foldid))) {
    held  <- foldid == k
    fit_k <- do.call(balnet, c(list(X = X[!held, , drop = FALSE], W = W[!held], standardize = ".inplace"),
                               settings, list(sample.weights = sw[!held])))
    args <- list(fit_k, X[held, , drop = FALSE], W[held], sw[held], lambdas)
    cv$cv.bloss[[k]] <- do.call(balnet:::get_balance_loss, args)
    cv$cv.smd[[k]]   <- do.call(balnet:::get_imbalance, c(args, list(X.stats = x_stats, type.measure = "imbalance.mean")))
    cv$cv.inf[[k]]   <- do.call(balnet:::get_imbalance, c(args, list(X.stats = x_stats, type.measure = "imbalance.inf")))
  }

  # boot selectors score the full fit on B random half samples
  w_hat <- predict(full, X, lambda = lambdas, type = "response", .simplify = FALSE)
  bootstrap <- function(measure) {
    halves <- replicate(B, sample.int(length(W), length(W) %/% 2), simplify = FALSE)
    lapply(halves, function(h)
      balnet:::get_imbalance(full, X[h, , drop = FALSE], W[h], sw[h], lambdas, X.stats = x_stats,
                             type.measure = measure, W.hat = lapply(w_hat, function(m) m[h, , drop = FALSE])))
  }
  boot.smd <- bootstrap("imbalance.mean")        # mean SMD half samples are drawn first
  boot.inf <- bootstrap("imbalance.inf")

  list(full = full,
       lambda = c(lapply(cv, lowest), list(boot.smd = lowest(boot.smd), boot.inf = lowest(boot.inf))))
}

#' every study 2 arm on one draw. Y has one column per outcome and strength (share_strengths())
#' maxit 1e4 = same paths as balnet's default on 24 draws at c = 0.05 to 0.5 in 1/10 the time
estimate_study2 <- function(data, lambda_grid, nfolds = 10, alpha = 1, max_imbalance = 1e-4,
                            min_ratio = if (alpha == 0) 1e-4 else 1e-2,
                            B = 500, outcome_folds = 5, maxit = 1e4, tol = 1e-7) {
  Y <- as.matrix(data$Y); W <- data$W; X <- data$X
  xbar <- colMeans(X)
  sdx  <- apply(X, 2, sd) * sqrt((nrow(X) - 1) / nrow(X))

  foldid   <- sample(rep(seq_len(nfolds), length.out = nrow(X)))
  selected <- do.call(select_lambdas, c(list(X, W, foldid, B = B, maxit = maxit, tol = tol),
                                        path_settings(alpha, max_imbalance, min_ratio)))
  om   <- fit_outcome_models(X, Y, W, outcome_folds)
  full <- selected$full
  lambda_sets <- c(selected$lambda,
                   list(fixed = list(control = min(full$lambda$control), treated = min(full$lambda$treated))))
  weight_sets <- lapply(lambda_sets, function(l) balweights(full, lambda = l))

  both <- function(w1, w0) rbind(ate_weighted(Y, w1, w0), ate_augmented(Y, w1, w0, om$m1, om$m0))
  by_set <- Map(function(w, name) `rownames<-`(both(w$treated, w$control), paste0(c("ipw_", "abw_"), name)),
                weight_sets, names(weight_sets))

  on_grid <- balweights(full, lambda = lambda_grid)
  k <- seq_along(lambda_grid)
  ipw_path <- t(sapply(k, function(j) ate_weighted(Y, on_grid$treated[, j], on_grid$control[, j])))
  abw_path <- t(sapply(k, function(j) ate_augmented(Y, on_grid$treated[, j], on_grid$control[, j], om$m1, om$m0)))
  rownames(ipw_path) <- sprintf("ipw_path_%02d", k)
  rownames(abw_path) <- sprintf("abw_path_%02d", k)

  diagnostics <- unlist(unname(Map(function(w, l, name) {
    d1 <- weight_diagnostics(w$treated, X, xbar, sdx)
    d0 <- weight_diagnostics(w$control, X, xbar, sdx)
    setNames(c(l$treated, l$control, d1[["smd"]], d0[["smd"]], d1[["ess"]], d0[["ess"]]),
             paste0(c("lam1", "lam0", "smd1", "smd0", "ess1", "ess0"), "_", name))
  }, weight_sets, lambda_sets, names(weight_sets))))

  out <- rbind(do.call(rbind, by_set), ipw_path, abw_path,
               matrix(diagnostics, length(diagnostics), ncol(Y), dimnames = list(names(diagnostics), NULL)))
  colnames(out) <- colnames(Y)
  out
}

#' Y at every signal strength from one draw. kappa m + the same noise. cols like linear_s0.5
share_strengths <- function(data, strengths = c(1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0)) {
  eps <- data$Y - data$m
  data$Y <- do.call(cbind, lapply(strengths, function(kappa)
    `colnames<-`(kappa * data$m + eps, paste0(colnames(data$m), "_s", kappa))))
  data
}
