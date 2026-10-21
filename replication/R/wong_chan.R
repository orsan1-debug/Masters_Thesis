# Wong and Chan check (Methods 3.5 App H). their 2018 design as run by Wang and Zubizarreta 2020
# draws in the order of R/dgp.R (dgp_wc_overlap) and runs/cv/wz_replication/wc_helpers.R

library(balnet)
library(glmnet)

wc_max_imbalance <- 1e-4
wc_solver <- list(maxit = 1e4, tol = 1e-5)
wc_master_seed <- 20260903

#' one draw. overlap multiplies the logit. sigma = noise sd (both added to match our misspec grid)
#' outcome A has individual effects 1.5 g (tau_i) and B has no effect
generate_wong_chan <- function(n, overlap = 1, sigma = 1) {
  Z     <- matrix(rnorm(n * 10), n, 10)
  logit <- overlap * (-Z[, 1] - 0.1 * Z[, 4])
  e     <- plogis(logit)
  W     <- rbinom(n, 1, e)
  g     <- 27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]
  means <- list(A = 210 + (1.5 * W - 0.5) * g,
                B = Z[, 1] * Z[, 2]^3 * Z[, 3]^2 * Z[, 4] + Z[, 4] * abs(Z[, 1])^0.5)
  eps <- sigma * rnorm(n)
  Y   <- vapply(c("A", "B"), function(o) means[[o]] + eps, numeric(n))
  X <- Z
  X[, 1] <- exp(Z[, 1] / 2)
  X[, 2] <- Z[, 2] / (1 + exp(Z[, 1]))
  X[, 3] <- (Z[, 1] * Z[, 3] / 25 + 0.6)^3
  X[, 4] <- (Z[, 2] + Z[, 4] + 20)^2
  list(Y = Y, W = W, X = X, X_true = Z, tau = 0, tau_i = 1.5 * g, e = e, logit = logit, e0 = 1 - e)
}

#' the observed covariates and their squares
wc_basis <- function(X) cbind(X, X^2)

#' a cell's common lambda grid. 100 log spaced points from the first rep's lambda_max to 1e-4
wc_lambda_grid <- function(X, W) {
  path <- balnet(X, W, target = "control", max.imbalance = wc_max_imbalance)$lambda
  exp(seq(log(max(path)), log(wc_max_imbalance), length.out = 100))
}

#' Algorithm 1 criterion (Wang and Zubizarreta 2020). mean abs std imbalance of weighted controls
#' vs the treated means over B subsamples of a share frac of the controls
wc_algorithm1_criterion <- function(w, basis, W, target, s, B = 10, frac = 0.1) {
  controls <- which(W == 0)
  m <- round(frac * length(controls))
  total <- 0
  for (b in seq_len(B)) {
    i  <- sample(controls, m, replace = TRUE)
    total <- total + mean(abs(colSums(w[i] * basis[i, , drop = FALSE]) / sum(w[i]) - target) / s)
  }
  total / B
}

#' path + five selectors (cv.balnet() w/ 5 folds and cv.boot.balnet()) + Algorithm 1 on one rep
wc_fit <- function(basis, W, tolerances) {
  fit <- function(f, ...) f(basis, W, target = "control", max.imbalance = wc_max_imbalance,
                            maxit = wc_solver$maxit, tol = wc_solver$tol, ...)
  path <- fit(balnet)
  selectors <- list(cv.bloss = fit(cv.balnet, nfolds = 5, type.measure = "balance.loss"),
                    cv.smd   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.mean"),
                    cv.inf   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.inf"),
                    boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
                    boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  on_tolerances <- balweights(path, lambda = tolerances)
  s <- apply(basis, 2, sd)
  criterion <- vapply(seq_along(tolerances), function(k)
    wc_algorithm1_criterion((on_tolerances[, k] - 1) * (1 - W), basis, W, colMeans(basis[W == 1, ]), s),
    numeric(1))
  list(path = path, selectors = selectors, on_tolerances = on_tolerances, algorithm1 = which.min(criterion))
}

#' one rep. ATT of both outcomes from every arm. satt = sample ATT of outcome A
#' Algorithm 1 and the MLE arms are kept as run but not reported in the thesis
estimate_wong_chan <- function(d, basis, lambda_grid, tolerances, clip = TRUE) {
  W <- d$W
  r <- wc_fit(basis, W, tolerances)
  # balnet targets the controls so w - 1 reweights the controls to the treated
  ybar_treated <- colMeans(d$Y[W == 1, , drop = FALSE])
  att_from    <- function(g) drop(ybar_treated - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
  att_weights <- function(w) att_from((w - 1) * (1 - W))
  att_score   <- function(e) att_from(cbind((1 - W) * e / (1 - e)))

  cv_fit <- cv.glmnet(basis, W, family = "binomial", nfolds = 5)
  e_cv   <- drop(predict(cv_fit, basis, s = "lambda.min", type = "response"))
  e_glm  <- suppressWarnings(fitted(glm(W ~ basis, family = binomial)))   # separates at c = 4
  cap    <- if (clip) function(e) pmin(e, 1 - 1e-6) else identity

  out <- list(
    est_path   = att_weights(balweights(r$path, lambda = lambda_grid)),
    lam_end    = min(r$path$lambda),
    est_sel    = cbind(vapply(r$selectors, function(m) att_weights(balweights(m)), numeric(2)),
                       alg1 = att_weights(r$on_tolerances[, r$algorithm1, drop = FALSE])),
    lam_sel    = c(vapply(r$selectors, `[[`, numeric(1), "lambda.min"), alg1 = tolerances[r$algorithm1]),
    est_mle    = cbind(glmnet = att_score(cap(e_cv)), glm = att_score(cap(e_glm)),
                       naive = ybar_treated - colMeans(d$Y[W == 0, , drop = FALSE])),
    lam_glmnet = cv_fit$lambda.min)
  if (clip) out$clipped <- c(glmnet = mean(e_cv[W == 0] > 1 - 1e-6), glm = mean(e_glm[W == 0] > 1 - 1e-6))
  out$satt <- mean(d$tau_i[W == 1])
  out
}

#' every cell of a grid. one .rds each. rep i of any cell uses stream i of the master seed
run_wong_chan_cells <- function(cells, draw, tolerances, replications, out_dir, file_of, clip = TRUE,
                                workers = 1) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  streams <- seed_streams(wc_master_seed, replications)
  one_rep <- function(i, cell, lambda_grid) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    started <- proc.time()[["elapsed"]]
    out <- tryCatch({
      d <- draw(cell)
      estimate_wong_chan(d, wc_basis(d$X), lambda_grid, tolerances, clip = clip)
    }, error = function(e) list(err = conditionMessage(e)))
    out$time_sec <- proc.time()[["elapsed"]] - started
    out
  }
  cluster <- NULL
  if (workers > 1) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(root)
      for (f in c("wong_chan.R", "simulate.R")) source(file.path(root, "replication", "R", f)),
      root = here::here())
    parallel::clusterExport(cluster, c("streams", "draw", "tolerances", "clip", "one_rep"), envir = environment())
  }
  for (j in seq_len(nrow(cells))) {
    cell <- cells[j, , drop = FALSE]
    out_file <- file.path(out_dir, file_of(cell))
    if (file.exists(out_file)) { message("skip ", basename(out_file)); next }
    assign(".Random.seed", streams[[1]], envir = globalenv())
    first <- draw(cell)
    lambda_grid <- wc_lambda_grid(wc_basis(first$X), first$W)
    reps <- if (is.null(cluster)) lapply(seq_len(replications), one_rep, cell = cell, lambda_grid = lambda_grid)
            else parallel::parLapply(cluster, seq_len(replications), one_rep, cell = cell, lambda_grid = lambda_grid)
    saveRDS(list(cell = cell, master_seed = wc_master_seed, n = nrow(first$X), n_rep = replications,
                 lam = lambda_grid, grid = tolerances,
                 balnet = as.character(packageVersion("balnet")), glmnet = as.character(packageVersion("glmnet")),
                 r = R.version.string, date = Sys.Date(), res = reps), out_file)
    message(basename(out_file), ": ", sum(vapply(reps, function(r) !is.null(r$err), logical(1))), " failed")
  }
  invisible(NULL)
}
