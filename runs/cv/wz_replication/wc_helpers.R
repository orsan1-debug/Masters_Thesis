# shared steps of the Wong & Chan run files. needs R/dgp.R estimators_cv.R registry.R path_floor.R
# rep i always draws from seeds[[i]] so its numbers dont depend on cell file or worker count

#' master seed of every Wong & Chan run
wc_master_seed <- 20260903

#' one L'Ecuyer-CMRG stream per rep. 2^127 draws apart and the first from set.seed(master_seed)
wc_seeds <- function(n_rep, master_seed = wc_master_seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(master_seed)
  seeds <- vector("list", n_rep)
  seeds[[1]] <- .Random.seed
  for (i in seq_len(n_rep)[-1])
    seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
  seeds
}

#' the paper's basis X and X^2 (K = 20). cell unused
wc_basis_sq <- function(X, cell) cbind(X, X^2)

#' common lambda grid of a cell. 100 log spaced points from rep 1's lambda_max to the floor
wc_lam_grid <- function(X, W) {
  l0 <- balnet(X, W, target = "control", max.imbalance = wc_floor)$lambda
  exp(seq(log(max(l0)), log(wc_floor), length.out = 100))
}

#' path plus the 5 tuning rules and Wang & Zubizarreta Algorithm 1 for one rep
#' alg 1 uses their cstat() settings (B = 10 subsamples of n / 10). rng order cv folds boot alg 1
wc_fit_rules <- function(X, W, grid) {
  fit <- function(f, ...) f(X, W, target = "control", max.imbalance = wc_floor,
                            maxit = wc_solver$maxit, tol = wc_solver$tol, ...)
  bl  <- fit(balnet)
  sel <- list(
    cv.bloss = fit(cv.balnet, nfolds = 5, type.measure = "balance.loss"),
    cv.smd   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.mean"),
    cv.inf   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.inf"),
    boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
    boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  wg <- balweights(bl, lambda = grid)
  s  <- apply(X, 2, sd)
  cs <- vapply(seq_along(grid), function(k)
    cstat((wg[, k] - 1) * (1 - W), X, W, arms = 0,
          target = colMeans(X[W == 1, ]), s), numeric(1))
  list(bl = bl, sel = sel, wg = wg, k = which.min(cs))
}

#' ATT per outcome col from control weights g. treated mean minus g weighted control mean (Hajek)
wc_att_g <- function(Y, W) {
  ybar1 <- colMeans(Y[W == 1, , drop = FALSE])
  function(g) drop(ybar1 - sweep(crossprod(Y, g), 2, colSums(g), "/"))
}

#' one rep. balnet ATT weights gamma = w - 1 = e / (1 - e) (balnet Remark 2) plus MLE arms
#' MLE fits after all balnet rng draws. clip caps e at 1 - 1e-6
wc_estimate <- function(d, X, lam, grid, clip = TRUE) {
  W     <- d$W
  r     <- wc_fit_rules(X, W, grid)
  att_g <- wc_att_g(d$Y, W)
  att   <- function(w) att_g((w - 1) * (1 - W))
  att_e <- function(e) att_g(cbind((1 - W) * e / (1 - e)))

  cvg  <- cv.glmnet(X, W, family = "binomial", nfolds = 5)
  e_cv <- drop(predict(cvg, X, s = "lambda.min", type = "response"))
  e_ml <- fitted(glm(W ~ X, family = binomial))
  cap  <- if (clip) function(e) pmin(e, 1 - 1e-6) else identity

  out <- list(
    est_path   = att(balweights(r$bl, lambda = lam)),
    lam_end    = min(r$bl$lambda),
    est_sel    = cbind(vapply(r$sel, function(m) att(balweights(m)), numeric(2)),
                       alg1 = att(r$wg[, r$k, drop = FALSE])),
    lam_sel    = c(vapply(r$sel, `[[`, numeric(1), "lambda.min"),
                   alg1 = grid[r$k]),
    est_mle    = cbind(glmnet = att_e(cap(e_cv)), glm = att_e(cap(e_ml)),
                       naive = colMeans(d$Y[W == 1, , drop = FALSE]) -
                         colMeans(d$Y[W == 0, , drop = FALSE])),
    lam_glmnet = cvg$lambda.min)
  if (clip)
    out$clipped <- c(glmnet = mean(e_cv[W == 0] > 1 - 1e-6),
                     glm    = mean(e_ml[W == 0] > 1 - 1e-6))
  out$satt <- mean(d$tau_i[W == 1])
  out
}

#' one rep on its own stream seeds[[i]]. an error is kept as list(err = message)
wc_rep <- function(seeds, draw, basis, lam, grid, cell, estimate = wc_estimate,
                   ...) {
  args <- list(...)
  force(seeds); force(draw); force(basis); force(lam); force(grid)
  force(cell); force(estimate)
  function(i) {
    assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    t0  <- proc.time()[["elapsed"]]
    out <- tryCatch({
      d <- draw(cell)
      do.call(estimate, c(list(d, basis(d$X, cell), lam, grid), args))
    }, error = function(e) list(err = conditionMessage(e)))
    out$time_sec <- proc.time()[["elapsed"]] - t0
    out
  }
}

#' load a worker with the code a rep needs
wc_worker_setup <- function(files) {
  RNGkind("L'Ecuyer-CMRG")
  for (f in files) source(here::here(f))
  library(balnet)
  library(glmnet)
  NULL
}

#' PSOCK workers one fewer than the cores
wc_cluster <- function(also = character()) {
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterCall(cl, wc_worker_setup,
                        c("R/dgp.R", "R/estimators_cv.R", also,
                          "runs/cv/wz_replication/path_floor.R",
                          "runs/cv/wz_replication/wc_helpers.R"))
  cl
}

#' every cell of a Wong & Chan grid to one rds per cell. skips cells whose rds exists
#' cells is a data frame or a vector of logit scales c. batch goes in registry.csv if a cell ran
run_wc_cells <- function(cells, draw, grid, n_rep, res_dir, file_of, glob,
                         script, basis = wc_basis_sq, clip = TRUE) {
  cell_at <- function(j) if (is.data.frame(cells)) cells[j, ] else cells[[j]]
  todo <- Filter(function(j) !file.exists(file.path(res_dir, file_of(cell_at(j)))),
                 seq_len(NROW(cells)))
  if (!length(todo)) {
    message("every cell of ", glob, " exists: skipped")
    return(invisible(NULL))
  }
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  seeds <- wc_seeds(n_rep)
  cl    <- wc_cluster()
  on.exit(parallel::stopCluster(cl), add = TRUE)

  t0 <- Sys.time()
  for (j in todo) {
    cell     <- cell_at(j)
    out_file <- file.path(res_dir, file_of(cell))
    assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
    d0  <- draw(cell)
    lam <- wc_lam_grid(basis(d0$X, cell), d0$W)
    res <- parallel::parLapply(cl, seq_len(n_rep),
                               wc_rep(seeds, draw, basis, lam, grid, cell,
                                      clip = clip))
    saveRDS(list(cell = cell, master_seed = wc_master_seed,
                 n = as.numeric(nrow(d0$X)), n_rep = n_rep, lam = lam,
                 grid = grid, balnet = as.character(packageVersion("balnet")),
                 glmnet = as.character(packageVersion("glmnet")),
                 r = R.version.string, date = Sys.Date(), res = res), out_file)
    message(out_file, "  failed ",
            sum(vapply(res, function(r) !is.null(r$err), logical(1))),
            "  ", format(Sys.time() - t0, digits = 3))
  }
  register_run(component = "cv", dgp = "dgp_wc_overlap",
               estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd",
                              "boot.inf", "alg1", "glmnet", "glm", "naive"),
               n_sim = n_rep, seed = wc_master_seed, script = script,
               result_file = file.path(res_dir, glob))
  invisible(NULL)
}
