# run_wc_basis.R -----------------------------------------------------------
# Basis-size x overlap grid on the Wong & Chan design, ATT, n = 5000:
# K in {10, 20, 65, 125} x logit scale c in {1, 3} via dgp_wc_overlap().
# K = 10: first moments of X; 20: plus squares (the paper's run); 65: plus
# all pairwise products; 125: plus 30 pure-noise N(0,1) covariates and
# their squares (estimator-side only; the DGP is unchanged). Tests Wang &
# Zubizarreta Theorem 4: balancing more functions than needed costs little
# under approximate balance, while exact balance degrades or fails. Per
# cell as run_wc_overlap_n.R. Same master seed and stream per rep; noise
# covariates are drawn after the DGP within the rep's stream.
# Output: results/wc_basis_K<K>_c<c>.rds; finished cells are skipped.

library(parallel)

res_dir     <- Sys.getenv("OUT_DIR", here::here("results", "cv", "wz_replication"))   # OUT_DIR=<tmp> for smoke runs
run_start   <- Sys.time()   # register_run() only if this run writes a cell
n           <- 5000
n_rep       <- as.integer(Sys.getenv("N_SIM", "1000"))   # N_SIM=2 for a smoke run
master_seed <- 20260903
cells       <- expand.grid(K = c(10, 20, 65, 125), c = c(1, 3))
grid        <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)

source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
source(here::here("R", "registry.R"))          # dgp_wc_overlap(), cstat()
library(balnet); library(glmnet)

# Basis of size K from the 10 observed covariates ----
basis <- function(X, K) {
  B <- X
  if (K >= 20) B <- cbind(B, X^2)
  if (K >= 65) B <- cbind(B, do.call(cbind, combn(10, 2, \(ij) X[, ij[1]] * X[, ij[2]],
                                                  simplify = FALSE)))
  if (K >= 125) { N <- matrix(rnorm(nrow(X) * 30), nrow(X), 30); B <- cbind(B, N, N^2) }
  B
}

# Seeds ----
RNGkind("L'Ecuyer-CMRG")
set.seed(master_seed)
seeds <- vector("list", n_rep)
seeds[[1]] <- .Random.seed
for (i in seq_len(n_rep)[-1]) seeds[[i]] <- nextRNGStream(seeds[[i - 1]])

# One replication (cell and lam exported per cell) ----
one_rep <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  t0  <- proc.time()[["elapsed"]]
  out <- tryCatch({
    d <- dgp_wc_overlap(n, overlap = cell$c)
    W <- d$W
    X <- basis(d$X, cell$K)
    fit <- \(f, ...) f(X, W, target = "control", max.imbalance = 1e-4,
                       maxit = 1e4, tol = 1e-5, ...)
    bl  <- fit(balnet)
    sel <- list(
      cv.bloss = fit(cv.balnet, nfolds = 5, type.measure = "balance.loss"),
      cv.smd   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.mean"),
      cv.inf   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.inf"),
      boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
      boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
    
    ybar1 <- colMeans(d$Y[W == 1, , drop = FALSE])
    att_w <- \(g) drop(ybar1 - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
    att   <- \(w) att_w((w - 1) * (1 - W))
    att_e <- \(e) att_w(cbind((1 - W) * e / (1 - e)))
    
    wg <- balweights(bl, lambda = grid)
    s  <- apply(X, 2, sd)
    cs <- vapply(seq_along(grid), \(k)
                 cstat((wg[, k] - 1) * (1 - W), X, W, arms = 0,
                       target = colMeans(X[W == 1, ]), s), numeric(1))
    k  <- which.min(cs)
    
    cvg  <- cv.glmnet(X, W, family = "binomial", nfolds = 5)
    e_cv <- drop(predict(cvg, X, s = "lambda.min", type = "response"))
    e_ml <- fitted(glm(W ~ X, family = binomial))
    clip <- \(e) pmin(e, 1 - 1e-6)
    
    list(est_path = att(balweights(bl, lambda = lam)),
         lam_end  = min(bl$lambda),
         est_sel  = cbind(vapply(sel, \(m) att(balweights(m)), numeric(2)),
                          alg1 = att(wg[, k, drop = FALSE])),
         lam_sel  = c(vapply(sel, `[[`, numeric(1), "lambda.min"),
                      alg1 = grid[k]),
         est_mle  = cbind(glmnet = att_e(clip(e_cv)), glm = att_e(clip(e_ml)),
                          naive = ybar1 - colMeans(d$Y[W == 0, , drop = FALSE])),
         lam_glmnet = cvg$lambda.min,
         clipped  = c(glmnet = mean(e_cv[W == 0] > 1 - 1e-6),
                      glm = mean(e_ml[W == 0] > 1 - 1e-6)),
         satt     = mean(d$tau_i[W == 1]))
  }, error = \(e) list(err = conditionMessage(e)))
  out$time_sec <- proc.time()[["elapsed"]] - t0
  out
}

# Run ----
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("seeds", "n", "grid", "basis", "one_rep"))
invisible(clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
  library(balnet); library(glmnet)
}))
t0 <- Sys.time()
for (j in seq_len(nrow(cells))) {
  cell <- cells[j, ]
  out_file <- file.path(res_dir, paste0("wc_basis_K", cell$K, "_c", cell$c, ".rds"))
  if (file.exists(out_file)) next
  assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
  d0  <- dgp_wc_overlap(n, overlap = cell$c)
  l0  <- balnet(basis(d0$X, cell$K), d0$W, target = "control",
                max.imbalance = 1e-4)$lambda
  lam <- exp(seq(log(max(l0)), log(1e-4), length.out = 100))
  clusterExport(cl, c("cell", "lam"))
  res <- parLapply(cl, seq_len(n_rep), one_rep)
  saveRDS(list(cell = cell, master_seed = master_seed, n = n, n_rep = n_rep,
               lam = lam, grid = grid,
               balnet = as.character(packageVersion("balnet")),
               glmnet = as.character(packageVersion("glmnet")),
               r = R.version.string, date = Sys.Date(), res = res), out_file)
  message(out_file, "  failed ",
          sum(vapply(res, \(r) !is.null(r$err), logical(1))),
          "  ", format(Sys.time() - t0, digits = 3))
}
stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, "^wc_basis_K\\d+_c\\d+\\.rds$", full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "dgp_wc_overlap", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf", "alg1", "glmnet", "glm", "naive"),
               n_sim = n_rep, seed = master_seed, script = "runs/cv/wz_replication/run_wc_basis.R",
               result_file = file.path(res_dir, "wc_basis_K*_c*.rds"))
}