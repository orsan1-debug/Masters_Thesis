# run_wc_overlap_n.R -------------------------------------------------------
# Overlap x sample-size grid on the Wong & Chan design, ATT: logit scale
# c in {1, 2, 3, 4} x n in {1000, 2500, 5000, 10000, 20000} via
# dgp_wc_overlap() (c = 1, n = 5000 is the paper's design). Tests Wang &
# Zubizarreta claims (1) approximate beats exact, especially under limited
# overlap, (2) their selector lands near the optimum, (3) the tuned
# tolerance shrinks with n. Per cell: balnet path (target = "control"),
# six tuning rules, Algorithm 1 (grid to 0.5), MLE ATT arms (cv.glmnet
# lambda.min and unpenalised logistic, Hajek on controls, propensities
# clipped at 1 - 1e-6 with clipped share stored), raw difference in means.
# Same master seed and stream per rep. Truth is the sample ATT.
# Output: results/wc_on_n<n>_c<c>.rds; finished cells are skipped.

library(parallel)

res_dir     <- Sys.getenv("OUT_DIR", here::here("results", "cv", "wz_replication"))   # OUT_DIR=<tmp> for smoke runs
run_start   <- Sys.time()   # register_run() only if this run writes a cell
n_rep       <- as.integer(Sys.getenv("N_SIM", "1000"))   # N_SIM=2 for a smoke run
master_seed <- 20260903
cells       <- expand.grid(n = c(1000, 2500, 5000, 10000, 20000), c = c(1, 2, 3, 4))
grid        <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)

source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
source(here::here("R", "registry.R"))          # dgp_wc_overlap(), cstat()
library(balnet); library(glmnet)

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
    d <- dgp_wc_overlap(cell$n, overlap = cell$c)
    W <- d$W
    X <- cbind(d$X, d$X^2)
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
clusterExport(cl, c("seeds", "grid", "one_rep"))
invisible(clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
  library(balnet); library(glmnet)
}))
t0 <- Sys.time()
for (j in seq_len(nrow(cells))) {
  cell <- cells[j, ]
  out_file <- file.path(res_dir, paste0("wc_on_n", cell$n, "_c", cell$c, ".rds"))
  if (file.exists(out_file)) next
  assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
  d0  <- dgp_wc_overlap(cell$n, overlap = cell$c)
  l0  <- balnet(cbind(d0$X, d0$X^2), d0$W, target = "control",
                max.imbalance = 1e-4)$lambda
  lam <- exp(seq(log(max(l0)), log(1e-4), length.out = 100))
  clusterExport(cl, c("cell", "lam"))
  res <- parLapply(cl, seq_len(n_rep), one_rep)
  saveRDS(list(cell = cell, master_seed = master_seed, n = cell$n,
               n_rep = n_rep, lam = lam, grid = grid,
               balnet = as.character(packageVersion("balnet")),
               glmnet = as.character(packageVersion("glmnet")),
               r = R.version.string, date = Sys.Date(), res = res), out_file)
  message(out_file, "  failed ",
          sum(vapply(res, \(r) !is.null(r$err), logical(1))),
          "  ", format(Sys.time() - t0, digits = 3))
}
stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, "^wc_on_n\\d+_c\\d+\\.rds$", full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "dgp_wc_overlap", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf", "alg1", "glmnet", "glm", "naive"),
               n_sim = n_rep, seed = master_seed, script = "runs/cv/wz_replication/run_wc_overlap_n.R",
               result_file = file.path(res_dir, "wc_on_n*_c*.rds"))
}