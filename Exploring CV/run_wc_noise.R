# run_wc_noise.R -----------------------------------------------------------
# Noise x overlap grid on the Wong & Chan design, ATT, n = 5000: outcome
# noise SD sigma in {1, 10, 30, 100} x logit scale c in {1, 2, 3, 4} via
# dgp_wc_overlap() (sigma = 1, c = 1 is the paper's design). Tests whether
# an interior RMSE-optimal lambda appears on model A as noise grows
# (Ben-Michael et al. 2021 §7.3; Bruns-Smith et al. 2025 §7.1.1) and
# whether any rule follows it. Per cell as run_wc_overlap_n.R, Algorithm 1
# grid to 1. Same master seed and stream per rep. Truth is the sample ATT.
# Output: results/wc_noise_s<sigma>_c<c>.rds; finished cells are skipped.

library(parallel)

dir         <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/CV Extension/"
n           <- 5000
n_rep       <- 1000
master_seed <- 20260903
cells       <- expand.grid(sigma = c(1, 10, 30, 100), c = c(1, 2, 3, 4))
grid        <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)

source(paste0(dir, "sbw_wc.R"))          # dgp_wc_overlap(), cstat()
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
    d <- dgp_wc_overlap(n, overlap = cell$c, sigma = cell$sigma)
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
clusterExport(cl, c("dir", "seeds", "n", "grid", "one_rep"))
invisible(clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  source(paste0(dir, "sbw_wc.R"))
  library(balnet); library(glmnet)
}))
t0 <- Sys.time()
for (j in seq_len(nrow(cells))) {
  cell <- cells[j, ]
  out_file <- paste0(dir, "results/wc_noise_s", cell$sigma, "_c", cell$c, ".rds")
  if (file.exists(out_file)) next
  assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
  d0  <- dgp_wc_overlap(n, overlap = cell$c, sigma = cell$sigma)
  l0  <- balnet(cbind(d0$X, d0$X^2), d0$W, target = "control",
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