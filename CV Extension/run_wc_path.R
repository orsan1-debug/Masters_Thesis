# run_wc_path.R ------------------------------------------------------------
# RMSE along the balnet lambda path on the Wong & Chan design (Wang &
# Zubizarreta 2020, Appendix D.4): ATE, basis X and X^2 (K = 20), outcome
# models A and B sharing draws. Tuning rules as in Erik's tuning
# experiments (cv.balnet: balance.loss, imbalance.mean, imbalance.inf;
# cv.boot.balnet: imbalance.mean, imbalance.inf) plus Wang & Zubizarreta's
# Algorithm 1 run over the same path with their Wong & Chan settings
# (B = 10 subsamples of size n / 10). Same master seed and stream per rep
# as wc_sbw_v2, so draws pair with the sbw run. PATE = 0 for both models;
# SATE stored per rep. Output: results/<batch_id>.rds with lam (common
# evaluation grid), grid (Algorithm 1 grid) and one list per rep.
# Assumes balweights(fit, lambda = <vector>) returns list(treated, control),
# each n x length(lambda) and zero off-arm, as for a scalar lambda.

library(parallel)

dir         <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/CV Extension/"
batch_id    <- "wc_path_v1"
n           <- 5000
n_rep       <- 200
master_seed <- 20260903
grid        <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)  # Algorithm 1

source(paste0(dir, "sbw_wc.R"))          # dgp_wc(), cstat()
library(balnet)

# Seeds (identical streams to wc_sbw_v2) ----
RNGkind("L'Ecuyer-CMRG")
set.seed(master_seed)
seeds <- vector("list", n_rep)
seeds[[1]] <- .Random.seed
for (i in seq_len(n_rep)[-1]) seeds[[i]] <- nextRNGStream(seeds[[i - 1]])

# Common lambda grid for the path curve, from the rep-1 draw ----
assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
d0  <- dgp_wc(n)
lam <- balnet(cbind(d0$X, d0$X^2), d0$W, max.imbalance = 1e-4)$lambda

# One replication ----
one_rep <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  t0  <- proc.time()[["elapsed"]]
  out <- tryCatch({
    d <- dgp_wc(n)
    W <- d$W
    X <- cbind(d$X, d$X^2)
    fit <- \(f, ...) f(X, W, max.imbalance = 1e-4, maxit = 1e4, tol = 1e-5, ...)
    bl  <- fit(balnet)
    sel <- list(
      cv.bloss = fit(cv.balnet, nfolds = 5, type.measure = "balance.loss"),
      cv.smd   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.mean"),
      cv.inf   = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.inf"),
      boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
      boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
    
    # ATE per outcome column: rows A, B; one column per lambda
    tau <- \(w) drop(crossprod(d$Y, w$treated - w$control)) / n
    
    # Algorithm 1 over the path
    wg <- balweights(bl, lambda = grid)
    s  <- apply(X, 2, sd)
    cs <- vapply(seq_along(grid), \(k)
                 cstat(wg$treated[, k] + wg$control[, k], X, W, arms = c(0, 1),
                       target = colMeans(X), s), numeric(1))
    k  <- which.min(cs)
    wk <- list(treated = wg$treated[, k], control = wg$control[, k])
    
    list(est_path = tau(balweights(bl, lambda = lam)),
         lam_end  = min(bl$lambda),
         est_sel  = cbind(vapply(sel, \(m) tau(balweights(m)), numeric(2)),
                          alg1 = tau(wk)),
         lam_sel  = c(vapply(sel, `[[`, numeric(1), "lambda.min"),
                      alg1 = grid[k]),
         sate     = mean(d$tau_i))
  }, error = \(e) list(err = conditionMessage(e)))
  out$time_sec <- proc.time()[["elapsed"]] - t0
  out
}

# Run ----
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("dir", "seeds", "n", "lam", "grid", "one_rep"))
invisible(clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  source(paste0(dir, "sbw_wc.R"))
  library(balnet)
}))
res <- parLapply(cl, seq_len(n_rep), one_rep)
stopCluster(cl)

# Write ----
saveRDS(list(batch_id = batch_id, master_seed = master_seed, n = n,
             n_rep = n_rep, lam = lam, grid = grid,
             balnet = as.character(packageVersion("balnet")),
             r = R.version.string, date = Sys.Date(), res = res),
        paste0(dir, "results/", batch_id, ".rds"))
table(failed = vapply(res, \(r) !is.null(r$err), logical(1)))
summary(vapply(res, `[[`, numeric(1), "time_sec"))