# run_wc_path.R ------------------------------------------------------------
# RMSE along the balnet lambda path on the Wong & Chan design (Wang &
# Zubizarreta 2020, Appendix D.4): ATE, basis X and X^2 (K = 20), outcome
# models A and B sharing draws. Tuning rules as in Erik's tuning
# run_wc_att.R -------------------------------------------------------------
# RMSE along the balnet lambda path on the Wong & Chan design (Wang &
# Zubizarreta 2020, Appendix D.4): ATT, basis X and X^2 (K = 20), outcome
# models A and B sharing draws. One-arm fit with target = "control"; ATT
# weights on controls are gamma = w - 1 = e/(1 - e) (balnet Remark 2), so
# the ATT-scale tolerance is (n / n1) lambda. Tuning rules as in Erik's
# tuning experiments (cv.balnet: balance.loss, imbalance.mean,
# imbalance.inf; cv.boot.balnet: imbalance.mean, imbalance.inf) plus Wang &
# Zubizarreta's Algorithm 1 over the same path with their Wong & Chan
# settings (B = 10 subsamples of size n / 10) evaluated on the ATT weights
# against treated means. Same master seed and stream per rep as wc_sbw_v2.
# Truth is the sample ATT (satt), stored per rep; model B has none.

library(parallel)

res_dir     <- Sys.getenv("OUT_DIR", here::here("results", "cv", "wz_replication"))   # OUT_DIR=<tmp> for smoke runs
run_start   <- Sys.time()   # register_run() only if this run writes a cell
batch_id    <- "wc_att_abw_v1"   # new id: keeps wc_att_v1.rds from being overwritten
n           <- 5000
n_rep       <- as.integer(Sys.getenv("N_SIM", "1000"))   # N_SIM=2 for a smoke run
master_seed <- 20260903
grid        <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)  # Algorithm 1

source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
source(here::here("R", "registry.R"))          # dgp_wc(), cstat()
library(balnet)

# Seeds (identical streams to wc_sbw_v2) ----
RNGkind("L'Ecuyer-CMRG")
set.seed(master_seed)
seeds <- vector("list", n_rep)
seeds[[1]] <- .Random.seed
for (i in seq_len(n_rep)[-1]) seeds[[i]] <- nextRNGStream(seeds[[i - 1]])

# Common evaluation grid: log-spaced from the rep-1 lambda_max to 1e-4 ----
assign(".Random.seed", seeds[[1]], envir = .GlobalEnv)
d0  <- dgp_wc(n)
l0  <- balnet(cbind(d0$X, d0$X^2), d0$W, target = "control",
              max.imbalance = 1e-4)$lambda
lam <- exp(seq(log(max(l0)), log(1e-4), length.out = 100))

# One replication ----
one_rep <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  t0  <- proc.time()[["elapsed"]]
  out <- tryCatch({
    d <- dgp_wc(n)
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
    
    # ATT per outcome column from control-arm weights w (n x L, zero on
    # treated): mean(Y | W = 1) minus the gamma-weighted control mean
    ybar1 <- colMeans(d$Y[W == 1, , drop = FALSE])
    att <- \(w) {
      g <- (w - 1) * (1 - W)
      drop(ybar1 - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
    }
    
    # Algorithm 1 over the path, on the ATT weights against treated means
    wg <- balweights(bl, lambda = grid)
    s  <- apply(X, 2, sd)
    cs <- vapply(seq_along(grid), \(k)
                 cstat((wg[, k] - 1) * (1 - W), X, W, arms = 0,
                       target = colMeans(X[W == 1, ]), s), numeric(1))
    k  <- which.min(cs)
    
    # control-arm lasso OR on the same basis X, after every other RNG use so
    # the rows above keep their draws
    ctrl <- W == 0
    m0 <- vapply(seq_len(ncol(d$Y)), \(j) {
      f <- cv.glmnet(X[ctrl, ], d$Y[ctrl, j], family = "gaussian", nfolds = 5)
      predict(f, newx = X, s = "lambda.min")[, 1]
    }, numeric(nrow(X)))
    # ABW ATT for each column of control-arm weights w
    abw <- \(w) drop(apply(as.matrix(w), 2, \(wk)
                           att_aug(d$Y, W, (wk - 1) * (1 - W), m0)))
    
    list(est_path = att(balweights(bl, lambda = lam)),
         abw_path = abw(balweights(bl, lambda = lam)),
         lam_end  = min(bl$lambda),
         est_sel  = cbind(vapply(sel, \(m) att(balweights(m)), numeric(2)),
                          alg1 = att(wg[, k, drop = FALSE])),
         abw_sel  = cbind(vapply(sel, \(m) abw(balweights(m)), numeric(2)),
                          alg1 = abw(wg[, k, drop = FALSE])),
         lam_sel  = c(vapply(sel, `[[`, numeric(1), "lambda.min"),
                      alg1 = grid[k]),
         satt     = mean(d$tau_i[W == 1]))
  }, error = \(e) list(err = conditionMessage(e)))
  out$time_sec <- proc.time()[["elapsed"]] - t0
  out
}

# Run ----
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("seeds", "n", "lam", "grid", "one_rep"))
invisible(clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  source(here::here("R", "dgp.R")); source(here::here("R", "estimators_cv.R"))
  source(here::here("R", "estimators_ipw.R"))   # att_aug()
  library(balnet)
  library(glmnet)
}))
res <- parLapply(cl, seq_len(n_rep), one_rep)
stopCluster(cl)

# Write ----
saveRDS(list(batch_id = batch_id, master_seed = master_seed, n = n,
             n_rep = n_rep, lam = lam, grid = grid,
             balnet = as.character(packageVersion("balnet")),
             r = R.version.string, date = Sys.Date(), res = res),
        file.path(res_dir, paste0(batch_id, ".rds")))
table(failed = vapply(res, \(r) !is.null(r$err), logical(1)))

# *** Registry ***
register_run(component = "cv", dgp = "dgp_wc",
             estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf", "alg1"),
             n_sim = n_rep, seed = master_seed, script = "runs/cv/wz_replication/run_wc_path.R",
             result_file = file.path(res_dir, paste0(batch_id, ".rds")))
summary(vapply(res, `[[`, numeric(1), "time_sec"))