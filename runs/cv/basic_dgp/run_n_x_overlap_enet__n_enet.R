# run_n_x_overlap_enet__n_enet.R: n {250, 500, 2000} x overlap {bad, awful}, elastic net alpha = 0.5,
# sigma_y {1, 3, 10} on shared fits, p = 100, 200 reps.
# From balnet_mini_experiment_edited.R lines 922-1011 (lines 862-921 were a
# byte-identical duplicate of the run block and are not carried over).
# Writes results/runs/n_enet_<ii>_r200.rds and figures/n_by_overlap_enet_r200.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
source(here::here("R", "registry.R"))
res_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "basic_dgp"))   # OUT_DIR=<tmp> for smoke runs
run_start <- Sys.time()   # register_run() only if this run writes a cell
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
p <- 100                  # was a global left over from the exploration
                          # section of the old file

# n sweep under elastic net, bad/awful overlap, sigma_y on shared fits ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- as.integer(Sys.getenv("N_SIM", "200"))
sigmas <- c(1, 3, 10)                      # evaluated on the same fits
grid <- expand.grid(overlap = c("bad", "awful"), n = c(250, 500, 2000),
                    stringsAsFactors = FALSE)      # add 5000 for overnight
gen_cell <- function(g) {                   # sparse outcome; awful = c_prop 4
  dat <- gen_data(g$n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
  if (g$overlap == "awful") {
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(g$n, 1, plogis(eta))
  }
  mu <- as.numeric(dat$X[, 1:5] %*% rep(1 / sqrt(5), 5))
  dat$Y <- mu + outer(dat$Y - mu, sigmas)   # n x 3: same eps, each sigma_y
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam, sigmas exported
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = 0.5,
                            maxit = 1e4, tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  w_path <- as.matrix(balweights(path, lambda = lam))
  list(est_path = crossprod(w_path, dat$Y) / g$n,      # lambda x sigma
       lam_end  = min(path$lambda),
       est_sel  = vapply(sel, \(m) colMeans(as.numeric(balweights(m)) * dat$Y),
                         numeric(length(sigmas))),      # sigma x selector
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(res_dir, sprintf("n_enet_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1400 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = 0.5)$lambda
  saveRDS(list(cell = g, lam = lam, sigmas = sigmas,
               res = run_par(seq_len(n_rep), one_rep)), out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, sprintf("^n_enet_\\d+_r%d\\.rds$", n_rep), full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "gen_data", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
               n_sim = n_rep, seed = "1400 + i", script = "runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R",
               result_file = file.path(res_dir, sprintf("n_enet_*_r%d.rds", n_rep)))
}

#  plot: one figure, rows sigma_y x overlap, cols n; path shown only where
#  at least 95% of replicates reached that lambda ----
alpha_run <- 0.5
files <- list.files(res_dir, sprintf("^n_enet_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
runs <- lapply(files, readRDS)
cells <- do.call(rbind, lapply(runs, `[[`, "cell"))
runs <- runs[order(cells$overlap == "awful", cells$n)]  # bad then awful, n up
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("n_by_overlap_enet_r%d.png", n_rep)), 1500, 2400,
    res = 110)
par(mfrow = c(6, 3), mar = c(4, 4, 3, 1))     # rows sigma_y x overlap, cols n
for (s in seq_along(sigmas)) {
  for (x in runs) {
    ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
    if (!all(ok)) message(x$cell$overlap, " n=", x$cell$n, ": ", sum(!ok),
                          " failed reps dropped")
    lam_end <- vapply(x$res[ok], `[[`, numeric(1), "lam_end")
    reached <- colMeans(outer(lam_end, x$lam, "<="))   # share of reps at lambda
    rmse_path <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_path[, s])))
    rmse_path[reached < 0.95] <- NA
    rmse_sel <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_sel[s, ])))
    plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
         xlab = "lambda (log scale)", ylab = "RMSE",
         main = sprintf(
           "alpha = %g: %s overlap, n = %d, sigma_y = %g (R2 = %.2f)",
           alpha_run, x$cell$overlap, x$cell$n, sigmas[s], r2(sigmas[s])))
    abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
    abline(v = median(lam_end), lty = 3, col = "red")  # median floor reached
    abline(h = rmse_sel, col = cols, lty = ltys)
  }
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()
