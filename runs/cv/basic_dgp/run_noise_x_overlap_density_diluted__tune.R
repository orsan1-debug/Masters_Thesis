# run_noise_x_overlap_density_diluted__tune.R: overlap x s_y x sigma_y grid, lasso, n = 1000, p = 100, 500 reps
# From balnet_mini_experiment_edited.R lines 393-455 ("overlap + SNR").
# Writes results/runs/tune<ii>_r500.rds and figures/tune_grid_r500.png.

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
n     <- 1000             # these three were globals left over from the
p     <- 100              # exploration section of the old file
alpha <- 1

# overlap + SNR ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- as.integer(Sys.getenv("N_SIM", "500"))   
grid <- expand.grid(overlap = c("good", "moderate", "bad"), s_y = c(5, 100),
                    sigma_y = c(1, 5), stringsAsFactors = FALSE)

one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- do.call(gen_data, c(list(n = n, p = p), g))
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha,
                            maxit = 1e4, tol = 1e-5, ...)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  list(est_path = colMeans(balweights(fit(balnet), lambda = lam) * dat$Y),
       est_sel  = vapply(sel, \(m) mean(balweights(m) * dat$Y), numeric(1)),
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}

for (i in seq_len(nrow(grid))) {
  out_file <- file.path(res_dir, sprintf("tune%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(100 + i)
  dat <- do.call(gen_data, c(list(n = n, p = p), g))
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}      #started at 1:39AM
parallel::stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, sprintf("^tune\\d+_r%d\\.rds$", n_rep), full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "gen_data", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
               n_sim = n_rep, seed = "100 + i", script = "runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R",
               result_file = file.path(res_dir, sprintf("tune*_r%d.rds", n_rep)))
}

#  plot of grid  ----
files <- list.files(res_dir, sprintf("^tune\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s SNR: %s outcome, %s overlap",
                                   if (cell$sigma_y == 1) "high" else "low",
                                   if (cell$s_y == 5) "sparse" else "dense",
                                   cell$overlap)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)

png(file.path(fig_dir, sprintf("tune_grid_r%d.png", n_rep)), 1600, 1000, res = 110)
par(mfcol = c(3, 4), mar = c(4, 4, 3, 1))
for (f in files) {
  x <- readRDS(f)
  bind <- function(k) do.call(rbind, lapply(x$res, `[[`, k))
  rmse_path <- rmse(bind("est_path"))
  rmse_sel <- rmse(bind("est_sel"))
  plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
       xlab = "lambda (log scale)", ylab = "RMSE", main = title_of(x$cell))
  abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
  abline(h = rmse_sel, col = cols, lty = ltys)
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()
