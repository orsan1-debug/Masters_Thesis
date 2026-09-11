# run_noise_x_overlap_density_fixed__tune2.R: overlap x density x sigma_y, fixed confounding (5 confounders at
# 1/sqrt(5) + 45 outcome-only predictors at 1/sqrt(45)), lasso, n = 1000, p = 100,
# 500 reps. From balnet_mini_experiment_edited.R lines 523-590.
# Writes results/runs/tune2_<ii>_r500.rds and figures/tune2_grid_r500.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n     <- 1000             # these three were globals left over from the
p     <- 100              # exploration section of the old file
alpha <- 1

# SNR grid, fixed confounding ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
grid <- expand.grid(overlap = c("good", "moderate", "bad"),
                    dense = c(FALSE, TRUE), sigma_y = c(1, 5),   # EDIT
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # EDIT: confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = g$overlap, s_y = 5, sigma_y = g$sigma_y)
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45)
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)                                             # EDIT
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
  out_file <- file.path(res_dir, sprintf("tune2_%02d_r%d.rds", i, n_rep))  # EDIT
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(500 + i)                                              # EDIT
  dat <- gen_cell(g)                                             # EDIT
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}
parallel::stopCluster(cl)

#  plot of grid  ----
files <- list.files(res_dir, sprintf("^tune2_\\d+_r%d\\.rds$", n_rep),   # EDIT
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s SNR: %s outcome, %s overlap",
                                   if (cell$sigma_y == 1) "high" else "low",
                                   if (cell$dense) "dense" else "sparse",  # EDIT
                                   cell$overlap)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("tune2_grid_r%d.png", n_rep)), 1600, 1000,   # EDIT
    res = 110)
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
