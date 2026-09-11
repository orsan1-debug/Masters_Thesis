# run_ridge_enet_x_overlap_x_density__tunea3.R: ridge and elastic net (alpha 0, 0.5) x overlap (good..awful) x
# density (fixed confounding), n = 1000, p = 100, sigma_y = 1, 500 reps.
# From balnet_mini_experiment_Machine_2.R lines 602-678.
# Writes results/runs/tunea3_<ii>_r500.rds and figures/alpha_overlap_grid_r500.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n <- 1000                 # these were globals left over from the
p <- 100                  # exploration section of the old file

# ridge + elastic net, overlap + worse , 500 reps, high SNR ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
grid <- expand.grid(overlap = c("good", "moderate", "bad", "awful"),
                    dense = c(FALSE, TRUE), alpha = c(0, 0.5),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = if (g$overlap == "awful") "bad"
                  else g$overlap, s_y = 5)
  if (g$overlap == "awful") {               # c_prop 2.5 -> 4, SD(eta) ~ 6
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(n, 1, plogis(eta))
  }
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45)
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
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
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(res_dir, sprintf("tunea3_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(800 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot of extended-overlap grid ----
files <- list.files(res_dir, sprintf("^tunea3_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s: %s overlap, %s outcome",
                                   if (cell$alpha == 0) "ridge (a=0)"
                                   else "elastic net (a=0.5)",
                                   cell$overlap,
                                   if (cell$dense) "dense" else "sparse")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("alpha_overlap_grid_r%d.png", n_rep)),
    1600, 1300, res = 110)
par(mfcol = c(4, 4), mar = c(4, 4, 3, 1))
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
