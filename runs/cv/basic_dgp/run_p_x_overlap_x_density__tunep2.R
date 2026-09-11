# run_p_x_overlap_x_density__tunep2.R: overlap x density x p {500, 1000, 2000}, lasso, n = 1000, 50 reps
# (smoke-test rep count; the planned 200 was never run).
# From balnet_mini_experiment_edited.R lines 458-521 ("overlap + p regimes").
# Writes results/runs/tunep2_<ii>_r50.rds and figures/tunep2_grid_r50.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n     <- 1000             # these two were globals left over from the
alpha <- 1                # exploration section of the old file

# overlap + p regimes (sparse/dense outcome) ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 50                                  # smoke test, then 200
grid <- expand.grid(overlap = c("good", "moderate", "bad"),
                    dense = c(FALSE, TRUE), p = c(500, 1000, 2000),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   #
  dat <- gen_data(n, p = g$p, overlap = g$overlap, s_y = 5)
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45) # + 45 covs fixed confounding s_y = 5
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
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
  out_file <- file.path(res_dir, sprintf("tunep2_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(200 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}
parallel::stopCluster(cl)

#  plot of p grid ----
files <- list.files(res_dir, sprintf("^tunep2_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s overlap, p = %d, %s",
                                   cell$overlap, cell$p,
                                   if (cell$dense) "dense" else "sparse")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("tunep2_grid_r%d.png", n_rep)), 2400, 1000,
    res = 110)
par(mfcol = c(3, 6), mar = c(4, 4, 3, 1))
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
dev.off()                                        # start 3:02 end 3:13
