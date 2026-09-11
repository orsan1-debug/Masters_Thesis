# run_noise_x_overlap_lasso__tune4.R: overlap {good, moderate, bad, awful} x sigma_y {1, 2, 3, 5, 10},
# sparse outcome, lasso, n = 1000, p = 100, 500 reps; awful = c_prop 4.
# From balnet_mini_experiment_edited.R lines 591-673.
# Writes results/runs/tune4_<ii>_r500.rds and figures/SNR_by_overlap_lasso_r500.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n <- 1000                 # these were globals left over from the
p <- 100                  # exploration section of the old file

# degrading SNR x overlap to 4, sparse outcome ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
alpha_run <- 1                              # 1 = lasso; rerun with 0.5 later
grid <- expand.grid(overlap = c("good", "moderate", "bad", "awful"),
                    sigma_y = c(1, 2, 3, 5, 10), stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # sparse outcome; awful = c_prop 4
  dat <- gen_data(n, p = p, overlap = if (g$overlap == "awful") "bad"
                  else g$overlap, s_y = 5, sigma_y = g$sigma_y)
  if (g$overlap == "awful") {
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(n, 1, plogis(eta))
  }
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha_run,
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
  out_file <- file.path(res_dir, sprintf("tune4_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1000 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha_run)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot ----
files <- list.files(res_dir, sprintf("^tune4_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s overlap, sigma_y = %d",
                                   cell$overlap, cell$sigma_y)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("SNR_by_overlap_lasso_r%d.png", n_rep)), 2000, 1600,
    res = 110)
par(mfcol = c(4, 5), mar = c(4, 4, 3, 1))
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

#COMMENT for SNR, both with diluted density and with independent density.

# interior optimum only exists and higher outcome noise
# with low enough SNR, CV can beat endpoint

#with confounding fixed
######## OuverlapGood>Bad: with high SNR overlap deg improves relative endpoint performance, with low SNR it shrinks gap (different CV perform differently, cv bloos and boot inf)
######## outcome density has little/no effect




### ELASTIC NET 
