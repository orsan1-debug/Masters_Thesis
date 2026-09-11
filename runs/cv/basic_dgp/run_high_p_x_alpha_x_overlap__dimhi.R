# run_high_p_x_alpha_x_overlap__dimhi.R: p {1000, 2000} x alpha {0.5, 0.75, 1} x overlap {moderate, bad},
# n = 1000, sigma_y = 1, 200 reps. From balnet_mini_experiment_Machine_2.R
# lines 835-906. Writes results/runs/dimhi_<ii>_r200.rds and
# figures/dimhi_by_alpha_r200.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n <- 1000                 # was a global left over from the exploration
                          # section of the old file

# high dimension p x alpha x overlap, n = 1000, sigma_y = 1, 200 reps ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
grid <- expand.grid(overlap = c("moderate", "bad"), alpha = c(0.5, 0.75, 1),
                    p = c(1000, 2000), stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # 5 confounders, coef 1/sqrt(5)
  gen_data(n, p = g$p, overlap = g$overlap, s_y = 5, sigma_y = 1)
}
one_rep <- function(rep_i) tryCatch({       # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
                            maxit = 1e4, tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  list(est_path = colMeans(balweights(path, lambda = lam) * dat$Y),
       lam_end  = min(path$lambda),
       est_sel  = vapply(sel, \(m) mean(balweights(m) * dat$Y), numeric(1)),
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(res_dir, sprintf("dimhi_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(2000 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot: rows overlap, cols alpha within each p block ----
files <- list.files(res_dir, sprintf("^dimhi_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("dimhi_by_alpha_r%d.png", n_rep)), 3000, 800,
    res = 110)
par(mfcol = c(2, 6), mar = c(4, 4, 3, 1))     # file order fills column-wise
for (f in files) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  if (!all(ok)) message(basename(f), ": ", sum(!ok), " failed reps dropped")
  lam_end <- vapply(x$res[ok], `[[`, numeric(1), "lam_end")
  reached <- colMeans(outer(lam_end, x$lam, "<="))     # share of reps at lambda
  bind <- function(k) do.call(rbind, lapply(x$res[ok], `[[`, k))
  rmse_path <- rmse(bind("est_path"))
  rmse_path[reached < 0.95] <- NA
  rmse_sel <- rmse(bind("est_sel"))
  plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
       xlab = "lambda (log scale)", ylab = "RMSE",
       main = sprintf("alpha = %g: %s overlap, p = %d, sigma_y = 1",
                      x$cell$alpha, x$cell$overlap, x$cell$p))
  abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
  abline(v = median(lam_end), lty = 3, col = "red")    # median floor reached
  abline(h = rmse_sel, col = cols, lty = ltys)
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()
