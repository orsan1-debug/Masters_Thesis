# run_alpha_sweep_bad_overlap__alpha_bad.R: alpha {0, 0.25, 0.75, 1} at bad overlap, sigma_y {1, 5, 10}
# evaluated on shared fits, n = 1000, p = 100, 200 reps (alpha 1 uses maxit 1e5).
# From balnet_mini_experiment_edited.R lines 784-861.
# Writes results/runs/alpha_bad_<ii>_r200.rds and figures/alpha_sweep_bad_r200.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n <- 1000                 # these were globals left over from the
p <- 100                  # exploration section of the old file

# alpha sweep at bad overlap, sigma_y on shared fits, n = 1000, 200 reps ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
sigmas <- c(1, 5, 10)                      # evaluated on the same fits
grid <- data.frame(alpha = c(0, 0.25, 0.75, 1))  # 0.5 = snr_enet; 1 = maxit 1e5
gen_cell <- function() {                    # bad overlap, sparse outcome
  dat <- gen_data(n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
  mu <- as.numeric(dat$X[, 1:5] %*% rep(1 / sqrt(5), 5))
  dat$Y <- mu + outer(dat$Y - mu, sigmas)   # n x 3: same eps, each sigma_y
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam, sigmas exported
  dat <- gen_cell()
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
                            maxit = if (g$alpha == 1) 1e5 else 1e4,
                            tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  w_path <- as.matrix(balweights(path, lambda = lam))
  list(est_path = crossprod(w_path, dat$Y) / n,        # lambda x sigma
       lam_end  = min(path$lambda),
       est_sel  = vapply(sel, \(m) colMeans(as.numeric(balweights(m)) * dat$Y),
                         numeric(length(sigmas))),      # sigma x selector
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(res_dir, sprintf("alpha_bad_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, , drop = FALSE]
  set.seed(1200 + i)
  dat <- gen_cell()
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha,
                maxit = if (g$alpha == 1) 1e5 else 1e4)$lambda
  saveRDS(list(cell = g, lam = lam, sigmas = sigmas,
               res = run_par(seq_len(n_rep), one_rep)), out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot: rows alpha, cols sigma_y ----
files <- list.files(res_dir, sprintf("^alpha_bad_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("alpha_sweep_bad_r%d.png", n_rep)), 1500, 1600,
    res = 110)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1))     # rows alpha, cols sigma_y
for (f in files) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  if (!all(ok)) message(basename(f), ": ", sum(!ok), " failed reps dropped")
  lam_end <- median(vapply(x$res[ok], `[[`, numeric(1), "lam_end"))
  for (s in seq_along(sigmas)) {
    rmse_path <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_path[, s])))
    rmse_sel <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_sel[s, ])))
    plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
         xlab = "lambda (log scale)", ylab = "RMSE",
         main = sprintf("alpha = %g, bad overlap, sigma_y = %g (R2 = %.2f)",
                        x$cell$alpha, sigmas[s], r2(sigmas[s])))
    abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
    abline(v = lam_end, lty = 3, col = "red")     # median floor reached
    abline(h = rmse_sel, col = cols, lty = ltys)
  }
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()
