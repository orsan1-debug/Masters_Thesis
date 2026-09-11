# run_noise_x_overlap_enet__snr_enet.R: same grid as tune4, elastic net alpha = 0.5, 500 reps.
# From balnet_mini_experiment_edited.R lines 674-783.
# Writes results/runs/snr_enet_<ii>_r500.rds, figures/snr_by_overlap_enet_r500.png
# and figures/ladder_enet_r500.png.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))
n <- 1000                 # these were globals left over from the
p <- 100                  # exploration section of the old file

# SNR x overlap grid, elastic net, n = 1000, 500 reps ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
alpha_run <- 0.5
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
one_rep <- function(rep_i) tryCatch({       # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha_run,
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
  out_file <- file.path(res_dir, sprintf("snr_enet_%02d_r%d.rds", i, n_rep))
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

#  plot: rows overlap, cols SNR ----
files <- list.files(res_dir, sprintf("^snr_enet_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
title_of <- function(cell) sprintf("%s overlap, sigma_y = %d (R2 = %.2f)",
                                   cell$overlap, cell$sigma_y, r2(cell$sigma_y))
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("snr_by_overlap_enet_r%d.png", n_rep)), 2000, 1600,
    res = 110)
par(mfcol = c(4, 5), mar = c(4, 4, 3, 1))
for (f in files) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  if (!all(ok)) message(basename(f), ": ", sum(!ok), " failed reps dropped")
  bind <- function(k) do.call(rbind, lapply(x$res[ok], `[[`, k))
  rmse_path <- rmse(bind("est_path"))
  rmse_sel <- rmse(bind("est_sel"))
  plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
       xlab = "lambda (log scale)", ylab = "RMSE", main = title_of(x$cell))
  abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
  abline(v = median(bind("lam_end")), lty = 3, col = "red")  # median floor
  abline(h = rmse_sel, col = cols, lty = ltys)
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()

#  ladder: median selected lambda vs sigma_y, against the RMSE-optimal lambda ----
sels <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
lad <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  lam_sel <- do.call(rbind, lapply(x$res[ok], `[[`, "lam_sel"))
  rmse_path <- rmse(do.call(rbind, lapply(x$res[ok], `[[`, "est_path")))
  data.frame(x$cell, selector = c(sels, "optimal"),
             lam = c(apply(lam_sel[, sels], 2, median),
                     x$lam[which.min(rmse_path)]))
}))
png(file.path(fig_dir, sprintf("ladder_enet_r%d.png", n_rep)), 1800, 500,
    res = 110)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
for (ov in c("good", "moderate", "bad", "awful")) {
  d <- lad[lad$overlap == ov, ]
  plot(NA, xlim = range(d$sigma_y), ylim = range(d$lam), log = "xy",
       xlab = "sigma_y (log scale)", ylab = "lambda (log scale)",
       main = paste(ov, "overlap"))
  for (k in seq_along(sels)) {
    dk <- d[d$selector == sels[k], ]
    dk <- dk[order(dk$sigma_y), ]
    lines(dk$sigma_y, dk$lam, col = cols[k], lty = ltys[k], lwd = 2)
  }
  dk <- d[d$selector == "optimal", ]
  dk <- dk[order(dk$sigma_y), ]
  lines(dk$sigma_y, dk$lam, col = "gray40", lwd = 3)
}
legend("topleft", c(sels, "RMSE-optimal"), col = c(cols, "gray40"),
       lty = c(ltys, 1), lwd = 2, bg = "white")
dev.off()
