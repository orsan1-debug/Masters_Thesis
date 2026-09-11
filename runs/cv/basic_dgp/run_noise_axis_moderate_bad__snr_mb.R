# run_noise_axis_moderate_bad__snr_mb.R: noise axis at moderate and bad overlap, sigma_y {4, 7}, lasso,
# n = 1000, p = 100, 500 reps.
# From balnet_mini_experiment_edited.R lines 1261-1326 (its private copies of
# gen_data() and run_par() were identical to R/dgp_cv.R and are not repeated).
# Writes results/runs/snr_mb_<ii>_r500.rds.

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))

# noise axis at moderate and bad overlap: sigma_y 4 and 7, 500 reps ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n <- 1000
p <- 100
n_rep <- 500
alpha_run <- 1
grid <- expand.grid(overlap = c("moderate", "bad"), sigma_y = c(4, 7),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) gen_data(n, p = p, overlap = g$overlap, s_y = 5,
                                 sigma_y = g$sigma_y)
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
  out_file <- file.path(res_dir, sprintf("snr_mb_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(2600 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha_run)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)
