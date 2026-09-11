# run_overlap_axis_maxit1e5__ov_s1_m5.R: overlap axis, numeric c_prop {8, 16}, sigma_y = 1, lasso,
# maxit 1e5, 500 reps. From balnet_mini_experiment_edited.R lines 1097-1145.
# Writes results/runs/ov_s1_m5_<ii>_r500.rds (no figure; plotted by plot_rmse_grids_for_qmd.R).

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
res_dir <- here::here("results", "cv", "basic_dgp")
fig_dir <- here::here("output", "cv", "basic_dgp", "figures")
stopifnot(dir.exists(res_dir), dir.exists(fig_dir))

# overlap axis, maxit 1e5 so the path floor is reached, sigma_y = 1 ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n <- 1000
p <- 100
n_rep <- 500
alpha_run <- 1
maxit_run <- 1e5                            # ov_s1 used 1e4
grid <- data.frame(c_prop = c(8, 16), sigma_y = 1)
gen_cell <- function(g) {                   # numeric c_prop, as the awful cells
  dat <- gen_data(n, p = p, overlap = "good", s_y = 5, sigma_y = g$sigma_y)
  eta <- as.numeric(dat$X[, 1:5] %*% rep(g$c_prop / sqrt(5), 5))
  dat$prop <- plogis(eta)
  dat$W <- rbinom(n, 1, dat$prop)
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha_run,
                            maxit = maxit_run, tol = 1e-5, ...)
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
  out_file <- file.path(res_dir, sprintf("ov_s1_m5_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(2400 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha_run,
                maxit = maxit_run)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)
