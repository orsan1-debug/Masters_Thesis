# run_overlap_axis__ov_s1.R: overlap axis, numeric c_prop {1, 2, 3, 5, 6, 8}, sigma_y = 1,
# lasso, maxit 1e4, 500 reps. From balnet_mini_experiment_Machine_2.R
# lines 907-954. Writes results/runs/ov_s1_<ii>_r500.rds (no figure;
# plotted by plot_rmse_grids_for_qmd.R).

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

# overlap axis at sigma_y = 1, lasso, n = 1000, 500 reps (Mac) ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n <- 1000
p <- 100
n_rep <- as.integer(Sys.getenv("N_SIM", "500"))
alpha_run <- 1
grid <- data.frame(c_prop = c(1, 2, 3, 5, 6, 8), sigma_y = 1)
# 0.7, 1.5, 2.5, 4 already in tune4 (good, moderate, bad, awful)
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
  out_file <- file.path(res_dir, sprintf("ov_s1_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(2300 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha_run)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, sprintf("^ov_s1_\\d+_r%d\\.rds$", n_rep), full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "gen_data", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
               n_sim = n_rep, seed = "2300 + i", script = "runs/cv/basic_dgp/run_overlap_axis__ov_s1.R",
               result_file = file.path(res_dir, sprintf("ov_s1_*_r%d.rds", n_rep)))
}
