# run_bad_overlap_noise10_1000reps__tune4_19.R: bad overlap, sigma_y = 10 rerun with 1000 reps, maxit 1e4;
# same seed as cell 19 of tune4, replaces tune4_19_r500.
# From balnet_mini_experiment_edited.R lines 1196-1258 (its private copies of
# gen_data() and run_par() were identical to R/dgp_cv.R and are not repeated).
# Writes results/runs/tune4_19_r1000.rds.

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

# bad overlap, sigma_y = 10 rerun, 1000 reps, maxit 1e4, replaces tune4_19 ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n <- 1000
p <- 100
n_rep <- as.integer(Sys.getenv("N_SIM", "1000"))
alpha_run <- 1
g <- data.frame(overlap = "bad", sigma_y = 10, stringsAsFactors = FALSE)
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
out_file <- file.path(res_dir, sprintf("tune4_19_r%d.rds", n_rep))
t0 <- Sys.time()
if (!file.exists(out_file)) {
  set.seed(1000 + 19)                       # same seed grid as tune4_19
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha_run)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, sprintf("^tune4_19_r%d\\.rds$", n_rep), full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "gen_data", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
               n_sim = n_rep, seed = "1000 + 19", script = "runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R",
               result_file = file.path(res_dir, sprintf("tune4_19_r%d.rds", n_rep)))
}
