# run_confounder_spread_x_penalty__spread.R: confounder spread s {5, 20, 50} x overlap (good..awful) x
# alpha {0, 0.5, 1}, n = 1000, p = 100, sigma_y = 1, 200 reps.
# From balnet_mini_experiment_edited.R lines 1012-1096.
# Writes results/runs/spread_<ii>_r200.rds and figures/spread_by_penalty_r200.png.

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
n <- 1000                 # these were globals left over from the
p <- 100                  # exploration section of the old file

# confounder spread s x overlap x penalty, n = 1000, sigma_y = 1, 200 reps ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- as.integer(Sys.getenv("N_SIM", "200"))
grid <- expand.grid(overlap = c("good", "moderate", "bad", "awful"),
                    s = c(5, 20, 50), alpha = c(0, 0.5, 1),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # s confounders, coef 1/sqrt(s)
  dat <- gen_data(n, p = p, overlap = if (g$overlap == "awful") "bad"
                  else g$overlap, s_prop = g$s, s_y = g$s, sigma_y = 1)
  if (g$overlap == "awful") {               # c_prop 2.5 -> 4
    eta <- as.numeric(dat$X[, seq_len(g$s)] %*% rep(4 / sqrt(g$s), g$s))
    dat$W <- rbinom(n, 1, plogis(eta))
  }
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
                            lambda.min.ratio = if (g$alpha == 0) 1e-4 else 1e-2,
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
  out_file <- file.path(res_dir, sprintf("spread_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1600 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha,
                lambda.min.ratio = if (g$alpha == 0) 1e-4 else 1e-2)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

# *** Registry (only if this run wrote at least one cell) ***
written <- list.files(res_dir, sprintf("^spread_\\d+_r%d\\.rds$", n_rep), full.names = TRUE)
if (any(file.mtime(written) >= run_start)) {
  register_run(component = "cv", dgp = "gen_data", estimators = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
               n_sim = n_rep, seed = "1600 + i", script = "runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R",
               result_file = file.path(res_dir, sprintf("spread_*_r%d.rds", n_rep)))
}

#  plot: rows overlap, cols s within each penalty (ridge | EN | lasso) ----
files <- list.files(res_dir, sprintf("^spread_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
pen <- c("0" = "ridge", "0.5" = "elastic net", "1" = "lasso")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(fig_dir, sprintf("spread_by_penalty_r%d.png", n_rep)), 4500, 1600,
    res = 110)
par(mfcol = c(4, 9), mar = c(4, 4, 3, 1))     # file order fills column-wise
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
       main = sprintf("%s (a=%g): %s overlap, %d confounders",
                      pen[as.character(x$cell$alpha)], x$cell$alpha,
                      x$cell$overlap, x$cell$s))
  abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
  abline(v = median(lam_end), lty = 3, col = "red")    # median floor reached
  abline(h = rmse_sel, col = cols, lty = ltys)
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()
