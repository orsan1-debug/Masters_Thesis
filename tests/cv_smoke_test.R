# cv_smoke_test.R: does the shared machinery run on this machine?
# gen_data() and run_par() from R/dgp_cv.R, balnet / cv.balnet / cv.boot.balnet
# from the cv-experimental branch, through the same one_rep() as
# run_noise_x_overlap_density_diluted__tune.R,
# on two cells of the tune grid with n_rep = 2. Writes to a temporary folder,
# never to results/runs, so nothing here can be mistaken for a real run.
# Takes a few minutes. Passing means: paths resolve on this OS, the cluster
# starts, every selector returns, and result shapes match what the plot and
# summary scripts expect. It says nothing about reproducing old rds files;
# that needs a full-length cell (run_noise_x_overlap_density_diluted__tune.R,
# then all.equal against the old file).

rm(list = ls())
library(balnet)
# *** Setup ***
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))      # gen_data(), run_par(); run_par() needs cl below
stopifnot(exists("gen_data"), exists("run_par"),
          exists("cv.boot.balnet"))   # FALSE = wrong balnet branch installed
out_dir <- tempfile("cv_smoke_")
dir.create(out_dir)

n     <- 1000
p     <- 100
alpha <- 1
n_rep <- 2                # smoke test only; run scripts use 200 to 1000
grid <- expand.grid(overlap = c("good", "bad"), s_y = 5, sigma_y = 1,
                    stringsAsFactors = FALSE)   # two cells of the tune grid
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))

one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- do.call(gen_data, c(list(n = n, p = p), g))
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

# *** Simulation ***
sels <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  set.seed(100 + i)
  dat <- do.call(gen_data, c(list(n = n, p = p), g))
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  res <- run_par(seq_len(n_rep), one_rep)
  saveRDS(list(cell = g, lam = lam, res = res),
          file.path(out_dir, sprintf("smoke%02d_r%d.rds", i, n_rep)))
  # shapes the plot and summary scripts rely on
  stopifnot(length(res) == n_rep,
            all(vapply(res, \(r) length(r$est_path) == length(lam), logical(1))),
            all(vapply(res, \(r) identical(names(r$est_sel), sels), logical(1))),
            all(vapply(res, \(r) all(is.finite(r$est_sel)), logical(1))))
  message(sprintf("cell %d (%s overlap) ok   %s", i, g$overlap,
                  format(Sys.time() - t0, digits = 3)))
}
parallel::stopCluster(cl)
message("smoke test passed; scratch files in ", out_dir)
