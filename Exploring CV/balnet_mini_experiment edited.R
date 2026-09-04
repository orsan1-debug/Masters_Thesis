# Some tuning experiments using the balnet experimental branch (has "cv.boot.balnet")
# Estimating E[Y(1)]

rm(list = ls())
set.seed(42)

library(balnet)
#' added in for parallel runs on windows
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
run_par <- function(X, FUN, ...) {
  parallel::clusterSetRNGStream(cl)     
  parallel::clusterExport(cl, setdiff(ls(globalenv()), "cl"))
  parallel::parLapply(cl, X, FUN) # X reps per regieme
}

#' High-dimensional DGP for estimating a counterfactual treated mean.
#'
#' Features are AR(1)-correlated with parameter rho. The first `s_prop`
#' features drive the propensity; the first `s_y` features drive the
#' outcome. They overlap on the first min(s_prop, s_y) coordinates (the
#' "confounders"); any remaining propensity-only coordinates are
#' instruments and outcome-only coordinates are pure predictors.
#'
#' Target: mu1 = E[Y(1)]. By construction E[X] = 0, so mu1 = 0.
#' Only Y for treated units is used by mu1 estimators; Y for controls
#' is returned for convenience but plays no role in the target.
#'
#' Overlap is controlled via the propensity coefficient magnitude:
#'   Var(eta) grows quadratically in the scale; large Var(eta) drives
#'   propensities toward 0 and 1.
gen_data <- function(n = 1000, p = 100, rho = 0.5,
                     s_prop = 5, s_y = 5,
                     overlap = c("bad", "moderate", "good"),
                     sigma_y = 1,
                     seed = NULL) {
  overlap <- match.arg(overlap)
  if (!is.null(seed)) set.seed(seed)
  
  # AR(1) covariance across features
  Sigma <- rho ^ abs(outer(1:p, 1:p, `-`))
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
  
  # Propensity coefficient scale controls overlap
  c_prop <- switch(overlap,
                   good     = 0.7,   # SD(eta) ~ 1,  most props in [0.2, 0.8]
                   moderate = 1.5,   # SD(eta) ~ 2,  some near 0/1
                   bad      = 2.5)   # SD(eta) ~ 3.7, many near 0/1
  # bad      = 4)   # SD(eta) ~ 3.7, many near 0/1
  
  beta_prop <- rep(0, p)
  beta_prop[seq_len(s_prop)] <- c_prop / sqrt(s_prop)
  
  beta_y <- rep(0, p)
  beta_y[seq_len(s_y)] <- 1 / sqrt(s_y)
  
  eta  <- as.numeric(X %*% beta_prop)
  prop <- plogis(eta)
  W    <- rbinom(n, 1, prop)
  
  # Y is Y(1): potential outcome under treatment. mu1 = E[Y(1)] = 0.
  Y <- as.numeric(X %*% beta_y) + rnorm(n, sd = sigma_y)
  
  list(X = X, W = W, Y = Y, prop = prop, true = 0,
       active_prop = which(beta_prop != 0),
       active_y    = which(beta_y != 0))
}



n = 1000
true = 0 # always zero

## sparse lasso with lowest imbalance = best
alpha = 1
sy = 5
p = 100
sigma = 1

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- run_par(seq_len(10), function(rep.i) {  # EDIT: was parallel::mclapply
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W
  
  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)
  
  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)
  
  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.mean",
                     maxit = 1e4, tol = 1e-5)
  
  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.inf",
                     maxit = 1e4, tol = 1e-5)
  
  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)
  
  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)
  
  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,
    
    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,
    
    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,
    
    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,
    
    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,
    
    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
  # }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)



## dense lasso with interior min
alpha = 1
sy = 100
p = 100
sigma = 1

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- run_par(seq_len(25), function(rep.i) {  # EDIT: was parallel::mclapply
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W
  
  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)
  
  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)
  
  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.mean",
                     maxit = 1e4, tol = 1e-5)
  
  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.inf",
                     maxit = 1e4, tol = 1e-5)
  
  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)
  
  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)
  
  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,
    
    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,
    
    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,
    
    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,
    
    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,
    
    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
  # }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)



## RMSE path "reversal" depending on outcome noise
## dense lasso with interior min
alpha = 1
sy = 100
p = 100
sigma = 5

dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)

m = balnet(dat$X, dat$W, target = "treated", alpha = alpha)
lam = m$lambda

# res = replicate(10, {
res <- run_par(seq_len(25), function(rep.i) {  # EDIT: was parallel::mclapply
  dat = gen_data(n, overlap = "bad", p = p, s_y = sy, sigma_y = sigma)
  X = dat$X
  Y = dat$Y
  W = dat$W
  
  bl = balnet(X, W, target = "treated", alpha = alpha, maxit = 1e4, tol = 1e-5)
  
  cv.bloss = cv.balnet(X, W, target = "treated", alpha = alpha,
                       type.measure = "balance.loss",
                       maxit = 1e4, tol = 1e-5)
  
  cv.smd = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.mean",
                     maxit = 1e4, tol = 1e-5)
  
  cv.inf = cv.balnet(X, W, target = "treated", alpha = alpha,
                     type.measure = "imbalance.inf",
                     maxit = 1e4, tol = 1e-5)
  
  boot.smd = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.mean",
                            maxit = 1e4, tol = 1e-5)
  
  boot.inf = cv.boot.balnet(X, W, target = "treated", alpha = alpha,
                            type.measure = "imbalance.inf",
                            maxit = 1e4, tol = 1e-5)
  
  list(
    est = colMeans(balweights(bl, lambda = lam) * Y),
    lam = bl$lambda,
    
    est.cv.bloss = mean(balweights(cv.bloss) * Y),
    lam.cv.bloss = cv.bloss$lambda.min,
    
    est.cv.smd = mean(balweights(cv.smd) * Y),
    lam.cv.smd = cv.smd$lambda.min,
    
    est.cv.inf = mean(balweights(cv.inf) * Y),
    lam.cv.inf = cv.inf$lambda.min,
    
    est.boot.smd = mean(balweights(boot.smd) * Y),
    lam.boot.smd = boot.smd$lambda.min,
    
    est.boot.inf = mean(balweights(boot.inf) * Y),
    lam.boot.inf = boot.inf$lambda.min
  )
  # }, simplify = FALSE)
}, mc.cores = 5, mc.set.seed = TRUE, mc.preschedule = FALSE)

est.mat = do.call(rbind, lapply(res, `[[`, "est"))
rmse = apply(est.mat, 2, function(col) sqrt(mean((col - true)^2)))

rmse.cv.bloss = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.bloss")) - true)^2))
lam.cv.bloss = median(do.call(rbind, lapply(res, `[[`, "lam.cv.bloss")))
rmse.cv.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.smd")) - true)^2))
lam.cv.smd = median(do.call(rbind, lapply(res, `[[`, "lam.cv.smd")))
rmse.cv.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.cv.inf")) - true)^2))
lam.cv.inf = median(do.call(rbind, lapply(res, `[[`, "lam.cv.inf")))

rmse.boot.smd = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.smd")) - true)^2))
lam.boot.smd = median(do.call(rbind, lapply(res, `[[`, "lam.boot.smd")))
rmse.boot.inf = sqrt(mean((do.call(rbind, lapply(res, `[[`, "est.boot.inf")) - true)^2))
lam.boot.inf = median(do.call(rbind, lapply(res, `[[`, "lam.boot.inf")))

plot(
  lam,
  rmse,
  xlim = rev(range(lam)),
  type = 'l',
  xlab = "Log(lambda)",
  log = 'x'
)
abline(v = lam[which.min(rmse)], lty = 3, col = "gray40")
# abline(v = lam.cv.bloss, col = 'blue', lty = 4)
# abline(v = lam.cv.smd, col = 'blue')
# abline(v = lam.cv.inf, col = 'blue', lty = 2)
# abline(v = lam.boot.smd, col = 'red')
# abline(v = lam.boot.inf, col = 'red', lty = 2)

abline(h = rmse.cv.bloss, col = 'blue', lty = 4)
abline(h = rmse.cv.smd, col = 'blue')
abline(h = rmse.cv.inf, col = 'blue', lty = 2)
abline(h = rmse.boot.smd, col = 'red')
abline(h = rmse.boot.inf, col = 'red', lty = 2)
legend("topright",
       legend = c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf"),
       col = c("blue", "blue", "blue", "red", "red"),
       lty = c(4, 1, 2, 1, 2), lwd = 2)






parallel::stopCluster(cl)  # close the cluster


# overlap + SNR ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500                                
grid <- expand.grid(overlap = c("good", "moderate", "bad"), s_y = c(5, 100),
                    sigma_y = c(1, 5), stringsAsFactors = FALSE)

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

for (i in seq_len(nrow(grid))) {
  out_file <- file.path("C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis",
                        "Exploring CV", sprintf("tune%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(100 + i)
  dat <- do.call(gen_data, c(list(n = n, p = p), g))
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}      #started at 1:39AM
parallel::stopCluster(cl)

#  plot of grid  ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
files <- list.files(dir, sprintf("^tune\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s SNR: %s outcome, %s overlap",
                                   if (cell$sigma_y == 1) "high" else "low",
                                   if (cell$s_y == 5) "sparse" else "dense",
                                   cell$overlap)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)

png(file.path(dir, sprintf("tune_grid_r%d.png", n_rep)), 1600, 1000, res = 110)
par(mfcol = c(3, 4), mar = c(4, 4, 3, 1))
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



# overlap + p regimes (sparse/dense outcome) ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"     # Mac path
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 50                                  # smoke test, then 200
grid <- expand.grid(overlap = c("good", "moderate", "bad"),
                    dense = c(FALSE, TRUE), p = c(500, 1000, 2000),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   #
  dat <- gen_data(n, p = g$p, overlap = g$overlap, s_y = 5)
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45) # + 45 covs fixed confounding s_y = 5
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
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
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(dir, sprintf("tunep2_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(200 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}
parallel::stopCluster(cl)

#  plot of p grid ----
files <- list.files(dir, sprintf("^tunep2_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s overlap, p = %d, %s",
                                   cell$overlap, cell$p,
                                   if (cell$dense) "dense" else "sparse")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("tunep2_grid_r%d.png", n_rep)), 2400, 1000,
    res = 110)
par(mfcol = c(3, 6), mar = c(4, 4, 3, 1))
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
dev.off()                                        # start 3:02 end 3:13

# SNR grid, fixed confounding ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"  # EDIT: set once
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
grid <- expand.grid(overlap = c("good", "moderate", "bad"),
                    dense = c(FALSE, TRUE), sigma_y = c(1, 5),   # EDIT
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # EDIT: confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = g$overlap, s_y = 5, sigma_y = g$sigma_y)
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45)
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)                                             # EDIT
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
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(dir, sprintf("tune2_%02d_r%d.rds", i, n_rep))  # EDIT
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(500 + i)                                              # EDIT
  dat <- gen_cell(g)                                             # EDIT
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}
parallel::stopCluster(cl)

#  plot of grid  ----
files <- list.files(dir, sprintf("^tune2_\\d+_r%d\\.rds$", n_rep),   # EDIT
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s SNR: %s outcome, %s overlap",
                                   if (cell$sigma_y == 1) "high" else "low",
                                   if (cell$dense) "dense" else "sparse",  # EDIT
                                   cell$overlap)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("tune2_grid_r%d.png", n_rep)), 1600, 1000,   # EDIT
    res = 110)
par(mfcol = c(3, 4), mar = c(4, 4, 3, 1))
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



# degrading SNR x overlap to 4, sparse outcome ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV" 
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
  out_file <- file.path(dir, sprintf("tune4_%02d_r%d.rds", i, n_rep))
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
files <- list.files(dir, sprintf("^tune4_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s overlap, sigma_y = %d",
                                   cell$overlap, cell$sigma_y)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("SNR_by_overlap_lasso_r%d.png", n_rep)), 2000, 1600,
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

# SNR x overlap grid, elastic net, n = 1000, 500 reps ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
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
  out_file <- file.path(dir, sprintf("snr_enet_%02d_r%d.rds", i, n_rep))
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
files <- list.files(dir, sprintf("^snr_enet_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
title_of <- function(cell) sprintf("%s overlap, sigma_y = %d (R2 = %.2f)",
                                   cell$overlap, cell$sigma_y, r2(cell$sigma_y))
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("snr_by_overlap_enet_r%d.png", n_rep)), 2000, 1600,
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
png(file.path(dir, sprintf("ladder_enet_r%d.png", n_rep)), 1800, 500,
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




# alpha sweep at bad overlap, sigma_y on shared fits, n = 1000, 200 reps ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
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
  out_file <- file.path(dir, sprintf("alpha_bad_%02d_r%d.rds", i, n_rep))
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
files <- list.files(dir, sprintf("^alpha_bad_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("alpha_sweep_bad_r%d.png", n_rep)), 1500, 1600,
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


# n sweep under elastic net, bad/awful overlap, sigma_y on shared fits ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
sigmas <- c(1, 3, 10)                      # evaluated on the same fits
grid <- expand.grid(overlap = c("bad", "awful"), n = c(250, 500, 2000),
                    stringsAsFactors = FALSE)      # add 5000 for overnight
gen_cell <- function(g) {                   # sparse outcome; awful = c_prop 4
  dat <- gen_data(g$n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
  if (g$overlap == "awful") {
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(g$n, 1, plogis(eta))
  }
  mu <- as.numeric(dat$X[, 1:5] %*% rep(1 / sqrt(5), 5))
  dat$Y <- mu + outer(dat$Y - mu, sigmas)   # n x 3: same eps, each sigma_y
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam, sigmas exported
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = 0.5,
                            maxit = 1e4, tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  w_path <- as.matrix(balweights(path, lambda = lam))
  list(est_path = crossprod(w_path, dat$Y) / g$n,      # lambda x sigma
       lam_end  = min(path$lambda),
       est_sel  = vapply(sel, \(m) colMeans(as.numeric(balweights(m)) * dat$Y),
                         numeric(length(sigmas))),      # sigma x selector
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(dir, sprintf("n_enet_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1400 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = 0.5)$lambda
  saveRDS(list(cell = g, lam = lam, sigmas = sigmas,
               res = run_par(seq_len(n_rep), one_rep)), out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)











# n sweep under elastic net, bad/awful overlap, sigma_y on shared fits ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
sigmas <- c(1, 3, 10)                      # evaluated on the same fits
grid <- expand.grid(overlap = c("bad", "awful"), n = c(250, 500, 2000),
                    stringsAsFactors = FALSE)      # add 5000 for overnight
gen_cell <- function(g) {                   # sparse outcome; awful = c_prop 4
  dat <- gen_data(g$n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
  if (g$overlap == "awful") {
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(g$n, 1, plogis(eta))
  }
  mu <- as.numeric(dat$X[, 1:5] %*% rep(1 / sqrt(5), 5))
  dat$Y <- mu + outer(dat$Y - mu, sigmas)   # n x 3: same eps, each sigma_y
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam, sigmas exported
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = 0.5,
                            maxit = 1e4, tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  w_path <- as.matrix(balweights(path, lambda = lam))
  list(est_path = crossprod(w_path, dat$Y) / g$n,      # lambda x sigma
       lam_end  = min(path$lambda),
       est_sel  = vapply(sel, \(m) colMeans(as.numeric(balweights(m)) * dat$Y),
                         numeric(length(sigmas))),      # sigma x selector
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(dir, sprintf("n_enet_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1400 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = 0.5)$lambda
  saveRDS(list(cell = g, lam = lam, sigmas = sigmas,
               res = run_par(seq_len(n_rep), one_rep)), out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot: one figure, rows sigma_y x overlap, cols n; path shown only where
#  at least 95% of replicates reached that lambda ----
alpha_run <- 0.5
files <- list.files(dir, sprintf("^n_enet_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
runs <- lapply(files, readRDS)
cells <- do.call(rbind, lapply(runs, `[[`, "cell"))
runs <- runs[order(cells$overlap == "awful", cells$n)]  # bad then awful, n up
r2 <- function(sigma_y) 2.23 / (2.23 + sigma_y^2)   # signal var beta'Sigma beta
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("n_by_overlap_enet_r%d.png", n_rep)), 1500, 2400,
    res = 110)
par(mfrow = c(6, 3), mar = c(4, 4, 3, 1))     # rows sigma_y x overlap, cols n
for (s in seq_along(sigmas)) {
  for (x in runs) {
    ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
    if (!all(ok)) message(x$cell$overlap, " n=", x$cell$n, ": ", sum(!ok),
                          " failed reps dropped")
    lam_end <- vapply(x$res[ok], `[[`, numeric(1), "lam_end")
    reached <- colMeans(outer(lam_end, x$lam, "<="))   # share of reps at lambda
    rmse_path <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_path[, s])))
    rmse_path[reached < 0.95] <- NA
    rmse_sel <- rmse(do.call(rbind, lapply(x$res[ok], \(r) r$est_sel[s, ])))
    plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
         xlab = "lambda (log scale)", ylab = "RMSE",
         main = sprintf(
           "alpha = %g: %s overlap, n = %d, sigma_y = %g (R2 = %.2f)",
           alpha_run, x$cell$overlap, x$cell$n, sigmas[s], r2(sigmas[s])))
    abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
    abline(v = median(lam_end), lty = 3, col = "red")  # median floor reached
    abline(h = rmse_sel, col = cols, lty = ltys)
  }
}
legend("topright", names(rmse_sel), col = cols, lty = ltys, lwd = 2,
       bg = "white")
dev.off()


# confounder spread s x overlap x penalty, n = 1000, sigma_y = 1, 200 reps ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
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
  out_file <- file.path(dir, sprintf("spread_%02d_r%d.rds", i, n_rep))
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

#  plot: rows overlap, cols s within each penalty (ridge | EN | lasso) ----
files <- list.files(dir, sprintf("^spread_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
pen <- c("0" = "ridge", "0.5" = "elastic net", "1" = "lasso")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("spread_by_penalty_r%d.png", n_rep)), 4500, 1600,
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