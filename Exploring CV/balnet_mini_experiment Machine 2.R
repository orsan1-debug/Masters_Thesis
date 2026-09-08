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


##### overlap + SNR ----
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
title_of <- function(cell) sprintf("%s overlap, s_y = %d, sigma_y = %d",
                                   cell$overlap, cell$s_y, cell$sigma_y)
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

######  plot of p grid ----
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


#dimensionality does not have an interior solution in sparse or dense regiemes.

# alpha grid (ridge to lasso), fixed confounding ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 50                                  # smoke test, then 50
grid <- expand.grid(overlap = c("good", "moderate", "bad"),
                    dense = c(FALSE, TRUE), alpha = c(0, 0.5, 1),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = g$overlap, s_y = 5)
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45)
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
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
  out_file <- file.path(dir, sprintf("tunea_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(600 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file)
}
parallel::stopCluster(cl)

#  plot of alpha grid ----
files <- list.files(dir, sprintf("^tunea_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s: %s overlap, %s outcome",
                                   switch(as.character(cell$alpha),
                                          "0" = "ridge (a=0)",
                                          "0.5" = "elastic net (a=0.5)",
                                          "1" = "lasso (a=1)"),
                                   cell$overlap,
                                   if (cell$dense) "dense" else "sparse")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("alpha_grid_ridge_to_lasso_r%d.png", n_rep)),
    2400, 1000, res = 110)
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
dev.off()








# ridge + elastic net, overlap + worse , 500 reps, high SNR ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 500
grid <- expand.grid(overlap = c("good", "moderate", "bad", "awful"),
                    dense = c(FALSE, TRUE), alpha = c(0, 0.5),
                    stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = if (g$overlap == "awful") "bad"
                  else g$overlap, s_y = 5)
  if (g$overlap == "awful") {               # c_prop 2.5 -> 4, SD(eta) ~ 6
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(n, 1, plogis(eta))
  }
  if (g$dense) dat$Y <- dat$Y + rowSums(dat$X[, 6:50]) / sqrt(45)
  dat
}
one_rep <- function(rep_i) {                # g, lam exported by run_par
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
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
  out_file <- file.path(dir, sprintf("tunea3_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(800 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot of extended-overlap grid ----
files <- list.files(dir, sprintf("^tunea3_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s: %s overlap, %s outcome",
                                   if (cell$alpha == 0) "ridge (a=0)"
                                   else "elastic net (a=0.5)",
                                   cell$overlap,
                                   if (cell$dense) "dense" else "sparse")
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("alpha_overlap_grid_r%d.png", n_rep)),
    1600, 1300, res = 110)
par(mfcol = c(4, 4), mar = c(4, 4, 3, 1))
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





# ridge + elastic net, lowered path floor, overlap + worse, high SNR ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200                           # 4-decade paths; 500 if time allows
ratio <- 1e-4                          # lambda.min.ratio, default 1e-2
nlam <- 200                            # keeps 50 points per decade
grid <- expand.grid(overlap = c("good", "moderate", "bad", "awful"),
                    alpha = c(0, 0.5), stringsAsFactors = FALSE)
gen_cell <- function(g) {                   # confounders stay at 1/sqrt(5)
  dat <- gen_data(n, p = p, overlap = if (g$overlap == "awful") "bad"
                  else g$overlap, s_y = 5)
  if (g$overlap == "awful") {               # c_prop 2.5 -> 4, SD(eta) ~ 6
    eta <- as.numeric(dat$X[, 1:5] %*% rep(4 / sqrt(5), 5))
    dat$W <- rbinom(n, 1, plogis(eta))
  }
  dat
}
one_rep <- function(rep_i) tryCatch({       # g, lam, ratio, nlam exported
  dat <- gen_cell(g)
  fit <- function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha,
                            lambda.min.ratio = ratio, nlambda = nlam,
                            maxit = 1e4, tol = 1e-5, ...)
  path <- fit(balnet)
  sel <- list(cv.bloss = fit(cv.balnet, type.measure = "balance.loss"),
              cv.smd   = fit(cv.balnet, type.measure = "imbalance.mean"),
              cv.inf   = fit(cv.balnet, type.measure = "imbalance.inf"),
              boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"),
              boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
  list(est_path = colMeans(balweights(path, lambda = lam) * dat$Y),
       lam_end  = min(path$lambda),          # floor actually reached
       n_lam    = length(path$lambda),
       est_sel  = vapply(sel, \(m) mean(balweights(m) * dat$Y), numeric(1)),
       lam_sel  = vapply(sel, \(m) m$lambda.min, numeric(1)))
}, error = \(e) list(err = conditionMessage(e)))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  out_file <- file.path(dir, sprintf("tunea4_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(900 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha,
                lambda.min.ratio = ratio, nlambda = nlam)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot of lowered-floor grid ----
files <- list.files(dir, sprintf("^tunea4_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
title_of <- function(cell) sprintf("%s: %s overlap, floor %g",
                                   if (cell$alpha == 0) "ridge (a=0)"
                                   else "elastic net (a=0.5)",
                                   cell$overlap, ratio)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("alpha_overlap_floor_grid_r%d.png", n_rep)),
    900, 1300, res = 110)
par(mfcol = c(4, 2), mar = c(4, 4, 3, 1))     # rows overlap, cols alpha
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




# dimension p x alpha x overlap, n = 1000, sigma_y = 1, 200 reps ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n_rep <- 200
grid <- expand.grid(overlap = c("moderate", "bad"), alpha = c(0.5, 0.75),
                    p = c(100, 500), stringsAsFactors = FALSE)  # p last = slow last
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
  out_file <- file.path(dir, sprintf("dima_%02d_r%d.rds", i, n_rep))
  if (file.exists(out_file)) next
  g <- grid[i, ]
  set.seed(1900 + i)
  dat <- gen_cell(g)
  lam <- balnet(dat$X, dat$W, target = "treated", alpha = g$alpha)$lambda
  saveRDS(list(cell = g, lam = lam, res = run_par(seq_len(n_rep), one_rep)),
          out_file)
  message(out_file, "  ", format(Sys.time() - t0, digits = 3))
}
parallel::stopCluster(cl)

#  plot: rows overlap, cols alpha within each p block ----
files <- list.files(dir, sprintf("^dima_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("dim_by_alpha_r%d.png", n_rep)), 2000, 800,
    res = 110)
par(mfcol = c(2, 4), mar = c(4, 4, 3, 1))     # file order fills column-wise
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



# high dimension p x alpha x overlap, n = 1000, sigma_y = 1, 200 reps ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
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
  out_file <- file.path(dir, sprintf("dimhi_%02d_r%d.rds", i, n_rep))
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
files <- list.files(dir, sprintf("^dimhi_\\d+_r%d\\.rds$", n_rep),
                    full.names = TRUE)
rmse <- function(m) sqrt(colMeans(m^2))          # true value is 0
cols <- c("blue", "blue", "blue", "red", "red")
ltys <- c(4, 1, 2, 1, 2)
png(file.path(dir, sprintf("dimhi_by_alpha_r%d.png", n_rep)), 3000, 800,
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


# overlap axis at sigma_y = 1, lasso, n = 1000, 500 reps (Mac) ----
dir <- "/Users/otis/Documents/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterEvalQ(cl, library(balnet))
n <- 1000
p <- 100
n_rep <- 500
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
  out_file <- file.path(dir, sprintf("ov_s1_%02d_r%d.rds", i, n_rep))
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


