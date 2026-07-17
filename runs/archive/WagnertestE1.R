# check_balnet0_floor.R -- diagnostic, not part of the simulation pipeline.
# Tests whether balweights(fit, lambda = 0) gives exact-balance weights
# or clamps to the path floor (candidate mechanism for flat balnet0 bias).

library(balnet)
library(MASS)

# --- DGP copied verbatim from run script (keep in sync manually) ---
generate_data = function(n, p, outcome = "linear", misspec = FALSE, covcor = "iid", overlap = 1) {
  if (covcor == "iid") {
    X <- matrix(rnorm(n * p), n, p)
  } else if (covcor == "ar1") {
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = stats::toeplitz(0.5^seq(0, p - 1)))
  }
  e = 1 / (1 + exp((X[, 1] - 0.5*X[, 2] + 0.25*X[,3] + 0.1*X[,4]) / overlap))
  W <- rbinom(n, 1, e)
  tau <- 0
  if (outcome == "linear") {
    Y <- tau * W + X[,1] + 0.5*X[,2] + 0.5*X[,3] + 0.5*X[,4] + rnorm(n)
  } else if (outcome == "quad1") {
    Y <- tau * W + rowSums(pmax(X[, 1:4], 0)^2) + rnorm(n)
  } else if (outcome == "exp") {
    Y <- tau * W + rowSums(exp(X[, 1:4] / 2)) + rnorm(n)
  }
  if (misspec) {
    Xreal <- X
    X[,1] <- exp(0.5 * Xreal[,1])
    X[,2] <- 10 + Xreal[,2] / (1 + exp(Xreal[,1]))
    X[,3] <- (0.04 * Xreal[,1] * Xreal[,3] + 0.6)^3
    X[,4] <- (Xreal[,2] + Xreal[,4] + 20)^2
    X[,1:4] <- scale(X[,1:4])
  }
  list(Y = Y, W = W, X = X, tau = tau, e = e)
}

n <- 5000; p <- 50
covcor <- "iid"     # matches the E1 file analysed; set "ar1" for the new run

# --- single-dataset checks ---
set.seed(101)
d   <- generate_data(n, p, outcome = "linear", covcor = covcor)
fit <- balnet(d$X, d$W)                 # same call as run_sim's fit_path

# 1. lambda grid floor per arm (floor > 0 => path never reaches 0)
sapply(fit[["_lambda"]], range)

# 2. is lambda = 0 clamped to the floor?
w0 <- balweights(fit, lambda = 0)
wt_floor <- balweights(fit, lambda = min(fit[["_lambda"]]$treated))$treated
wc_floor <- balweights(fit, lambda = min(fit[["_lambda"]]$control))$control
all.equal(w0$treated, wt_floor)
all.equal(w0$control, wc_floor)

# 3. attained imbalance of the lambda = 0 weights vs pooled means
xbar  <- colMeans(d$X)
imb_t <- as.numeric(crossprod(w0$treated, d$X)) / sum(w0$treated) - xbar
imb_c <- as.numeric(crossprod(w0$control, d$X)) / sum(w0$control) - xbar
max(abs(imb_t)); max(abs(imb_c))        # ~1e-8 exact; ~1e-2 floor
max(abs(imb_t) / apply(d$X, 2, sd))     # standardized (lambda units)

# 4. does residual imbalance reproduce the bias?
gam <- c(1, 0.5, 0.5, 0.5, rep(0, p - 4))
sum(gam * (imb_t - imb_c))

# --- check 4 over seeds ---
R <- 20
b <- vapply(seq_len(R), function(r) {
  set.seed(100 + r)
  d   <- generate_data(n, p, outcome = "linear", covcor = covcor)
  fit <- balnet(d$X, d$W)
  w   <- balweights(fit, lambda = 0)
  xbar  <- colMeans(d$X)
  imb_t <- as.numeric(crossprod(w$treated, d$X)) / sum(w$treated) - xbar
  imb_c <- as.numeric(crossprod(w$control, d$X)) / sum(w$control) - xbar
  sum(gam * (imb_t - imb_c))
}, numeric(1))
mean(b); sd(b) / sqrt(R)   # expect ~ -0.010 if the floor explains the bias

# ---confirmation of DGP ---
b_fix <- vapply(seq_len(R), function(r) {
  set.seed(100 + r)
  d   <- generate_data(n, p, outcome = "linear", covcor = covcor)
  fit <- balnet(d$X, d$W, max.imbalance = 1e-6)
  w   <- balweights(fit, lambda = 1e-6)
  xbar  <- colMeans(d$X)
  imb_t <- as.numeric(crossprod(w$treated, d$X)) / sum(w$treated) - xbar
  imb_c <- as.numeric(crossprod(w$control, d$X)) / sum(w$control) - xbar
  sum(gam * (imb_t - imb_c))
}, numeric(1))
mean(b_fix); sd(b_fix) / sqrt(R)   # expect ~0

