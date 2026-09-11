# Dropped variants (Phase 3)

Names defined with **different bodies** in more than one file (25 names from
`notes/duplicate_functions.csv`). Rule: keep the `R/` version where one exists; the
script-local variant is recorded here with a unified diff and is **not deleted yet**.
Where no `R/` version exists the variants are script-local helpers (same name,
different job per script); all are kept in place and listed for reference.
Files are given by their current (post-Phase 2/3) paths.

## `arity_ok`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,6 @@
+function(e) {
+    ok <- TRUE
+    walk(e, function(x) if (is.name(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=", "<<-") && length(x) != 3L) 
+        ok <<- FALSE)
+    ok
+}
```

## `att`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/cv/wz_replication/run_wc_path.R`
- variant 2: `runs/cv/wz_replication/run_wc_basis.R`, `runs/cv/wz_replication/run_wc_noise.R`, `runs/cv/wz_replication/run_wc_overlap.R`, `runs/cv/wz_replication/run_wc_overlap_n.R`

Diff of variant 1 vs variant 2:

```diff
@@ -1,4 +1 @@
-function(w) {
-    g <- (w - 1) * (1 - W)
-    drop(ybar1 - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
-}
+function(w) att_w((w - 1) * (1 - W))
```

## `axis_panels`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `analysis/cv/basic_dgp/plot_grids.R`
- variant 2: `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R`

Diff of variant 1 vs variant 2:

```diff
@@ -1,5 +1,5 @@
 function(family, key, keep, prefer = NULL, tag = NULL) {
-    files <- list.files(dir, sprintf("^%s_?\\d+_r\\d+\\.rds$", family), full.names = TRUE)
+    files <- list.files(res_dir, sprintf("^%s_?\\d+_r\\d+\\.rds$", family), full.names = TRUE)
     runs <- lapply(files, readRDS)
     cells <- do.call(rbind, lapply(runs, function(x) {
         cell <- x$cell
```

## `bind`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_machine_2.R`, `archive/exploring_cv/balnet_mini_experiment_windows.R`, `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`, `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`, `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`, `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- variant 2: `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`, `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`, `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`, `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`, `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- variant 3: `analysis/cv/wz_replication/plot_wc_att.R`, `analysis/cv/wz_replication/wc_grids.R`, `analysis/cv/wz_replication/wc_findings.qmd`

Diff of variant 1 vs variant 2 (3 variants in total):

```diff
@@ -1 +1 @@
-function(k) do.call(rbind, lapply(x$res, `[[`, k))
+function(k) do.call(rbind, lapply(x$res[ok], `[[`, k))
```

## `dgp_gen`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/ipw/PoverMDGP2.R`, `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/weak_instruments.R`
- variant 3: `runs/ipw/Sparsity.R`
- variant 4: `runs/ipw/sparsity_2.R`
- variant 5: `runs/ipw/cvfix.R`
- variant 6: `runs/ipw/correctspec_dgp2_n800k.R`
- variant 7: `runs/ipw/SNR.R`
- variant 8: `runs/ipw/treat_prop.R`
- variant 9: `runs/ipw/correct_misspec_dgp2.R`, `runs/ipw/overlap_misspec_dgp2.R`
- variant 10: `runs/ipw/DimensionDGP2.R`, `runs/ipw/DimensionDGP2N5K.R`

Diff of variant 1 vs variant 2 (10 variants in total):

```diff
@@ -0,0 +1 @@
+function(cell) dgp2(n = cell$n, p = 50, s = cell$s, signs = "pos", outcome = c("linear", "quad1"), outcome_set = "fixed4", decay_ps = cell$decay, covcor = "iid", misspec = FALSE, overlap = 1)
```

## `estimate_all`

Canonical (kept): `R/estimators_ipw.R`

Variant in `tests/smoke_simulate.R` (kept in place for now):

```diff
@@ -1,42 +1,6 @@
-function(data, lambdas = c(`0` = 0, `05` = 0.05, `10` = 0.1), nfolds = 5, max_imbalance = 1e-04, cv_curve = FALSE, lambda_grid = NULL) {
-    Y <- as.matrix(data$Y)
-    W <- data$W
-    X <- data$X
-    n <- nrow(X)
-    p <- ncol(X)
-    k <- ncol(Y)
-    e1 <- data$e
-    e0 <- if (is.null(data$e0)) 
-        1 - e1
-    else data$e0
-    fit_bal <- cv.balnet(X, W, nfolds = nfolds, max.imbalance = max_imbalance)
-    w_cv <- balweights(fit_bal)
-    w_fix <- lapply(lambdas, function(l) balweights(fit_bal, lambda = l))
-    w_rate <- balweights(fit_bal, lambda = sqrt(log(p)/n))
-    if (!is.null(lambda_grid)) {
-        w_grid <- lapply(lambda_grid, function(l) balweights(fit_bal, lambda = l))
-        tau_path <- do.call(rbind, lapply(w_grid, ate_bal, Y = Y))
-    }
-    fit_glm <- cv.glmnet(X, W, family = "binomial", nfolds = nfolds)
-    e_hat <- predict(fit_glm, newx = X, s = "lambda.min", type = "response")[, 1]
-    ate_fix <- do.call(rbind, lapply(w_fix, ate_bal, Y = Y))
-    rownames(ate_fix) <- paste0("balnet", names(lambdas))
-    xbar <- colMeans(X)
-    sdx <- apply(X, 2, sd) * sqrt((n - 1)/n)
-    cv1 <- cv_loss(fit_bal, "treated")
-    cv0 <- cv_loss(fit_bal, "control")
-    smd_fix <- numeric(0)
-    for (nm in names(lambdas)) {
-        smd_fix[paste0("smd1_", nm)] <- max_smd(w_fix[[nm]]$treated, X, xbar, sdx)
-        smd_fix[paste0("smd0_", nm)] <- max_smd(w_fix[[nm]]$control, X, xbar, sdx)
-    }
-    diags <- c(lam_balcv1 = fit_bal$lambda.min$treated, lam_balcv0 = fit_bal$lambda.min$control, lam_end1 = min(fit_bal$lambda$treated), lam_end0 = min(fit_bal$lambda$control), lam_glmcv = fit_glm$lambda.min, trunc05 = as.numeric(min(fit_bal$lambda$treated) > 0.05 | min(fit_bal$lambda$control) > 0.05), nnz_balcv1 = sum(coef(fit_bal)$treated[-1, ] != 0), nnz_balcv0 = sum(coef(fit_bal)$control[-1, ] != 0), nnz_glm = sum(coef(fit_glm, s = "lambda.min")[-1] != 0), smd1_cv = max_smd(w_cv$treated, X, xbar, 
-        sdx), smd0_cv = max_smd(w_cv$control, X, xbar, sdx), smd_fix, prev = mean(W), emin = min(e1), emax = max(e1), nout05 = sum(e1 < 0.05 | e1 > 0.95), nout01 = sum(e1 < 0.01 | e1 > 0.99), cvloss_cv1 = cv1[1], cvloss_cv0 = cv0[1], cvloss_end1 = cv1[2], cvloss_end0 = cv0[2])
-    out <- rbind(balnetcv = ate_bal(Y, w_cv), ate_fix, balnetrate = ate_bal(Y, w_rate), glmnetcv_ht = ate_ht(Y, W, e_hat, 1 - e_hat), glmnetcv_hajek = ate_hajek(Y, W, e_hat, 1 - e_hat), oracle_ht = ate_ht(Y, W, e1, e0), oracle_hajek = ate_hajek(Y, W, e1, e0), matrix(diags, length(diags), k, dimnames = list(names(diags), NULL)))
-    colnames(out) <- colnames(Y)
+function(dat, cv_curve = FALSE) {
+    m <- matrix(c(mean(dat$x), median(dat$x), sd(dat$x)), ncol = 1, dimnames = list(c("est_mean", "est_median", "est_sd"), "tau_hat"))
     if (cv_curve) 
-        attr(out, "cv_curve") <- list(lambda = fit_bal$lambda, cv.mean = fit_bal$`_cv.info`$cv.mean)
-    if (!is.null(lambda_grid)) 
-        attr(out, "tau_path") <- list(lambda = lambda_grid, tau = tau_path)
-    out
+        attr(m, "cv_curve") <- dat$x[1:3]
+    m
 }
```

## `f`

More than one `R/` file defines this name with different bodies:

- `R/dgp.R`, `archive/cv_extension/dgp_wc_overlap.R`
- `R/dgp.R`, `archive/tests/dgp2_validity_t1.R`, `tests/dgp2_validity.R`

## `fit`

Canonical (kept): `R/estimators_cv.R`

Variant in `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = 0.5, maxit = 10000, tol = 1e-05, ...)
```

Variant in `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_machine_2.R`, `archive/exploring_cv/balnet_mini_experiment_windows.R`, `tests/cv_smoke_test.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`, `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha, maxit = 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`, `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`, `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`, `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`, `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha_run, maxit = 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`, `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = alpha_run, maxit = maxit_run, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha, lambda.min.ratio = if (g$alpha == 0) 1e-04 else 0.01, maxit = 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha, lambda.min.ratio = ratio, nlambda = nlam, maxit = 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`, `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`, `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`, `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha, maxit = 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(dat$X, dat$W, target = "treated", alpha = g$alpha, maxit = if (g$alpha == 1) 1e+05 else 10000, tol = 1e-05, ...)
```

Variant in `runs/cv/wz_replication/run_wc_basis.R`, `runs/cv/wz_replication/run_wc_noise.R`, `runs/cv/wz_replication/run_wc_overlap.R`, `runs/cv/wz_replication/run_wc_overlap_n.R`, `runs/cv/wz_replication/run_wc_path.R` (kept in place for now):

```diff
@@ -1 +1 @@
-function(delta) tryCatch(fit_sbw(dat, colnames(X), delta, estimand, solver), error = function(e) conditionMessage(e))
+function(f, ...) f(X, W, target = "control", max.imbalance = 1e-04, maxit = 10000, tol = 1e-05, ...)
```

## `fmt`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/extension_cv/cv_summary.qmd`
- variant 2: `analysis/cv/basic_dgp/cv_summary.qmd`
- variant 3: `analysis/cv/wz_replication/wc_grids.R`

Diff of variant 1 vs variant 2 (3 variants in total):

```diff
@@ -1 +1,6 @@
-function(d) mutate(d, across(-1, function(v) sprintf("%.3f", v)))
+function(d, first) {
+    d <- d[order(d[[first]]), ]
+    out <- data.frame(d[[first]], sprintf("%.3f", d$lambda_opt), sprintf("%.3f", d$cv.bloss), sprintf("%.3f", d$boot.inf))
+    out[out == "NA"] <- ""
+    out
+}
```

## `gen_cell`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- variant 2: `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- variant 3: `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_machine_2.R`, `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- variant 4: `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R`, `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`, `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R`
- variant 5: `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- variant 6: `archive/exploring_cv/balnet_mini_experiment_windows.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- variant 7: `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- variant 8: `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- variant 9: `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- variant 10: `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`, `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- variant 11: `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`, `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- variant 12: `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`, `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`

Diff of variant 1 vs variant 2 (12 variants in total):

```diff
@@ -1,5 +1,9 @@
-function() {
-    dat <- gen_data(n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
+function(g) {
+    dat <- gen_data(g$n, p = p, overlap = "bad", s_y = 5, sigma_y = 1)
+    if (g$overlap == "awful") {
+        eta <- as.numeric(dat$X[, 1:5] %*% rep(4/sqrt(5), 5))
+        dat$W <- rbinom(g$n, 1, plogis(eta))
+    }
     mu <- as.numeric(dat$X[, 1:5] %*% rep(1/sqrt(5), 5))
     dat$Y <- mu + outer(dat$Y - mu, sigmas)
     dat
```

## `generate_data`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `tests/dgp2_validity.R`
- variant 2: `archive/tests/dgp2_validity_t1.R`

Diff of variant 1 vs variant 2:

```diff
@@ -3,9 +3,6 @@
     covcor <- match.arg(covcor)
     signs <- match.arg(signs)
     outcome_set <- match.arg(outcome_set)
-    .Sig <- function(idx, covcor) if (covcor == "iid") 
-        diag(length(idx))
-    else 0.5^abs(outer(idx, idx, "-"))
     if (is.null(idx_ps)) 
         idx_ps <- seq_len(s)
     else s <- length(idx_ps)
@@ -23,8 +20,7 @@
     alpha <- if (treat_prop == 0.5) 
         0
     else uniroot(function(a) integrate(function(z) plogis(a - z/overlap) * dnorm(z), -Inf, Inf)$value - treat_prop, c(-50, 50))$root
-    e <- plogis(alpha - eta)
-    e0 <- plogis(eta - alpha)
+    e <- 1/(1 + exp(eta - alpha))
     W <- rbinom(n, 1, e)
     Xo <- X[, idx_out, drop = FALSE]
     if (outcome_set == "fixed4") {
@@ -36,8 +32,7 @@
         a <- a/sqrt(drop(crossprod(a, .Sig(idx_out, covcor) %*% a)))
         f <- function(o) switch(o, linear = drop(Xo %*% a), quad1 = drop(pmax(Xo, 0)^2 %*% a) - 0.5 * sum(a), exp = drop(exp(Xo/2) %*% a) - exp(1/8) * sum(a))
     }
-    eps <- rnorm(n)
-    Y <- vapply(outcome, function(o) tau * W + f(o) + eps, numeric(n))
+    Y <- drop(tau * W + sapply(outcome, f) + rnorm(n))
     if (misspec) {
         j <- idx_ps[1:4]
         X[, j[1]] <- exp(0.5 * X_true[, j[1]])
@@ -46,5 +41,5 @@
         X[, j[4]] <- (X_true[, j[2]] + X_true[, j[4]] + 20)^2
         X[, j] <- scale(X[, j])
     }
-    list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta, e0 = e0)
+    list(Y = Y, W = W, X = X, X_true = X_true, tau = tau, e = e, eta = eta)
 }
```

## `grid_png`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `analysis/cv/basic_dgp/plot_grids.R`
- variant 2: `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R`

Diff of variant 1 vs variant 2:

```diff
@@ -1,5 +1,5 @@
 function(family, n_rep, rows, cols, out, keep = function(cell) TRUE, title = NULL) {
-    files <- list.files(dir, sprintf("^%s_?\\d+_r%d\\.rds$", family, n_rep), full.names = TRUE)
+    files <- list.files(res_dir, sprintf("^%s_?\\d+_r%d\\.rds$", family, n_rep), full.names = TRUE)
     if (length(files) == 0) 
         stop("no files for family ", family)
     panels <- unlist(lapply(files, function(f) {
@@ -22,7 +22,7 @@
     use <- use[do.call(order, cells[use, c(cols, rows), drop = FALSE])]
     nr <- nrow(unique(cells[use, rows, drop = FALSE]))
     nc <- length(use)/nr
-    png(file.path(dir, out), width = 480 * nc, height = 360 * nr + 60, res = 150)
+    png(file.path(fig_dir, out), width = 480 * nc, height = 360 * nr + 60, res = 150)
     par(mfcol = c(nr, nc), mar = c(3.5, 3.5, 2.5, 0.5), mgp = c(2.2, 0.7, 0), oma = c(0, 0, if (is.null(title)) 0 else 2, 0), cex = 0.7)
     for (i in use) {
         r <- panels[[i]]$rmse
```

## `num`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `analysis/ipw/correctspec.qmd`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1 @@
+function(z) as.numeric(unlist(z, use.names = FALSE))
```

## `one_rep`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/cv/wz_replication/run_wc_path.R`
- variant 2: `runs/cv/wz_replication/run_wc_overlap_n.R`
- variant 3: `runs/cv/wz_replication/run_wc_overlap.R`
- variant 4: `runs/cv/wz_replication/run_wc_basis.R`
- variant 5: `runs/cv/wz_replication/run_wc_noise.R`
- variant 6: `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_machine_2.R`, `archive/exploring_cv/balnet_mini_experiment_windows.R`, `tests/cv_smoke_test.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- variant 7: `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`, `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- variant 8: `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- variant 9: `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`, `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- variant 10: `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- variant 11: `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- variant 12: `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`, `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`, `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`, `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R`
- variant 13: `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`, `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R`
- variant 14: `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- variant 15: `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- variant 16: `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`, `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`

Diff of variant 1 vs variant 2 (16 variants in total):

```diff
@@ -2,22 +2,26 @@
     assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
     t0 <- proc.time()[["elapsed"]]
     out <- tryCatch({
-        d <- dgp_wc(n)
+        d <- dgp_wc_overlap(cell$n, overlap = cell$c)
         W <- d$W
         X <- cbind(d$X, d$X^2)
         fit <- function(f, ...) f(X, W, target = "control", max.imbalance = 1e-04, maxit = 10000, tol = 1e-05, ...)
         bl <- fit(balnet)
         sel <- list(cv.bloss = fit(cv.balnet, nfolds = 5, type.measure = "balance.loss"), cv.smd = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.mean"), cv.inf = fit(cv.balnet, nfolds = 5, type.measure = "imbalance.inf"), boot.smd = fit(cv.boot.balnet, type.measure = "imbalance.mean"), boot.inf = fit(cv.boot.balnet, type.measure = "imbalance.inf"))
         ybar1 <- colMeans(d$Y[W == 1, , drop = FALSE])
-        att <- function(w) {
-            g <- (w - 1) * (1 - W)
-            drop(ybar1 - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
-        }
+        att_w <- function(g) drop(ybar1 - sweep(crossprod(d$Y, g), 2, colSums(g), "/"))
+        att <- function(w) att_w((w - 1) * (1 - W))
+        att_e <- function(e) att_w(cbind((1 - W) * e/(1 - e)))
         wg <- balweights(bl, lambda = grid)
         s <- apply(X, 2, sd)
         cs <- vapply(seq_along(grid), function(k) cstat((wg[, k] - 1) * (1 - W), X, W, arms = 0, target = colMeans(X[W == 1, ]), s), numeric(1))
         k <- which.min(cs)
-        list(est_path = att(balweights(bl, lambda = lam)), lam_end = min(bl$lambda), est_sel = cbind(vapply(sel, function(m) att(balweights(m)), numeric(2)), alg1 = att(wg[, k, drop = FALSE])), lam_sel = c(vapply(sel, `[[`, numeric(1), "lambda.min"), alg1 = grid[k]), satt = mean(d$tau_i[W == 1]))
+        cvg <- cv.glmnet(X, W, family = "binomial", nfolds = 5)
+        e_cv <- drop(predict(cvg, X, s = "lambda.min", type = "response"))
+        e_ml <- fitted(glm(W ~ X, family = binomial))
+        clip <- function(e) pmin(e, 1 - 1e-06)
+        list(est_path = att(balweights(bl, lambda = lam)), lam_end = min(bl$lambda), est_sel = cbind(vapply(sel, function(m) att(balweights(m)), numeric(2)), alg1 = att(wg[, k, drop = FALSE])), lam_sel = c(vapply(sel, `[[`, numeric(1), "lambda.min"), alg1 = grid[k]), est_mle = cbind(glmnet = att_e(clip(e_cv)), glm = att_e(clip(e_ml)), naive = ybar1 - colMeans(d$Y[W == 0, , drop = FALSE])), lam_glmnet = cvg$lambda.min, clipped = c(glmnet = mean(e_cv[W == 0] > 1 - 1e-06), glm = mean(e_ml[W == 0] > 
+            1 - 1e-06)), satt = mean(d$tau_i[W == 1]))
     }, error = function(e) list(err = conditionMessage(e)))
     out$time_sec <- proc.time()[["elapsed"]] - t0
     out
```

## `panel`

Canonical (kept): `R/plots.R`

Variant in `analysis/cv/wz_replication/wc_grids.R` (kept in place for now):

```diff
@@ -1,14 +1,6 @@
-function(data, oc, x = n, ests = est_main, extra = NULL, title = oc_labs[oc]) {
-    d <- filter(data, outcome == oc)
-    l <- line_plot(d, rmse, "RMSE", x = {
-        {
-            x
-        }
-    }, ests = ests) + extra
-    r <- line_plot(d, abs(bias), "|Bias|", x = {
-        {
-            x
-        }
-    }, ests = ests) + extra
-    (l | r) + plot_layout(guides = "collect") + plot_annotation(title = title) & theme(legend.position = "bottom")
+function(z, o, main) {
+    m <- z[[o]]
+    plot(z$lam, m$r, log = "x", xlim = rev(range(z$lam)), type = "l", lwd = 1.5, ylim = range(c(m$r, m$rmse_sel), na.rm = TRUE), xlab = "lambda", ylab = "RMSE", main = main, cex.main = 0.95)
+    abline(v = z$lam[m$opt], lty = 3, col = "gray40")
+    abline(h = m$rmse_sel, col = cols[names(m$rmse_sel)], lty = ltys[names(m$rmse_sel)], lwd = 1.5)
 }
```

## `panel_rmse`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/cv/basic_dgp/make_cv_picks.R`
- variant 2: `analysis/cv/basic_dgp/plot_grids.R`, `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R`

Diff of variant 1 vs variant 2:

```diff
@@ -1,6 +1,17 @@
-function(x) {
-    res <- x$res[vapply(x$res, function(r) is.null(r$err), logical(1))]
-    path <- do.call(rbind, lapply(res, `[[`, "est_path"))
-    sel <- do.call(rbind, lapply(res, `[[`, "est_sel"))
-    list(lam = x$lam, path = sqrt(colMeans(path^2)), sel = sqrt(colMeans(sel^2)))
+function(x, s = NULL) {
+    ok <- vapply(x$res, function(r) is.null(r$err), logical(1))
+    res <- x$res[ok]
+    reached <- if (is.null(res[[1]]$lam_end)) 
+        1
+    else colMeans(outer(vapply(res, `[[`, numeric(1), "lam_end"), x$lam, "<="))
+    keep <- reached >= 0.95
+    if (is.null(s)) {
+        ep <- do.call(rbind, lapply(res, `[[`, "est_path"))
+        es <- do.call(rbind, lapply(res, `[[`, "est_sel"))
+    }
+    else {
+        ep <- do.call(rbind, lapply(res, function(r) r$est_path[, s]))
+        es <- do.call(rbind, lapply(res, function(r) r$est_sel[s, ]))
+    }
+    list(lam = x$lam[keep], path = sqrt(colMeans(ep^2))[keep], sel = sqrt(colMeans(es^2))[sels])
 }
```

## `patch_cv_balnet`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,68 @@
+function() {
+    target <- function(e, head) is.call(e) && identical(e[[1]], as.name(head)) && length(e) == 2L && is.call(e[[2]]) && identical(e[[2]][[1]], quote(matrix)) && length(e[[2]]) == 4L && identical(e[[2]][[if (head == "colMeans") 
+        3L
+    else 4L]], quote(nfolds))
+    walk <- function(e, fun) {
+        if (is.call(e)) {
+            fun(e)
+            for (i in seq_along(e)) {
+                if (is.null(e[[i]])) 
+                  next
+                if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) 
+                  next
+                walk(e[[i]], fun)
+            }
+        }
+    }
+    scan <- function(e, head) {
+        n <- 0L
+        walk(e, function(x) if (target(x, head)) 
+            n <<- n + 1L)
+        n
+    }
+    arity_ok <- function(e) {
+        ok <- TRUE
+        walk(e, function(x) if (is.name(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=", "<<-") && length(x) != 3L) 
+            ok <<- FALSE)
+        ok
+    }
+    rewrite <- function(e) {
+        if (is.call(e)) {
+            if (target(e, "colMeans")) {
+                e[[1]] <- quote(rowMeans)
+                m <- e[[2]]
+                tmp <- m[[3]]
+                m[[3]] <- m[[4]]
+                m[[4]] <- tmp
+                e[[2]] <- m
+                return(e)
+            }
+            for (i in seq_along(e)) {
+                if (is.null(e[[i]])) 
+                  next
+                if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) 
+                  next
+                e[[i]] <- rewrite(e[[i]])
+            }
+        }
+        e
+    }
+    f <- getFromNamespace("cv.balnet", "balnet")
+    b <- body(f)
+    if (!arity_ok(b)) 
+        stop("cv.balnet body is corrupted from an earlier patch attempt: restart R, then rerun.")
+    if (scan(b, "rowMeans") == 2L) 
+        return(invisible(TRUE))
+    if (scan(b, "colMeans") != 2L) 
+        stop("cv.balnet body unexpected (version change?): nothing modified.")
+    body(f) <- rewrite(b)
+    stopifnot(scan(body(f), "colMeans") == 0L, scan(body(f), "rowMeans") == 2L, arity_ok(body(f)))
+    assignInNamespace("cv.balnet", f, ns = "balnet")
+    if ("package:balnet" %in% search()) {
+        pe <- as.environment("package:balnet")
+        unlockBinding("cv.balnet", pe)
+        assign("cv.balnet", f, envir = pe)
+        lockBinding("cv.balnet", pe)
+    }
+    invisible(TRUE)
+}
```

## `rewrite`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,21 @@
+function(e) {
+    if (is.call(e)) {
+        if (target(e, "colMeans")) {
+            e[[1]] <- quote(rowMeans)
+            m <- e[[2]]
+            tmp <- m[[3]]
+            m[[3]] <- m[[4]]
+            m[[4]] <- tmp
+            e[[2]] <- m
+            return(e)
+        }
+        for (i in seq_along(e)) {
+            if (is.null(e[[i]])) 
+                next
+            if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) 
+                next
+            e[[i]] <- rewrite(e[[i]])
+        }
+    }
+    e
+}
```

## `rmse`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `analysis/cv/wz_replication/wc_grids.R`, `analysis/cv/wz_replication/wc_findings.qmd`
- variant 2: `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_machine_2.R`, `archive/exploring_cv/balnet_mini_experiment_windows.R`, `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`, `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`, `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`, `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`, `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`, `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`, `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`, `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`, `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`, `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`, `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- variant 3: `analysis/cv/wz_replication/plot_wc_att.R`

Diff of variant 1 vs variant 2 (3 variants in total):

```diff
@@ -1 +1 @@
-function(e) sqrt(colMeans(e^2))
+function(m) sqrt(colMeans(m^2))
```

## `scan`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,6 @@
+function(e, head) {
+    n <- 0L
+    walk(e, function(x) if (target(x, head)) 
+        n <<- n + 1L)
+    n
+}
```

## `selftest_cv_fix`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,21 @@
+function() {
+    set.seed(1)
+    n <- 400
+    p <- 10
+    X <- matrix(rnorm(n * p), n, p)
+    W <- rbinom(n, 1, plogis(X[, 1]))
+    nfolds <- 5
+    foldid <- sample(rep(seq_len(nfolds), length.out = n))
+    cvfit <- cv.balnet(X, W, nfolds = nfolds, foldid = foldid, max.imbalance = 0.01)
+    lam <- cvfit$`_lambda`
+    v <- unlist(lapply(seq_len(nfolds), function(k) {
+        tr <- foldid != k
+        f <- balnet(X[tr, ], W[tr], standardize = ".inplace", max.imbalance = 0.01)
+        balnet:::get_balance_loss(f, X[!tr, , drop = FALSE], W[!tr], rep(1, sum(!tr)), lam)$treated
+    }))
+    L <- length(lam$treated)
+    fix <- rowMeans(matrix(v, L, nfolds))
+    num <- function(z) as.numeric(unlist(z, use.names = FALSE))
+    stopifnot(isTRUE(all.equal(num(cvfit$`_cv.info`$cv.mean$treated), fix)), isTRUE(all.equal(num(cvfit$lambda.min$treated), lam$treated[which.min(fix)])))
+    message("cv.balnet patch verified")
+}
```

## `summarise_file`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/extension_cv/scripts/summarise_runs_to_csv.R`
- variant 2: `analysis/cv/basic_dgp/prep_cv_summary.R`

Diff of variant 1 vs variant 2:

```diff
@@ -27,12 +27,21 @@
         j_min <- which.min(rp)
         j_floor <- max(which(keep))
         rs <- sqrt(colMeans(es^2))[sels]
+        lm <- colMeans(do.call(rbind, lapply(res, function(r) {
+            v <- r$lam_sel
+            if (is.null(v)) 
+                return(c(cv.bloss = NA_real_, boot.inf = NA_real_))
+            if (!is.null(dim(v))) 
+                v <- v[s, ]
+            v[c("cv.bloss", "boot.inf")]
+        })))
         best <- sels[which.min(rs)]
         z_vs_floor <- function(k) {
             d <- es[, k]^2 - ep[, j_floor]^2
             mean(d)/(sd(d)/sqrt(length(d)))
         }
-        data.frame(file = basename(f), cell, sigma_y = sig[s], n_ok = sum(ok), lam_min = x$lam[j_min], lam_floor = x$lam[j_floor], rmse_min = rp[j_min], rmse_floor = rp[j_floor], gain_floor = rp[j_floor]/rp[j_min] - 1, as.list(setNames(rs/rp[j_min], paste0(sels, "_vs_min"))), best = best, rmse_best = rs[best], best_vs_floor = rs[best]/rp[j_floor] - 1, z_best = z_vs_floor(best), z_bootinf = z_vs_floor("boot.inf"), z_cvbloss = z_vs_floor("cv.bloss"))
+        data.frame(file = basename(f), cell, sigma_y = sig[s], n_ok = sum(ok), lam_min = x$lam[j_min], lam_floor = x$lam[j_floor], rmse_min = rp[j_min], rmse_floor = rp[j_floor], gain_floor = rp[j_floor]/rp[j_min] - 1, as.list(setNames(rs/rp[j_min], paste0(sels, "_vs_min"))), lam_bloss = unname(lm["cv.bloss"]), lam_bootinf = unname(lm["boot.inf"]), RMSE_cv = unname(rs["cv.bloss"]), RMSE_boot = unname(rs["boot.inf"]), CV_excess_pct = 100 * (unname(rs["cv.bloss"])/rp[j_min] - 1), Boot_excess_pct = 100 * 
+            (unname(rs["boot.inf"])/rp[j_min] - 1), best = best, rmse_best = rs[best], best_vs_floor = rs[best]/rp[j_floor] - 1, z_best = z_vs_floor(best), z_bootinf = z_vs_floor("boot.inf"), z_cvbloss = z_vs_floor("cv.bloss"))
     })
     do.call(rbind, by_sigma)
 }
```

## `target`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1 @@
+function(e, head) is.call(e) && identical(e[[1]], as.name(head)) && length(e) == 2L && is.call(e[[2]]) && identical(e[[2]][[1]], quote(matrix)) && length(e[[2]]) == 4L && identical(e[[2]][[if (head == "colMeans") 3L else 4L]], quote(nfolds))
```

## `title_of`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- variant 2: `archive/exploring_cv/balnet_mini_experiment_machine_2.R`
- variant 3: `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- variant 4: `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- variant 5: `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- variant 6: `archive/exploring_cv/balnet_mini_experiment_edited.R`, `archive/exploring_cv/balnet_mini_experiment_windows.R`, `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- variant 7: `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- variant 8: `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- variant 9: `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`

Diff of variant 1 vs variant 2 (9 variants in total):

```diff
@@ -1 +1 @@
-function(cell) sprintf("%s overlap, p = %d, %s", cell$overlap, cell$p, if (cell$dense) "dense" else "sparse")
+function(cell) sprintf("%s overlap, s_y = %d, sigma_y = %d", cell$overlap, cell$s_y, cell$sigma_y)
```

## `walk`

No `R/` version. Script-local helper; variants kept in place:

- variant 1: `archive/tests/balnet_fix.R`
- variant 2: `runs/ipw/cvfix.R`, `tests/balnet_cvfix_gate.R`

Diff of variant 1 vs variant 2:

```diff
@@ -0,0 +1,12 @@
+function(e, fun) {
+    if (is.call(e)) {
+        fun(e)
+        for (i in seq_along(e)) {
+            if (is.null(e[[i]])) 
+                next
+            if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) 
+                next
+            walk(e[[i]], fun)
+        }
+    }
+}
```

