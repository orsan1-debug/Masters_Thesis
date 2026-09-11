# plot_wc_att.R ------------------------------------------------------------
# Path plots in the style of Erik's tuning experiments, Wong & Chan design:
# RMSE of the balnet ATT along lambda (reversed log axis, exact balance to
# the right), vertical lines at the RMSE-optimal lambda (grey) and the
# median path endpoint (red), horizontal lines for each tuning rule, and
# grey reference lines for sbw exact / approx from wc_sbw_v2 (same draws,
# restricted to the reps present here). Truth is the sample ATT for model A
# and 0 for model B. Points where fewer than 95% of paths reached lambda
# are dropped. Then: selected-lambda quantiles per rule, and each rule's
# RMSE split into level (path RMSE at its selected lambda, averaged over
# reps) and spread (what selecting a different lambda each rep adds).

x   <- readRDS(here::here("results", "cv", "wz_replication", "wc_att_v2.rds"))
sbw <- data.table::fread(file = here::here("results", "cv", "wz_replication", "wc_sbw_v2.csv.gz"))

ok  <- vapply(x$res, \(r) is.null(r$err), logical(1))
res <- x$res[ok]
if (!all(ok)) message(sum(!ok), " failed reps dropped")
sbw <- subset(sbw, rep %in% which(ok))
lam_end <- vapply(res, `[[`, numeric(1), "lam_end")
reached <- colMeans(outer(lam_end, x$lam, "<="))     # share of reps at lambda
satt    <- vapply(res, `[[`, numeric(1), "satt")
lam_sel <- do.call(rbind, lapply(res, `[[`, "lam_sel"))

bind  <- \(o, k) do.call(rbind, lapply(res, \(r) r[[k]][o, ]))
truth <- \(o) if (o == "A") satt else 0
rmse  <- \(m, tau) sqrt(colMeans((m - tau)^2))
sbw_rmse <- \(o, b) {
  z <- subset(sbw, estimand == "att" & outcome == o & balance == b)
  sqrt(mean((z$tau_hat - if (o == "A") z$satt else 0)^2))
}

cols <- c(cv.bloss = "blue", cv.smd = "blue", cv.inf = "blue",
          boot.smd = "red", boot.inf = "red", alg1 = "black")
ltys <- c(4, 1, 2, 1, 2, 1)

# Figure and RMSE tables ----
png(here::here("output", "cv", "wz_replication", "figures", "wc_att_v2.png"), 1800, 800, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (o in c("A", "B")) {
  rmse_path <- rmse(bind(o, "est_path"), truth(o))
  rmse_path[reached < 0.95] <- NA
  rmse_sel  <- rmse(bind(o, "est_sel"), truth(o))
  sbw_ref   <- c(exact = sbw_rmse(o, "exact"), approx = sbw_rmse(o, "approx"))
  plot(x$lam, rmse_path, log = "x", xlim = rev(range(x$lam)), type = "l",
       lwd = 2, ylim = range(c(rmse_path, rmse_sel, sbw_ref), na.rm = TRUE),
       xlab = "lambda (control-arm max SMD, log scale)",
       ylab = "RMSE vs sample ATT",
       main = sprintf("ATT, outcome model %s (n = %d, %d reps)",
                      o, x$n, sum(ok)))
  abline(v = x$lam[which.min(rmse_path)], lty = 3, col = "gray40")
  abline(v = median(lam_end), lty = 3, col = "red")
  abline(h = rmse_sel, col = cols, lty = ltys)
  abline(h = sbw_ref, col = "gray60", lty = c(2, 1), lwd = 2)
  cat("\nOutcome", o, "\n")
  print(round(c(path_min = min(rmse_path, na.rm = TRUE),
                lam_opt  = x$lam[which.min(rmse_path)],
                path_end = rmse_path[max(which(!is.na(rmse_path)))],
                rmse_sel, sbw_exact = sbw_ref[["exact"]],
                sbw_approx = sbw_ref[["approx"]]), 3))
}
legend("topright", c(names(cols), "sbw exact", "sbw approx"),
       col = c(cols, "gray60", "gray60"), lty = c(ltys, 2, 1), lwd = 2,
       bg = "white")
dev.off()

# Selected lambda per rule ----
cat("\nSelected lambda\n")
print(signif(rbind(median = apply(lam_sel, 2, median),
                   q10 = apply(lam_sel, 2, quantile, 0.1),
                   q90 = apply(lam_sel, 2, quantile, 0.9)), 3))
cat("median path endpoint:", median(lam_end), "\n")

# Level / spread decomposition ----
for (o in c("A", "B")) {
  r <- rmse(bind(o, "est_path"), truth(o)); r[reached < 0.95] <- NA
  level <- apply(lam_sel, 2, \(l)
                 mean(approx(log(x$lam[!is.na(r)]), r[!is.na(r)], xout = log(l),
                             rule = 2)$y))
  rule <- rmse(bind(o, "est_sel"), truth(o))
  cat("\nOutcome", o, "\n")
  print(round(rbind(rule = rule, level = level, spread = rule - level), 3))
}