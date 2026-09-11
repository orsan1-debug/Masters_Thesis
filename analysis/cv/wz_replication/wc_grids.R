# wc_grids.R ---------------------------------------------------------------
# Functions for the path figures and tables on the Wong & Chan design in the
# layout of cv_summary.pdf (Sections 2.1, 2.3, 2.4; Figures 2 and 3).
# Sourced by plot_wc_grids.R (writes png to output/cv/wz_replication/figures,
# csv to output/cv/wz_replication/tables) and by the qmd (draws inline).
# read_grid() and write_tabs() take root = here::here() (Phase 3b); the batches
# live in results/cv/wz_replication/.
# Batches: results/wc_noise_s<sigma>_c<c>.rds (noise x overlap, n = 5000),
# results/wc_on_n<n>_c<c>.rds (n x overlap, sigma = 1) and
# results/wc_basis_K<K>_c<c>.rds (basis size x overlap, n = 5000, sigma = 1).
# The sigma = 1 row of the noise grid is the overlap axis (identical draws
# and estimates to results/wc_overlap_c<c>.rds).
# Per cell and outcome model: RMSE of the balnet ATT along lambda (reversed
# log axis, exact balance to the right; points reached by fewer than 95% of
# reps are dropped, so the right end of the curve is the path floor), a
# dotted line at the RMSE-optimal lambda, and a horizontal line per tuning
# rule (cv.balnet: balance loss, mean SMD, max SMD; cv.boot.balnet: mean
# SMD, max SMD; Wang & Zubizarreta Algorithm 1). Truth is the sample ATT
# for model A and 0 for model B. Then the median selected lambda per rule
# against the RMSE-optimal lambda along each factor axis, and tables of
# lambda and RMSE (MCSE, delta method; vs_floor = RMSE / path-floor RMSE)
# per cell. Colours as in plot_wc_att.R / the qmd figure: cv.balnet rules
# blue (balance loss dot-dash, mean SMD solid, max SMD dashed),
# cv.boot.balnet rules red (mean SMD solid, max SMD dashed), alg1 black.

rmse <- \(e) sqrt(colMeans(e^2))
mcse <- \(e) apply(e^2, 2, sd) / sqrt(nrow(e)) / (2 * rmse(e))
cols <- c(cv.bloss = "blue", cv.smd = "blue", cv.inf = "blue",
          boot.smd = "red", boot.inf = "red", alg1 = "black")
ltys <- c(cv.bloss = 4, cv.smd = 1, cv.inf = 2, boot.smd = 1, boot.inf = 2,
          alg1 = 1)

# One cell: RMSE along the path and per rule, both outcome models ----
read_cell <- function(file) {
  x   <- readRDS(file)
  ok  <- vapply(x$res, \(r) is.null(r$err), logical(1))
  res <- x$res[ok]
  lam_end <- vapply(res, `[[`, numeric(1), "lam_end")
  reached <- colMeans(outer(lam_end, x$lam * (1 + 1e-6), "<="))
  satt    <- vapply(res, `[[`, numeric(1), "satt")
  lam_sel <- do.call(rbind, lapply(res, `[[`, "lam_sel"))
  bind    <- \(o, k) do.call(rbind, lapply(res, \(r) r[[k]][o, ]))
  model   <- \(o) {
    truth  <- if (o == "A") satt else 0
    e_path <- bind(o, "est_path") - truth
    r <- rmse(e_path); r[reached < 0.95] <- NA
    floor <- max(which(!is.na(r))); opt <- which.min(r)
    e_sel <- bind(o, "est_sel") - truth
    list(r = r, floor = floor, opt = opt,
         rmse_path = r[c(opt, floor)], mcse_path = mcse(e_path)[c(opt, floor)],
         rmse_sel = rmse(e_sel), mcse_sel = mcse(e_sel))
  }
  list(cell = x$cell, n = x$n, n_rep = length(res), lam = x$lam,
       reach = mean(lam_end <= 1e-4 * (1 + 1e-6)),
       lam_med = apply(lam_sel, 2, median),
       A = model("A"), B = model("B"))
}

# Grid of cells: rows x columns of read_cell() results ----
read_grid <- function(prefix, rows, columns, root = here::here()) {
  lapply(rows, \(r) lapply(columns, \(k)
                           read_cell(file.path(root, "results", "cv", "wz_replication",
                                               sprintf("%s%s_c%s.rds", prefix, r, k)))))
}

# One panel, cv_summary.pdf style ----
panel <- function(z, o, main) {
  m <- z[[o]]
  plot(z$lam, m$r, log = "x", xlim = rev(range(z$lam)), type = "l",
       lwd = 1.5, ylim = range(c(m$r, m$rmse_sel), na.rm = TRUE),
       xlab = "lambda", ylab = "RMSE", main = main, cex.main = 0.95)
  abline(v = z$lam[m$opt], lty = 3, col = "gray40")
  abline(h = m$rmse_sel, col = cols[names(m$rmse_sel)],
         lty = ltys[names(m$rmse_sel)], lwd = 1.5)
}

# Grid figure: one panel per cell, legend strip below. file = NULL draws on
# the current device (knitr); a path writes a png ----
grid_fig <- function(cells, rows, columns, row_lab, col_lab, o,
                     title = NULL, file = NULL) {
  nr <- length(rows); nc <- length(columns)
  if (!is.null(file)) {
    png(file, width = max(400 * nc, 1100), height = 340 * nr + 80, res = 105)
  }
  layout(rbind(matrix(seq_len(nr * nc), nr, nc, byrow = TRUE),
               rep(nr * nc + 1, nc)),
         heights = c(rep(1, nr), 0.2))
  par(mar = c(4, 4, 2.5, 1), oma = c(0, 0, if (is.null(title)) 0 else 3, 0))
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    panel(cells[[i]][[j]], o,
          sprintf("%s = %s, %s = %s", row_lab, rows[i], col_lab, columns[j]))
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", names(cols), col = cols, lty = ltys, lwd = 1.5,
         horiz = TRUE, bty = "n")
  if (!is.null(title)) mtext(title, outer = TRUE, cex = 1.3, font = 2)
  if (!is.null(file)) dev.off()
}

# Median selected lambda per rule against the RMSE-optimal lambda (Fig. 3) ----
# axis: values of the factor on the x axis; cells: list of read_cell() along it.
pick_panel <- function(axis, cells, xlab, main, log_x = TRUE) {
  opt <- sapply(c("A", "B"), \(o) sapply(cells, \(z) z$lam[z[[o]]$opt]))
  flo <- sapply(cells, \(z) z$lam[z$A$floor])            # same for A and B
  med <- sapply(cells, `[[`, "lam_med")                   # rules x axis
  plot(axis, opt[, "A"], log = if (log_x) "xy" else "y", type = "b", pch = 19,
       lwd = 2, col = "gray30", ylim = range(c(opt, med, flo)), xlab = xlab,
       ylab = "lambda", main = main, cex.main = 0.95)
  lines(axis, opt[, "B"], type = "b", pch = 1, lwd = 2, lty = 2, col = "gray30")
  lines(axis, flo, lty = 3, col = "gray50")
  for (k in rownames(med)) {
    lines(axis, med[k, ], type = "b", pch = 20, col = cols[k], lty = ltys[k])
  }
}

pick_fig <- function(cells, rows, columns, row_lab, col_lab, title = NULL,
                     file = NULL, log_rows = TRUE, log_cols = FALSE) {
  nr <- length(rows); nc <- length(columns); w <- max(nr, nc)
  if (!is.null(file)) png(file, width = max(380 * w, 1100), height = 760,
                          res = 105)
  layout(rbind(c(seq_len(nc), rep(0, w - nc)),         # 0 = empty slot
               c(nc + seq_len(nr), rep(0, w - nr)),
               rep(nc + nr + 1, w)),
         heights = c(1, 1, 0.22))
  par(mar = c(4, 4, 2.5, 1), oma = c(0, 0, if (is.null(title)) 0 else 3, 0))
  for (j in seq_len(nc)) {                        # row factor on the x axis
    pick_panel(rows, lapply(cells, `[[`, j), row_lab,
               sprintf("%s = %s", col_lab, columns[j]), log_rows)
  }
  for (i in seq_len(nr)) {                        # column factor on the x axis
    pick_panel(columns, cells[[i]], col_lab,
               sprintf("%s = %s", row_lab, rows[i]), log_cols)
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", c("RMSE-optimal, A", "RMSE-optimal, B", "path floor",
                     names(cols)),
         col = c("gray30", "gray30", "gray50", cols), lty = c(1, 2, 3, ltys),
         pch = c(19, 1, NA, rep(20, 6)), lwd = c(2, 2, 1, rep(1, 6)),
         horiz = TRUE, bty = "n", cex = 0.9)
  if (!is.null(title)) mtext(title, outer = TRUE, cex = 1.3, font = 2)
  if (!is.null(file)) dev.off()
}

# Tables: lambda per cell; RMSE (MCSE) per cell and model ----
lambda_tab <- function(cells, rows, columns, row_lab, col_lab) {
  do.call(rbind, lapply(seq_along(rows), \(i) do.call(rbind, lapply(
    seq_along(columns), \(j) {
      z <- cells[[i]][[j]]
      data.frame(rows[i], columns[j], n_rep = z$n_rep, reach = z$reach,
                 floor = z$lam[z$A$floor], opt_A = z$lam[z$A$opt],
                 opt_B = z$lam[z$B$opt], t(z$lam_med),
                 check.names = FALSE) |>
        setNames(c(row_lab, col_lab, "n_rep", "reach_1e-4", "floor",
                   "opt_A", "opt_B", paste0("med_", names(z$lam_med))))
    }))))
}

rmse_tab <- function(cells, rows, columns, row_lab, col_lab, o) {
  fmt <- \(r, m) sprintf("%.3f (%.3f)", r, m)
  do.call(rbind, lapply(seq_along(rows), \(i) do.call(rbind, lapply(
    seq_along(columns), \(j) {
      m <- cells[[i]][[j]][[o]]
      data.frame(rows[i], columns[j],
                 t(fmt(c(m$rmse_path, m$rmse_sel), c(m$mcse_path, m$mcse_sel))),
                 t(round(c(m$rmse_path[1], m$rmse_sel) / m$rmse_path[2], 2)),
                 check.names = FALSE) |>
        setNames(c(row_lab, col_lab, "oracle", "floor", names(m$rmse_sel),
                   paste0("vs_floor_", c("oracle", names(m$rmse_sel)))))
    }))))
}

write_tabs <- function(cells, rows, columns, row_lab, col_lab, stem,
                       root = here::here()) {
  write.csv(lambda_tab(cells, rows, columns, row_lab, col_lab),
            file.path(root, "output", "cv", "wz_replication", "tables",
                      paste0(stem, "_lambda.csv")), row.names = FALSE)
  for (o in c("A", "B")) {
    write.csv(rmse_tab(cells, rows, columns, row_lab, col_lab, o),
              file.path(root, "output", "cv", "wz_replication", "tables",
                        paste0(stem, "_rmse_", o, ".csv")),
              row.names = FALSE)
  }
}