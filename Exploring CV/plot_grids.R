# Legible RMSE-path grids, base graphics as in Erik's script ----
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
ov_lev <- c("good", "moderate", "bad", "awful")
sels <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
sel_col <- c("dodgerblue", "dodgerblue", "dodgerblue", "red", "red")
sel_lty <- c(4, 1, 2, 1, 2)

#' RMSE path (to the 95 percent floor) and selector RMSE for one file,
#' column s of est_path when the file holds several sigma_y.
panel_rmse <- function(x, s = NULL) {
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  res <- x$res[ok]
  reached <- if (is.null(res[[1]]$lam_end)) 1 else
    colMeans(outer(vapply(res, `[[`, numeric(1), "lam_end"), x$lam, "<="))
  keep <- reached >= 0.95
  if (is.null(s)) {
    ep <- do.call(rbind, lapply(res, `[[`, "est_path"))
    es <- do.call(rbind, lapply(res, `[[`, "est_sel"))
  } else {
    ep <- do.call(rbind, lapply(res, \(r) r$est_path[, s]))
    es <- do.call(rbind, lapply(res, \(r) r$est_sel[s, ]))
  }
  list(lam = x$lam[keep], path = sqrt(colMeans(ep^2))[keep],
       sel = sqrt(colMeans(es^2))[sels])
}

#' One PNG: panels of family filtered by keep(cell), rows keyed by `rows`,
#' columns keyed by `cols`. `title` is printed across the top.
grid_png <- function(family, n_rep, rows, cols, out, keep = \(cell) TRUE,
                     title = NULL) {
  files <- list.files(dir, sprintf("^%s_?\\d+_r%d\\.rds$", family, n_rep),
                      full.names = TRUE)
  if (length(files) == 0) stop("no files for family ", family)
  panels <- unlist(lapply(files, function(f) {
    x <- readRDS(f)
    sig <- if (!is.null(x$sigmas)) x$sigmas else
      if (!is.null(x$cell$sigma_y)) x$cell$sigma_y else 1
    lapply(seq_along(sig), function(s) {
      cell <- x$cell
      cell$sigma_y <- sig[s]
      list(cell = cell,
           rmse = panel_rmse(x, if (is.null(x$sigmas)) NULL else s))
    })
  }), recursive = FALSE)
  cells <- do.call(rbind, lapply(panels, `[[`, "cell"))
  if (!is.null(cells$overlap)) cells$overlap <- factor(cells$overlap, ov_lev)
  use <- which(vapply(seq_len(nrow(cells)), \(i) keep(cells[i, ]),
                      logical(1)))
  use <- use[do.call(order, cells[use, c(cols, rows), drop = FALSE])]
  nr <- nrow(unique(cells[use, rows, drop = FALSE]))
  nc <- length(use) / nr
  png(file.path(dir, out), width = 480 * nc, height = 360 * nr + 60,
      res = 150)
  par(mfcol = c(nr, nc), mar = c(3.5, 3.5, 2.5, 0.5), mgp = c(2.2, 0.7, 0),
      oma = c(0, 0, if (is.null(title)) 0 else 2, 0), cex = 0.7)
  for (i in use) {
    r <- panels[[i]]$rmse
    cell <- cells[i, c(rows, cols)]
    plot(r$lam, r$path, log = "x", xlim = rev(range(r$lam)), type = "l",
         xlab = "lambda", ylab = "RMSE", cex.main = 0.85,
         main = paste(names(cell), vapply(cell, format, character(1)),
                      sep = " = ", collapse = ", "))
    abline(v = r$lam[which.min(r$path)], lty = 3, col = "gray40")
    abline(h = r$sel, col = sel_col, lty = sel_lty, lwd = 1.5)
  }
  legend("topright", sels, col = sel_col, lty = sel_lty, lwd = 2,
         bg = "white", cex = 0.8)
  if (!is.null(title)) mtext(title, outer = TRUE, cex = 0.9, font = 2)
  dev.off()
}

# Figures, at most three columns each ----
file.remove(list.files(dir, "^fig_", full.names = TRUE))
ov_pair <- list(strong = c("good", "moderate"), weak = c("bad", "awful"))
sig3 <- c(1, 3, 10)
for (k in names(ov_pair)) {
  grid_png("tune4", 500, "overlap", "sigma_y",
           sprintf("fig_snr_lasso_%s.png", k),
           keep = \(c) c$overlap %in% ov_pair[[k]] & c$sigma_y %in% sig3,
           title = "Lasso")
  grid_png("snr_enet", 500, "overlap", "sigma_y",
           sprintf("fig_snr_enet_%s.png", k),
           keep = \(c) c$overlap %in% ov_pair[[k]] & c$sigma_y %in% sig3,
           title = "Elastic net, alpha 0.5")
  for (sg in sig3)
    grid_png("tunen2", 500, "overlap", "n",
             sprintf("fig_n_lasso_s%d_%s.png", sg, k),
             keep = \(c) c$sigma_y == sg & c$overlap %in% ov_pair[[k]] &
               c$n != 500,
             title = sprintf("Lasso, sigma_y = %d", sg))
}
for (ov in c("bad", "awful"))
  grid_png("n_enet", 200, "sigma_y", "n", sprintf("fig_n_enet_%s.png", ov),
           keep = \(c) c$overlap == ov,
           title = sprintf("Elastic net, alpha 0.5, %s overlap", ov))
grid_png("tune", 500, "overlap", "sigma_y", "fig_density_sy5.png",
         keep = \(c) c$s_y == 5, title = "Sparse outcome, s_y = 5")
grid_png("tune", 500, "overlap", "sigma_y", "fig_density_sy100.png",
         keep = \(c) c$s_y == 100,
         title = "Diluted confounding, 100 outcome covariates at 0.10")
grid_png("tune2", 500, "overlap", "sigma_y", "fig_density_fixed_sparse.png",
         keep = \(c) !c$dense, title = "Fixed confounding, sparse")
grid_png("tune2", 500, "overlap", "sigma_y", "fig_density_fixed_dense.png",
         keep = \(c) c$dense, title = "Fixed confounding, dense")
grid_png("alpha_bad", 200, "alpha", "sigma_y", "fig_alpha_bad.png",
         title = "Alpha sweep, bad overlap")
grid_png("tunea4", 200, "overlap", "alpha", "fig_ridge_floor.png",
         title = "Path extended to 1e-4, sigma_y = 1")
grid_png("dima", 200, "overlap", "p", "fig_p_low.png",
         keep = \(c) c$alpha == 0.5, title = "Elastic net, alpha 0.5")
for (pp in c(1000, 2000))
  grid_png("dimhi", 200, "overlap", "alpha", sprintf("fig_p_high_%d.png", pp),
           keep = \(c) c$alpha != 0.75 & c$p == pp,
           title = sprintf("p = %d, sigma_y = 1", pp))
for (ov in c("moderate", "bad")) for (pp in c(500, 1000))
  grid_png("dimsnr", 200, "sigma_y", "alpha",
           sprintf("fig_p_snr_%s_%d.png", ov, pp),
           keep = \(c) c$overlap == ov & c$p == pp,
           title = sprintf("p = %d, %s overlap", pp, ov))
for (a in c(1, 0.5))
  grid_png("spread", 200, "overlap", "s", sprintf("fig_spread_a%g.png", a),
           keep = \(c) c$alpha == a,
           title = if (a == 1) "Lasso" else "Elastic net, alpha 0.5")


# Added figures ----
for (k in names(ov_pair)) {
  grid_png("tune4", 500, "overlap", "sigma_y",
           sprintf("fig_snr_lasso_%s_mid.png", k),
           keep = \(c) c$overlap %in% ov_pair[[k]] & c$sigma_y %in% c(2, 5),
           title = "Lasso")
  grid_png("snr_enet", 500, "overlap", "sigma_y",
           sprintf("fig_snr_enet_%s_mid.png", k),
           keep = \(c) c$overlap %in% ov_pair[[k]] & c$sigma_y %in% c(2, 5),
           title = "Elastic net, alpha 0.5")
}
grid_png("tunea3", 500, "overlap", "alpha", "fig_ridge_default.png",
         keep = \(c) !c$dense, title = "Default path floor, sigma_y = 1")
grid_png("spread", 200, "overlap", "s", "fig_spread_a0.png",
         keep = \(c) c$alpha == 0, title = "Ridge, alpha 0")



# Density side by side, one sigma_y per figure ----
for (sg in c(1, 5)) {
  grid_png("tune", 500, "overlap", "s_y",
           sprintf("fig_density_sy_sigma%d.png", sg),
           keep = \(c) c$sigma_y == sg,
           title = sprintf("Sparse (s_y = 5) vs diluted (s_y = 100), sigma_y = %d", sg))
  grid_png("tune2", 500, "overlap", "dense",
           sprintf("fig_density_fixed_sigma%d.png", sg),
           keep = \(c) c$sigma_y == sg,
           title = sprintf("Fixed confounding, sparse vs dense, sigma_y = %d", sg))
}


