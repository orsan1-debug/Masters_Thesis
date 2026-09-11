# analyse_path_end_slope.R: relative RMSE of the last four path points against the floor,
# lasso, weak-overlap cells of tune4 (and tunen2 if present). Console only.
# From balnet_mini_experiment_WINDOWS.R lines 1045-1067, the only block in
# that file not also in balnet_mini_experiment_edited.R.
# Note: the last line auto-prints only when run interactively; wrap it in
# print() if you source() this file.

rm(list = ls())
# *** Setup ***
res_dir <- here::here("results", "cv", "basic_dgp")
stopifnot(dir.exists(res_dir))

# End slope of the lasso RMSE path, weak-overlap cells ----
files <- list.files(res_dir, "^(tune4|tunen2)_\\d+_r\\d+\\.rds$",
                    full.names = TRUE)
end_slope <- function(f) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  res <- x$res[ok]
  cell <- x$cell
  sig <- if (!is.null(x$sigmas)) x$sigmas else cell$sigma_y
  cell$sigma_y <- NULL
  k <- length(x$lam)
  rows <- lapply(seq_along(sig), function(s) {
    ep <- if (is.null(x$sigmas)) do.call(rbind, lapply(res, `[[`, "est_path"))
    else do.call(rbind, lapply(res, \(r) r$est_path[, s]))
    rp <- sqrt(colMeans(ep^2))
    data.frame(file = basename(f), cell, sigma_y = sig[s],
               lam_floor = x$lam[k],
               rel = t(round(rp[(k - 3):k] / rp[k], 3)))
  })
  do.call(rbind, rows)
}
slope <- dplyr::bind_rows(lapply(files, end_slope))
slope[slope$overlap %in% c("bad", "awful"), ]
