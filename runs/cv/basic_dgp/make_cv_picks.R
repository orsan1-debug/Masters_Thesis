

dir     <- here::here("results", "cv", "basic_dgp")
sum_dir <- here::here("output", "cv", "basic_dgp", "summaries")

# panel_rmse: RMSE over reps along the path and at each selector's pick ----
panel_rmse <- function(x) {
  res <- x$res[vapply(x$res, \(r) is.null(r$err), logical(1))]
  path <- do.call(rbind, lapply(res, `[[`, "est_path"))
  sel <- do.call(rbind, lapply(res, `[[`, "est_sel"))
  list(lam = x$lam,
       path = sqrt(colMeans(path^2)),       # true = 0
       sel = sqrt(colMeans(sel^2)))
}

chk <- panel_rmse(readRDS(file.path(dir, "tune4_01_r500.rds")))
c(rmse_min = min(chk$path), lam_opt = chk$lam[which.min(chk$path)],
  bloss_vs_min = chk$sel[["cv.bloss"]] / min(chk$path))



# Rebuild cv_picks with all five selectors ----
sel_names <- c("cv.bloss", "cv.smd", "boot.smd", "cv.inf", "boot.inf")
c_of <- c(good = 0.7, moderate = 1.5, bad = 2.5, awful = 4)
files <- list.files(dir,
                    "^(tune4|snr_good|snr_mb|ov_s1|ov_1k)_\\d+_r\\d+\\.rds$",
                    full.names = TRUE)
hand <- do.call(rbind, lapply(files, \(f) {
  x <- readRDS(f)
  r <- panel_rmse(x)
  res <- x$res[vapply(x$res, \(r) is.null(r$err), logical(1))]
  ls <- do.call(rbind, lapply(res, `[[`, "lam_sel"))
  cell <- x$cell
  if (is.null(cell$c_prop)) cell$c_prop <- c_of[[cell$overlap]]
  picks <- as.list(apply(ls[, sel_names, drop = FALSE], 2, median))
  rmses <- as.list(r$sel[sel_names] / min(r$path))
  names(rmses) <- paste0(sel_names, "_vs_min")
  data.frame(family = sub("_?\\d+_r\\d+\\.rds$", "", basename(f)),
             c_prop = cell$c_prop, sigma_y = cell$sigma_y,
             reps = length(res),
             lam_floor = min(r$lam), lambda_opt = r$lam[which.min(r$path)],
             floor_vs_min = r$path[which.min(r$lam)] / min(r$path),
             picks, rmses, check.names = FALSE)
}))
hand <- hand[order(hand$c_prop, hand$sigma_y, -hand$reps), ]
hand <- hand[!duplicated(hand[c("c_prop", "sigma_y")]), ]   # largest reps
readr::write_csv(hand, file.path(sum_dir, "cv_picks.csv"))

noise_cols <- c("sigma_y", "lambda_opt", sel_names)
ov_cols <- c("c_prop", "lambda_opt", sel_names)
readr::write_csv(hand[hand$c_prop == 0.7, noise_cols],
                 file.path(sum_dir, "cv_picks_noise.csv"))
readr::write_csv(hand[hand$sigma_y == 1, ov_cols],
                 file.path(sum_dir, "cv_picks_overlap.csv"))
print(hand[hand$c_prop == 0.7, noise_cols], digits = 3, row.names = FALSE)
print(hand[hand$sigma_y == 1, ov_cols], digits = 3, row.names = FALSE)

